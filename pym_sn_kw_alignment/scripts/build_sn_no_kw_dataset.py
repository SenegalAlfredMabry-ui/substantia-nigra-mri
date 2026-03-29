#!/usr/bin/env python3
"""Build feature/label tables for no-KW SN correction modeling.

This script extracts subject-level geometry features from candidate SN masks and
anchor ROIs. If a manual mask directory is provided, it also computes target
delta vectors (manual COM - candidate COM) for supervised training.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import nibabel as nib
import numpy as np


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--subjects-file",
        type=Path,
        help="Optional subject list (one ID per line). If omitted, infer from candidate-dir.",
    )
    ap.add_argument(
        "--candidate-dir",
        type=Path,
        required=True,
        help="Directory containing candidate masks (sn_groupmask_in_<ID>.nii[.gz]).",
    )
    ap.add_argument(
        "--anchor-dir",
        type=Path,
        required=True,
        help="Directory containing anchor masks (<prefix>_<ID>.nii[.gz]).",
    )
    ap.add_argument(
        "--anchor-prefixes",
        nargs="+",
        default=["DC_mask_inTSE", "brainstem_mask_inTSE", "brainstem_4v_mask_inTSE"],
        help="Anchor prefixes to use (default: DC, brainstem, 4V brainstem).",
    )
    ap.add_argument(
        "--manual-dir",
        type=Path,
        help="Optional manual mask directory for training labels.",
    )
    ap.add_argument("--candidate-stem", default="sn_groupmask_in")
    ap.add_argument("--manual-stem", default="sn_groupmask_in")
    ap.add_argument(
        "--summary-tsv",
        type=Path,
        help="Optional summary TSV keyed by subject (for metadata columns).",
    )
    ap.add_argument(
        "--candidate-threshold",
        type=float,
        default=0.5,
        help="Binarization threshold for candidate masks (default: 0.5).",
    )
    ap.add_argument(
        "--manual-threshold",
        type=float,
        default=0.5,
        help="Binarization threshold for manual masks (default: 0.5).",
    )
    ap.add_argument(
        "--anchor-threshold",
        type=float,
        default=0.5,
        help="Binarization threshold for anchor masks (default: 0.5).",
    )
    ap.add_argument("--output-tsv", type=Path, required=True)
    return ap.parse_args()


def resolve_mask(mask_dir: Path, stem: str, sid: str) -> Path | None:
    p_gz = mask_dir / f"{stem}_{sid}.nii.gz"
    if p_gz.exists():
        return p_gz
    p_nii = mask_dir / f"{stem}_{sid}.nii"
    if p_nii.exists():
        return p_nii
    return None


def parse_subjects(subjects_file: Path | None, candidate_dir: Path, stem: str) -> List[str]:
    if subjects_file is not None:
        subjects: List[str] = []
        for line in subjects_file.read_text(encoding="utf-8").splitlines():
            tok = line.strip()
            if not tok:
                continue
            subjects.append(tok.split()[0])
        return sorted(set(subjects), key=lambda x: (len(x), x))

    out: List[str] = []
    for path in sorted(candidate_dir.glob(f"{stem}_*.nii*")):
        # Example name: sn_groupmask_in_117.nii.gz
        name = path.name
        base = name.replace(".nii.gz", "").replace(".nii", "")
        if not base.startswith(f"{stem}_"):
            continue
        sid = base[len(stem) + 1 :]
        if sid:
            out.append(sid)
    return sorted(set(out), key=lambda x: (len(x), x))


def load_mask(mask_path: Path, threshold: float) -> tuple[np.ndarray, np.ndarray]:
    img = nib.load(str(mask_path))
    data = img.get_fdata() > threshold
    return data, img.affine


def load_mask_safe(mask_path: Path, threshold: float) -> tuple[np.ndarray | None, np.ndarray | None, str]:
    try:
        data, aff = load_mask(mask_path, threshold)
        return data, aff, ""
    except Exception as exc:
        return None, None, str(exc)


def center_of_mass_mm(mask: np.ndarray, affine: np.ndarray) -> np.ndarray | None:
    idx = np.argwhere(mask)
    if idx.size == 0:
        return None
    com_vox = idx.mean(axis=0)
    return nib.affines.apply_affine(affine, com_vox)


def voxel_volume_mm3(affine: np.ndarray) -> float:
    return float(abs(np.linalg.det(affine[:3, :3])))


def bbox_span_mm(mask: np.ndarray, affine: np.ndarray) -> np.ndarray | None:
    idx = np.argwhere(mask)
    if idx.size == 0:
        return None
    coords_mm = nib.affines.apply_affine(affine, idx)
    mins = coords_mm.min(axis=0)
    maxs = coords_mm.max(axis=0)
    return maxs - mins


def shape_axes_mm(mask: np.ndarray, affine: np.ndarray) -> np.ndarray | None:
    idx = np.argwhere(mask)
    if idx.shape[0] < 3:
        return None
    coords_mm = nib.affines.apply_affine(affine, idx)
    c = coords_mm - coords_mm.mean(axis=0, keepdims=True)
    cov = np.cov(c.T)
    try:
        vals = np.linalg.eigvalsh(cov)
    except np.linalg.LinAlgError:
        return None
    vals = np.clip(vals, 0.0, None)
    axes = np.sqrt(vals[::-1])  # descending
    return axes


def lr_imbalance(mask: np.ndarray, affine: np.ndarray) -> float:
    idx = np.argwhere(mask)
    if idx.size == 0:
        return math.nan
    coords_mm = nib.affines.apply_affine(affine, idx)
    x = coords_mm[:, 0]
    mid = float(np.median(x))
    n_left = int(np.sum(x < mid))
    n_right = int(np.sum(x >= mid))
    total = n_left + n_right
    if total == 0:
        return math.nan
    return float(abs(n_left - n_right) / total)


def dice_coeff(a: np.ndarray, b: np.ndarray) -> float:
    inter = int(np.logical_and(a, b).sum())
    total = int(a.sum()) + int(b.sum())
    if total == 0:
        return math.nan
    return float((2.0 * inter) / total)


def read_summary(summary_tsv: Path | None) -> Dict[str, dict]:
    if summary_tsv is None or (not summary_tsv.exists()):
        return {}
    out: Dict[str, dict] = {}
    with summary_tsv.open(newline="", encoding="utf-8") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            sid = (row.get("subject") or "").strip()
            if sid:
                out[sid] = row
    return out


def as_str(v: float | int | str | None) -> str:
    if v is None:
        return ""
    if isinstance(v, str):
        return v
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if not np.isfinite(v):
            return ""
        return f"{v:.6f}"
    return str(v)


def summarize_anchor_errors(anchor_errors: Sequence[float]) -> tuple[float, float]:
    arr = np.array([x for x in anchor_errors if np.isfinite(x)], dtype=float)
    if arr.size == 0:
        return (math.nan, math.nan)
    return (float(arr.mean()), float(arr.max()))


def main() -> None:
    args = parse_args()

    subjects = parse_subjects(args.subjects_file, args.candidate_dir, args.candidate_stem)
    summary_map = read_summary(args.summary_tsv)
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Header assembly
    header = [
        "subject",
        "candidate_mask",
        "candidate_present",
        "candidate_load_error",
        "candidate_voxels",
        "candidate_volume_mm3",
        "candidate_cm_x_mm",
        "candidate_cm_y_mm",
        "candidate_cm_z_mm",
        "candidate_bbox_x_mm",
        "candidate_bbox_y_mm",
        "candidate_bbox_z_mm",
        "candidate_shape_pc1_mm",
        "candidate_shape_pc2_mm",
        "candidate_shape_pc3_mm",
        "candidate_lr_imbalance",
        "manual_mask",
        "manual_present",
        "manual_load_error",
        "manual_voxels",
        "manual_volume_mm3",
        "manual_cm_x_mm",
        "manual_cm_y_mm",
        "manual_cm_z_mm",
        "target_delta_x_mm",
        "target_delta_y_mm",
        "target_delta_z_mm",
        "target_delta_mag_mm",
        "manual_candidate_dice",
        "anchor_error_mean_mm",
        "anchor_error_max_mm",
    ]
    for anchor in args.anchor_prefixes:
        header.extend(
            [
                f"anchor_present__{anchor}",
                f"anchor_cm_x_mm__{anchor}",
                f"anchor_cm_y_mm__{anchor}",
                f"anchor_cm_z_mm__{anchor}",
                f"vec_x_mm__{anchor}",
                f"vec_y_mm__{anchor}",
                f"vec_z_mm__{anchor}",
                f"dist_mm__{anchor}",
                f"anchor_error_mm__{anchor}",
                f"anchor_load_error__{anchor}",
            ]
        )

    # Add summary columns (if provided) as "summary__<field>"
    summary_fields: List[str] = []
    if summary_map:
        sample = next(iter(summary_map.values()))
        summary_fields = [k for k in sample.keys() if k != "subject"]
        header.extend([f"summary__{k}" for k in summary_fields])

    rows: List[dict] = []
    for sid in subjects:
        row: Dict[str, str] = {"subject": sid}
        cand_path = resolve_mask(args.candidate_dir, args.candidate_stem, sid)
        row["candidate_mask"] = str(cand_path) if cand_path else ""
        row["candidate_present"] = "1" if cand_path else "0"
        row["candidate_load_error"] = ""

        manual_path = None
        if args.manual_dir is not None:
            manual_path = resolve_mask(args.manual_dir, args.manual_stem, sid)
        row["manual_mask"] = str(manual_path) if manual_path else ""
        row["manual_present"] = "1" if manual_path else "0"
        row["manual_load_error"] = ""

        cand_mask = cand_aff = None
        cand_cm = cand_bbox = cand_axes = None
        cand_vox = cand_vol = cand_lr = math.nan
        if cand_path is not None:
            cand_mask, cand_aff, cand_err = load_mask_safe(cand_path, args.candidate_threshold)
            if cand_err:
                row["candidate_present"] = "0"
                row["candidate_load_error"] = cand_err
            elif cand_mask is not None and cand_aff is not None:
                cand_vox = float(cand_mask.sum())
                cand_vol = cand_vox * voxel_volume_mm3(cand_aff)
                cand_cm = center_of_mass_mm(cand_mask, cand_aff)
                cand_bbox = bbox_span_mm(cand_mask, cand_aff)
                cand_axes = shape_axes_mm(cand_mask, cand_aff)
                cand_lr = lr_imbalance(cand_mask, cand_aff)

        row["candidate_voxels"] = as_str(cand_vox)
        row["candidate_volume_mm3"] = as_str(cand_vol)
        row["candidate_cm_x_mm"] = as_str(cand_cm[0] if cand_cm is not None else math.nan)
        row["candidate_cm_y_mm"] = as_str(cand_cm[1] if cand_cm is not None else math.nan)
        row["candidate_cm_z_mm"] = as_str(cand_cm[2] if cand_cm is not None else math.nan)
        row["candidate_bbox_x_mm"] = as_str(cand_bbox[0] if cand_bbox is not None else math.nan)
        row["candidate_bbox_y_mm"] = as_str(cand_bbox[1] if cand_bbox is not None else math.nan)
        row["candidate_bbox_z_mm"] = as_str(cand_bbox[2] if cand_bbox is not None else math.nan)
        row["candidate_shape_pc1_mm"] = as_str(cand_axes[0] if cand_axes is not None else math.nan)
        row["candidate_shape_pc2_mm"] = as_str(cand_axes[1] if cand_axes is not None else math.nan)
        row["candidate_shape_pc3_mm"] = as_str(cand_axes[2] if cand_axes is not None else math.nan)
        row["candidate_lr_imbalance"] = as_str(cand_lr)

        manual_mask = manual_aff = None
        manual_cm = None
        manual_vox = manual_vol = math.nan
        if manual_path is not None:
            manual_mask, manual_aff, manual_err = load_mask_safe(manual_path, args.manual_threshold)
            if manual_err:
                row["manual_present"] = "0"
                row["manual_load_error"] = manual_err
            elif manual_mask is not None and manual_aff is not None:
                manual_vox = float(manual_mask.sum())
                manual_vol = manual_vox * voxel_volume_mm3(manual_aff)
                manual_cm = center_of_mass_mm(manual_mask, manual_aff)

        row["manual_voxels"] = as_str(manual_vox)
        row["manual_volume_mm3"] = as_str(manual_vol)
        row["manual_cm_x_mm"] = as_str(manual_cm[0] if manual_cm is not None else math.nan)
        row["manual_cm_y_mm"] = as_str(manual_cm[1] if manual_cm is not None else math.nan)
        row["manual_cm_z_mm"] = as_str(manual_cm[2] if manual_cm is not None else math.nan)

        target_delta = np.array([math.nan, math.nan, math.nan], dtype=float)
        target_mag = math.nan
        if manual_cm is not None and cand_cm is not None:
            target_delta = manual_cm - cand_cm
            target_mag = float(np.linalg.norm(target_delta))
        row["target_delta_x_mm"] = as_str(float(target_delta[0]))
        row["target_delta_y_mm"] = as_str(float(target_delta[1]))
        row["target_delta_z_mm"] = as_str(float(target_delta[2]))
        row["target_delta_mag_mm"] = as_str(target_mag)

        dice = math.nan
        if manual_mask is not None and cand_mask is not None and manual_mask.shape == cand_mask.shape:
            dice = dice_coeff(manual_mask, cand_mask)
        row["manual_candidate_dice"] = as_str(dice)

        anchor_errors: List[float] = []
        for anchor in args.anchor_prefixes:
            anc_path = resolve_mask(args.anchor_dir, anchor, sid)
            present = anc_path is not None
            row[f"anchor_present__{anchor}"] = "1" if present else "0"
            anc_cm = None
            row[f"anchor_load_error__{anchor}"] = ""
            if anc_path is not None:
                anc_mask, anc_aff, anc_err = load_mask_safe(anc_path, args.anchor_threshold)
                if anc_err:
                    row[f"anchor_present__{anchor}"] = "0"
                    row[f"anchor_load_error__{anchor}"] = anc_err
                elif anc_mask is not None and anc_aff is not None:
                    anc_cm = center_of_mass_mm(anc_mask, anc_aff)
            row[f"anchor_cm_x_mm__{anchor}"] = as_str(anc_cm[0] if anc_cm is not None else math.nan)
            row[f"anchor_cm_y_mm__{anchor}"] = as_str(anc_cm[1] if anc_cm is not None else math.nan)
            row[f"anchor_cm_z_mm__{anchor}"] = as_str(anc_cm[2] if anc_cm is not None else math.nan)

            vec = np.array([math.nan, math.nan, math.nan], dtype=float)
            dist = math.nan
            if cand_cm is not None and anc_cm is not None:
                vec = cand_cm - anc_cm
                dist = float(np.linalg.norm(vec))
            row[f"vec_x_mm__{anchor}"] = as_str(float(vec[0]))
            row[f"vec_y_mm__{anchor}"] = as_str(float(vec[1]))
            row[f"vec_z_mm__{anchor}"] = as_str(float(vec[2]))
            row[f"dist_mm__{anchor}"] = as_str(dist)

            anc_err = math.nan
            if manual_cm is not None and cand_cm is not None and anc_cm is not None:
                manual_vec = manual_cm - anc_cm
                cand_vec = cand_cm - anc_cm
                anc_err = float(np.linalg.norm(manual_vec - cand_vec))
                anchor_errors.append(anc_err)
            row[f"anchor_error_mm__{anchor}"] = as_str(anc_err)

        anc_mean, anc_max = summarize_anchor_errors(anchor_errors)
        row["anchor_error_mean_mm"] = as_str(anc_mean)
        row["anchor_error_max_mm"] = as_str(anc_max)

        if sid in summary_map:
            summ = summary_map[sid]
            for key in summary_fields:
                row[f"summary__{key}"] = summ.get(key, "")
        else:
            for key in summary_fields:
                row[f"summary__{key}"] = ""

        rows.append(row)

    with args.output_tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"WROTE {args.output_tsv} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
