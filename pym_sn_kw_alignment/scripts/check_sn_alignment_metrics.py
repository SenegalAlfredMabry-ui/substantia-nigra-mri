#!/usr/bin/env python3
"""Compute overlap and centroid offsets between manual SN tracings and baseline masks.

Example:
    python Scripts/check_sn_alignment_metrics.py \
        --subjects 132 223 564 \
        --manual-dir MRI_data/TSE/probabilistic_masks_v9_manual/base \
        --baseline-dir MRI_data/TSE/probabilistic_masks_v9/base \
        --output tmp/sn_alignment_metrics.tsv
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import nibabel as nib
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare manual SN masks with baseline auto warps."
    )
    parser.add_argument(
        "--subjects",
        nargs="*",
        help="Subject IDs (space/comma separated). May be omitted if --subjects-file is set.",
    )
    parser.add_argument(
        "--subjects-file",
        help="Path to a text file containing one subject ID per line.",
    )
    parser.add_argument(
        "--manual-dir",
        required=True,
        help="Directory containing manual SN masks (sn_groupmask_in_<ID>.nii[.gz]).",
    )
    parser.add_argument(
        "--baseline-dir",
        required=True,
        help="Directory containing baseline SN masks in native space.",
    )
    parser.add_argument(
        "--output",
        help="Optional TSV output path. Defaults to stdout if not provided.",
    )
    parser.add_argument(
        "--mask-threshold",
        type=float,
        default=0.5,
        help="Binarization threshold for manual/baseline masks (default: 0.5).",
    )
    return parser.parse_args()


def load_subjects(subject_args: Optional[Sequence[str]], subjects_file: Optional[str]) -> List[str]:
    subjects: List[str] = []
    if subjects_file:
        path = Path(subjects_file)
        if not path.exists():
            raise FileNotFoundError(f"Subjects file not found: {path}")
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                token = line.strip()
                if token:
                    subjects.append(token)
    if subject_args:
        for token in subject_args:
            for part in token.split(","):
                part = part.strip()
                if part:
                    subjects.append(part)
    deduped = sorted(dict.fromkeys(subjects))
    if not deduped:
        raise ValueError("No subjects provided; set --subjects and/or --subjects-file.")
    return deduped


def load_mask(path: Path, threshold: float) -> Tuple[np.ndarray, np.ndarray]:
    img = nib.load(str(path))
    data = img.get_fdata()
    mask = data > threshold
    return mask, img.affine


def center_of_mass(mask: np.ndarray) -> Optional[np.ndarray]:
    coords = np.argwhere(mask)
    if coords.size == 0:
        return None
    return coords.mean(axis=0)


def voxel_volume_mm3(affine: np.ndarray) -> float:
    return abs(np.linalg.det(affine[:3, :3]))


def dice_coeff(a: np.ndarray, b: np.ndarray) -> float:
    inter = np.logical_and(a, b).sum()
    total = a.sum() + b.sum()
    if total == 0:
        return math.nan
    return (2.0 * inter) / total


def format_float(value: Optional[float]) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "nan"
    return f"{value:.6f}"


def collect_metrics(
    subjects: Sequence[str],
    manual_dir: Path,
    baseline_dir: Path,
    mask_threshold: float,
) -> List[List[str]]:
    rows: List[List[str]] = []
    for sub in subjects:
        manual_path = manual_dir / f"sn_groupmask_in_{sub}.nii.gz"
        if not manual_path.exists():
            manual_path = manual_dir / f"sn_groupmask_in_{sub}.nii"
        baseline_path = baseline_dir / f"sn_groupmask_in_{sub}.nii.gz"
        if not baseline_path.exists():
            baseline_path = baseline_dir / f"sn_groupmask_in_{sub}.nii"

        manual_mask = baseline_mask = manual_affine = baseline_affine = None
        notes: List[str] = []
        if manual_path.exists():
            manual_mask, manual_affine = load_mask(manual_path, mask_threshold)
        else:
            notes.append("manual_missing")
        if baseline_path.exists():
            baseline_mask, baseline_affine = load_mask(baseline_path, mask_threshold)
        else:
            notes.append("baseline_missing")

        dice = math.nan
        manual_vox = baseline_vox = math.nan
        manual_vol = baseline_vol = math.nan
        cm_manual_world = cm_baseline_world = None
        delta = (math.nan, math.nan, math.nan)
        delta_mag = math.nan

        if manual_mask is not None:
            manual_vox = float(manual_mask.sum())
            manual_vol = manual_vox * voxel_volume_mm3(manual_affine)
            cm_manual_vox = center_of_mass(manual_mask)
            if cm_manual_vox is not None:
                cm_manual_world = nib.affines.apply_affine(manual_affine, cm_manual_vox)

        if baseline_mask is not None:
            baseline_vox = float(baseline_mask.sum())
            baseline_vol = baseline_vox * voxel_volume_mm3(baseline_affine)
            cm_baseline_vox = center_of_mass(baseline_mask)
            if cm_baseline_vox is not None:
                cm_baseline_world = nib.affines.apply_affine(baseline_affine, cm_baseline_vox)

        if manual_mask is not None and baseline_mask is not None:
            if manual_mask.shape != baseline_mask.shape:
                notes.append(f"shape_mismatch_manual{manual_mask.shape}_baseline{baseline_mask.shape}")
            else:
                dice = dice_coeff(manual_mask, baseline_mask)
            if cm_manual_world is not None and cm_baseline_world is not None:
                delta_vec = cm_manual_world - cm_baseline_world
                delta = tuple(float(x) for x in delta_vec)
                delta_mag = float(np.linalg.norm(delta_vec))

        rows.append(
            [
                sub,
                format_float(dice),
                format_float(manual_vox),
                format_float(baseline_vox),
                format_float(manual_vol),
                format_float(baseline_vol),
                format_float(cm_manual_world[0] if cm_manual_world is not None else math.nan),
                format_float(cm_manual_world[1] if cm_manual_world is not None else math.nan),
                format_float(cm_manual_world[2] if cm_manual_world is not None else math.nan),
                format_float(cm_baseline_world[0] if cm_baseline_world is not None else math.nan),
                format_float(cm_baseline_world[1] if cm_baseline_world is not None else math.nan),
                format_float(cm_baseline_world[2] if cm_baseline_world is not None else math.nan),
                format_float(delta[0]),
                format_float(delta[1]),
                format_float(delta[2]),
                format_float(delta_mag),
                ";".join(notes) if notes else "",
            ]
        )
    return rows


def write_rows(rows: Iterable[List[str]], output: Optional[Path]) -> None:
    header = [
        "subject",
        "dice",
        "manual_voxels",
        "baseline_voxels",
        "manual_volume_mm3",
        "baseline_volume_mm3",
        "cm_manual_x_mm",
        "cm_manual_y_mm",
        "cm_manual_z_mm",
        "cm_baseline_x_mm",
        "cm_baseline_y_mm",
        "cm_baseline_z_mm",
        "delta_x_mm",
        "delta_y_mm",
        "delta_z_mm",
        "delta_mag_mm",
        "notes",
    ]
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        handle = output.open("w", newline="", encoding="utf-8")
    else:
        handle = sys.stdout
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(header)
    writer.writerows(rows)
    if output:
        handle.close()


def main() -> None:
    args = parse_args()
    subjects = load_subjects(args.subjects, args.subjects_file)
    manual_dir = Path(args.manual_dir)
    baseline_dir = Path(args.baseline_dir)
    if not manual_dir.is_dir():
        raise FileNotFoundError(f"Manual directory not found: {manual_dir}")
    if not baseline_dir.is_dir():
        raise FileNotFoundError(f"Baseline directory not found: {baseline_dir}")
    rows = collect_metrics(subjects, manual_dir, baseline_dir, args.mask_threshold)
    output_path = Path(args.output) if args.output else None
    write_rows(rows, output_path)


if __name__ == "__main__":
    main()
