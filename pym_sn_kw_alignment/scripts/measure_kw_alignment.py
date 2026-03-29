#!/usr/bin/env python3
"""Compute KW vs v9/v10 mask overlap metrics for Highhat subjects."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import nibabel as nib
import numpy as np


SUBJECTS: Iterable[str] = [
    "120",
    "128",
    "142",
    "151",
    "166",
    "201",
    "206",
    "210",
    "219",
    "226",
    "234",
    "239",
]

QA_ROOT = Path("QA") / "manual_alignment_review"
OUTPUT_CSV = QA_ROOT / "highhat_kw_alignment.csv"


def load_mask(path: Path) -> tuple[np.ndarray, np.ndarray]:
    img = nib.load(str(path))
    data = img.get_fdata()
    mask = data > 0
    return mask, img.affine


def centroid_mm(mask: np.ndarray, affine: np.ndarray) -> np.ndarray | None:
    coords = np.argwhere(mask)
    if coords.size == 0:
        return None
    vox_centroid = coords.mean(axis=0)
    return nib.affines.apply_affine(affine, vox_centroid)


def overlap_metrics(reference: np.ndarray, target: np.ndarray) -> tuple[int, int, int]:
    ref_count = int(reference.sum())
    tgt_count = int(target.sum())
    inter = int(np.logical_and(reference, target).sum())
    return ref_count, tgt_count, inter


def dice(ref_count: int, tgt_count: int, inter: int) -> float | None:
    denom = ref_count + tgt_count
    if denom == 0:
        return None
    return 2.0 * inter / denom


def jaccard(ref_count: int, tgt_count: int, inter: int) -> float | None:
    union = ref_count + tgt_count - inter
    if union == 0:
        return None
    return inter / union


def main() -> None:
    rows: list[list[str | float | int | None]] = []
    rows.append(
        [
            "subject",
            "kw_vox",
            "v10_vox",
            "v10_dice",
            "v10_jaccard",
            "v10_centroid_mm",
            "v10_centroid_distance_mm",
            "v9_vox",
            "v9_dice",
            "v9_jaccard",
            "v9_centroid_mm",
            "v9_centroid_distance_mm",
        ]
    )

    for sub in SUBJECTS:
        qa_dir = QA_ROOT / sub
        kw_path = qa_dir / "kw_tracings" / f"SN_ROI_KW_{sub}_ZP.nii"
        v10_path = qa_dir / "black_ant_masks_v10" / f"sn_groupmask_in_{sub}.nii.gz"
        v9_path = qa_dir / f"sn_groupmask_in_{sub}_v9.nii.gz"

        if not kw_path.exists():
            print(f"[WARN] Missing KW mask for {sub}: {kw_path}")
            continue
        if not v10_path.exists():
            print(f"[WARN] Missing v10 mask for {sub}: {v10_path}")
            continue
        if not v9_path.exists():
            print(f"[WARN] Missing v9 mask for {sub}: {v9_path}")
            continue

        kw_mask, kw_affine = load_mask(kw_path)
        kw_centroid = centroid_mm(kw_mask, kw_affine)

        v10_mask, v10_affine = load_mask(v10_path)
        v9_mask, v9_affine = load_mask(v9_path)

        v10_counts = overlap_metrics(kw_mask, v10_mask)
        v9_counts = overlap_metrics(kw_mask, v9_mask)

        v10_centroid = centroid_mm(v10_mask, v10_affine)
        v9_centroid = centroid_mm(v9_mask, v9_affine)

        def centroid_distance(other: np.ndarray | None) -> float | None:
            if kw_centroid is None or other is None:
                return None
            return float(np.linalg.norm(kw_centroid - other))

        rows.append(
            [
                sub,
                v10_counts[0],
                v10_counts[1],
                dice(*v10_counts),
                jaccard(*v10_counts),
                ";".join(f"{x:.2f}" for x in v10_centroid) if v10_centroid is not None else "",
                centroid_distance(v10_centroid),
                v9_counts[1],
                dice(*v9_counts),
                jaccard(*v9_counts),
                ";".join(f"{x:.2f}" for x in v9_centroid) if v9_centroid is not None else "",
                centroid_distance(v9_centroid),
            ]
        )

    with OUTPUT_CSV.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerows(rows)
    print(f"[INFO] Wrote alignment metrics to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
