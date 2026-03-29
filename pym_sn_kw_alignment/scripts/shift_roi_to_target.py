#!/usr/bin/env python3
"""Shift a ROI so its centroid matches a target voxel coordinate using sub-voxel interpolation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import nibabel as nib
import numpy as np
from scipy.ndimage import shift as nd_shift

EPS = 1e-3
MAX_ITERS = 5
PADDING = 10

def compute_centroid(data: np.ndarray) -> np.ndarray:
    coords = np.argwhere(data > 0)
    if coords.size == 0:
        raise RuntimeError("Input ROI is empty")
    return coords.mean(axis=0)

def pad_volume(data: np.ndarray, pad: int) -> np.ndarray:
    return np.pad(data, pad_width=pad, mode='constant', constant_values=0)

def crop_volume(data: np.ndarray, pad: int, shape: tuple) -> np.ndarray:
    slices = tuple(slice(pad, pad + dim) for dim in shape)
    return data[slices]

def determine_target_voxel(
    data: np.ndarray,
    affine: np.ndarray,
    target_arg: list[float],
    offset_json: Path | None,
    subject: str | None,
) -> np.ndarray:
    """Resolve the voxel target from either direct coords or per-subject offsets."""
    centroid = compute_centroid(data)
    if offset_json is not None:
        if subject is None:
            raise SystemExit("--subject is required when --offset-json is used")
        offsets = json.loads(offset_json.read_text())
        entry = offsets.get(subject)
        if entry is None:
            raise SystemExit(f"Subject {subject} missing from {offset_json}")
        delta_mm = np.array(entry.get("delta_mm"))
        if delta_mm.shape != (3,):
            raise SystemExit(f"Invalid delta_mm for {subject} in {offset_json}")
        delta_vox = np.linalg.solve(affine[:3, :3], delta_mm)
        return centroid + delta_vox
    if target_arg:
        if len(target_arg) != 3:
            raise SystemExit("Target voxel must have exactly 3 components")
        return np.array(target_arg, dtype=float)
    raise SystemExit("Provide either target voxel coordinates or --offset-json + --subject")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Input ROI NIfTI")
    parser.add_argument("output", help="Output NIfTI path")
    parser.add_argument("target", nargs="*", type=float, help="Target voxel coordinate (x y z)")
    parser.add_argument("--offset-json", dest="offset_json", help="JSON file containing per-subject delta_mm entries")
    parser.add_argument("--subject", help="Subject ID when using --offset-json")
    args = parser.parse_args()

    img = nib.load(args.input)
    data = img.get_fdata()
    offset_path = Path(args.offset_json) if args.offset_json else None
    target = determine_target_voxel(data, img.affine, args.target, offset_path, args.subject)
    pad = PADDING
    padded = pad_volume(data, pad)
    shifted = padded
    for _ in range(MAX_ITERS):
        centroid = compute_centroid(shifted)
        delta = target + pad - centroid
        if np.linalg.norm(delta) < EPS:
            break
        shifted = nd_shift(shifted, shift=delta, order=0, mode='constant', cval=0.0, prefilter=False)
    cropped = crop_volume(shifted, pad, data.shape)
    out_img = nib.Nifti1Image(cropped.astype(img.get_data_dtype()), img.affine, img.header)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    nib.save(out_img, args.output)

if __name__ == '__main__':
    main()
