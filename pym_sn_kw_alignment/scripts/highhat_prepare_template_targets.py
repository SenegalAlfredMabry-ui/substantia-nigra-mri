#!/usr/bin/env python3
"""Compute template-space ROI targets from KW manual masks for Highhat."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import tempfile
import nibabel as nib
import numpy as np

def read_subjects(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]

def compute_centroid_mm(mask_path: Path) -> np.ndarray:
    img = nib.load(str(mask_path))
    data = img.get_fdata()
    coords = np.argwhere(data > 0)
    if coords.size == 0:
        raise ValueError(f"Mask {mask_path} empty")
    centroid_vox = coords.mean(axis=0)
    centroid_mm = nib.affines.apply_affine(img.affine, centroid_vox)
    return centroid_mm

def load_affine(path: Path) -> np.ndarray:
    vals = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            vals.extend(float(x) for x in line.split())
    if len(vals) != 12:
        raise ValueError(f"Unexpected affine length in {path}")
    mat = np.eye(4)
    mat[:3, :4] = np.array(vals).reshape(3, 4)
    return mat


def run_ants_apply(points_mm: np.ndarray, transforms: list[str]) -> np.ndarray:
    with tempfile.TemporaryDirectory() as tmpdir:
        in_path = Path(tmpdir) / "points_in.csv"
        out_path = Path(tmpdir) / "points_out.csv"
        with in_path.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["x", "y", "z"])
            writer.writerow(points_mm.tolist())
        cmd = [
            "antsApplyTransformsToPoints",
            "-d", "3",
            "-i", str(in_path),
            "-o", str(out_path),
        ]
        for t in transforms:
            cmd.extend(["-t", t])
        subprocess.run(cmd, check=True)
        with out_path.open() as fh:
            reader = csv.DictReader(fh)
            # skip header row? DictReader handles header line.
            row = next(reader)
            return np.array([float(row["output_x"]), float(row["output_y"]), float(row["output_z"])])

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subjects", required=True, help="Text file with subject IDs")
    parser.add_argument("--manual-dir", required=True)
    parser.add_argument("--warps-root", required=True)
    parser.add_argument("--motip-tse-root", required=True)
    parser.add_argument("--template", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    manual_dir = Path(args.manual_dir)
    warps_root = Path(args.warps_root)
    motip_root = Path(args.motip_tse_root)
    template_img = nib.load(args.template)

    targets = {}
    for sub in read_subjects(Path(args.subjects)):
        mask_path = manual_dir / f"sn_groupmask_in_{sub}.nii.gz"
        if not mask_path.exists():
            print(f"[WARN] manual mask missing for {sub}")
            continue
        centroid_tse_mm = compute_centroid_mm(mask_path)
        affine_path = motip_root / f"{sub}_anat_toTSE_mat.aff12.1D"
        if not affine_path.exists():
            print(f"[WARN] Missing anat_toTSE affine for {sub}: {affine_path}")
            continue
        affine = load_affine(affine_path)
        try:
            inv_affine = np.linalg.inv(affine)
        except np.linalg.LinAlgError:
            print(f"[WARN] Singular affine for {sub}")
            continue
        centroid_anat_mm = nib.affines.apply_affine(inv_affine, centroid_tse_mm)
        transforms = [
            str(warps_root / sub / "anat_to_template_antpunk_1Warp.nii.gz"),
            str(warps_root / sub / "anat_to_template_antpunk_0GenericAffine.mat"),
        ]
        try:
            tpl_mm = run_ants_apply(centroid_anat_mm, transforms)
        except subprocess.CalledProcessError as exc:
            print(f"[WARN] antsApplyTransforms failed for {sub}: {exc}")
            continue
        tpl_vox = np.linalg.solve(template_img.affine[:3, :3], tpl_mm - template_img.affine[:3, 3])
        targets[sub] = {
            "template_mm": tpl_mm.tolist(),
            "template_vox": tpl_vox.tolist(),
            "kw_tse_mm": centroid_tse_mm.tolist(),
            "anat_mm": centroid_anat_mm.tolist(),
        }
    Path(args.output).write_text(json.dumps(targets, indent=2))

if __name__ == "__main__":
    main()
