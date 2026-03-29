#!/usr/bin/env python3
"""Build heuristic nigrosome templates from the SN probabilistic atlas.

The script encodes the qualitative anatomy from histology/7T MRI (e.g.,
Blazejewska et al., 2013; Damier et al., 1999) and splits the SN
probabilistic atlas into five per-hemisphere subregions. Each subregion is
assigned by combining dorsal/ventral, lateral/medial, and rostral/caudal
fractions so we can generate canonical N1–N5 masks without hand drawing.
"""

import argparse
import json
from pathlib import Path

import nibabel as nib
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate N1–N5 nigrosome templates from the SN probabilistic atlas.")
    parser.add_argument(
        "--sn-prob",
        required=True,
        type=Path,
        help="Path to SN probabilistic atlas (e.g., MRI_data/ROIs/SN_probabilistic_MNI152/SN_probabilistic.nii)"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory where N1/N2/... NIfTIs will be written"
    )
    parser.add_argument(
        "--sn-threshold",
        type=float,
        default=0.05,
        help="Probability threshold for the SN mask (default: 0.05)"
    )
    parser.add_argument(
        "--dorsal-high",
        type=float,
        default=0.65,
        help="Normalized dorsal fraction for the dorsal-most layer (default: 0.65)"
    )
    parser.add_argument(
        "--dorsal-mid",
        type=float,
        default=0.45,
        help="Normalized dorsal fraction for the mid layer boundary (default: 0.45)"
    )
    parser.add_argument(
        "--lateral-high",
        type=float,
        default=0.50,
        help="Normalized lateral fraction that separates lateral vs medial zones (default: 0.50)"
    )
    parser.add_argument(
        "--lateral-mid",
        type=float,
        default=0.35,
        help="Normalized lateral fraction for the ventral split (default: 0.35)"
    )
    parser.add_argument(
        "--min-rostral",
        type=float,
        default=0.15,
        help="Ignored rostral fraction below this value (default: 0.15)"
    )
    parser.add_argument(
        "--max-rostral",
        type=float,
        default=0.90,
        help="Ignored rostral fraction above this value (default: 0.90)"
    )
    parser.add_argument(
        "--prefix",
        default="SN_",
        help="Filename prefix for outputs (default: SN_)"
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        help="Optional path to write JSON summary of voxel counts"
    )
    return parser.parse_args()


def normalize(coords, min_val, max_val):
    if max_val <= min_val:
        return np.zeros_like(coords, dtype=np.float32)
    return (coords - min_val) / (max_val - min_val)


def compute_axes(indices):
    return np.min(indices), np.max(indices)


def assign_regions(sn_prob, sn_mask, hemi_mask, axes, params):
    assigned = np.zeros(sn_mask.shape, dtype=bool)
    labels = {}

    x_min, x_max = axes['x']
    y_min, y_max = axes['y']
    z_min, z_max = axes['z']
    mid_x = (x_min + x_max) / 2.0

    yy, xx, zz = np.indices(sn_mask.shape)
    hemi = hemi_mask & sn_mask
    if not np.any(hemi):
        return labels

    hemi_idxs = np.where(hemi)
    x = xx[hemi]
    y = yy[hemi]
    z = zz[hemi]

    if hemi_mask.mean() == 0:
        raise ValueError("Hemisphere mask empty")

    if hemi_mask[0, 0, 0] is np.ma.masked:
        pass

    lat = np.abs(x - mid_x)
    lat_norm = normalize(lat, 0, np.max(lat) if np.max(lat) > 0 else 1)
    dor_norm = normalize(z, z_min, z_max)
    ros_norm = normalize(y, y_min, y_max)

    valid_ros = (ros_norm >= params['min_rostral']) & (ros_norm <= params['max_rostral'])

    def make_label(name, condition):
        mask = np.zeros_like(sn_mask, dtype=bool)
        mask[hemi_idxs] = condition
        labels[name] = mask

    cond_n1 = (dor_norm >= params['dorsal_high']) & (lat_norm >= params['lateral_high']) & valid_ros
    make_label('N1', cond_n1)
    assigned |= labels['N1']

    cond_n2 = (~assigned[hemi_idxs]) & (dor_norm >= params['dorsal_mid']) & (lat_norm >= params['lateral_high']) & valid_ros
    make_label('N2', cond_n2)
    assigned |= labels['N2']

    cond_n3 = (~assigned[hemi_idxs]) & (dor_norm >= params['dorsal_mid']) & valid_ros
    cond_n3 &= (lat_norm < params['lateral_high'])
    make_label('N3', cond_n3)
    assigned |= labels['N3']

    cond_n4 = (~assigned[hemi_idxs]) & (dor_norm < params['dorsal_mid']) & (lat_norm >= params['lateral_mid']) & valid_ros
    make_label('N4', cond_n4)
    assigned |= labels['N4']

    cond_n5 = (~assigned[hemi_idxs])
    make_label('N5', cond_n5)
    assigned |= labels['N5']

    return labels


def main():
    args = parse_args()
    sn_img = nib.load(str(args.sn_prob))
    sn_data = sn_img.get_fdata().astype(np.float32)

    sn_mask = sn_data >= args.sn_threshold
    if not np.any(sn_mask):
        raise RuntimeError("SN probabilistic mask is empty after thresholding")

    coords = np.array(np.where(sn_mask))
    axes = {
        'x': (coords[1].min(), coords[1].max()),
        'y': (coords[0].min(), coords[0].max()),
        'z': (coords[2].min(), coords[2].max()),
    }

    mid_x = (axes['x'][0] + axes['x'][1]) / 2.0
    left_mask = np.zeros_like(sn_mask)
    right_mask = np.zeros_like(sn_mask)
    left_mask[:, :int(np.floor(mid_x)) + 1, :] = True
    right_mask[:, int(np.floor(mid_x)) + 1:, :] = True

    params = {
        'dorsal_high': args.dorsal_high,
        'dorsal_mid': args.dorsal_mid,
        'lateral_high': args.lateral_high,
        'lateral_mid': args.lateral_mid,
        'min_rostral': args.min_rostral,
        'max_rostral': args.max_rostral,
    }

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    summary = {}
    hemi_info = [('L', left_mask), ('R', right_mask)]
    for hemi, mask in hemi_info:
        labels = assign_regions(sn_data, sn_mask, mask, axes, params)
        hemi_summary = {}
        for name, region_mask in labels.items():
            prob_region = np.zeros_like(sn_data, dtype=np.float32)
            prob_region[region_mask] = sn_data[region_mask]
            out_path = output_dir / f"{args.prefix}{name}_{hemi}.nii.gz"
            nib.Nifti1Image(prob_region, sn_img.affine, sn_img.header).to_filename(str(out_path))
            hemi_summary[name] = {
                'voxels': int(region_mask.sum()),
                'mean_probability': float(sn_data[region_mask].mean() if np.any(region_mask) else 0.0)
            }
        summary[hemi] = hemi_summary

    meta = {
        'sn_prob': str(args.sn_prob),
        'sn_threshold': args.sn_threshold,
        'parameters': params,
        'summary': summary,
    }
    if args.metadata:
        args.metadata.parent.mkdir(parents=True, exist_ok=True)
        args.metadata.write_text(json.dumps(meta, indent=2))
    else:
        print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
