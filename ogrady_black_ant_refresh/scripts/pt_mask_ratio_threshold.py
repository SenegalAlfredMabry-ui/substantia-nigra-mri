#!/usr/bin/env python3
"""Apply ratio-based PT-control mask trimming for selected subjects."""
import argparse
import csv
from pathlib import Path

import nibabel as nib
import numpy as np


def load_mask(mask_dir: Path, prefix: str, sub: str) -> nib.Nifti1Image:
    for ext in ('.nii.gz', '.nii'):
        path = mask_dir / f"{prefix}_groupmask_in_{sub}{ext}"
        if path.exists():
            return nib.load(str(path))
    raise FileNotFoundError(f"Missing {prefix} mask for {sub}")


def ratio_trim(sub: str, mask_dir: Path, out_dir: Path, sn_prob: float,
               ratio: float, min_prob: float, max_cut_frac: float,
               min_prob_floor: float) -> dict:
    sn_img = load_mask(mask_dir, 'sn', sub)
    sn_data = sn_img.get_fdata()
    sn_mask = sn_data >= sn_prob
    sn_count = int(sn_mask.sum())
    if sn_count == 0:
        return {'subject': sub, 'status': 'no_sn_voxels'}
    slice_mask = np.any(np.any(sn_mask, axis=0), axis=0)
    if not np.any(slice_mask):
        return {'subject': sub, 'status': 'no_sn_slices'}
    pt_img = load_mask(mask_dir, 'ptctrl', sub)
    pt_data = pt_img.get_fdata()
    pt_slice = pt_data[:, :, slice_mask]
    pt_vals = pt_slice[pt_slice > 0]
    if pt_vals.size == 0:
        return {'subject': sub, 'status': 'pt_empty'}
    total_pt = int(pt_vals.size)
    min_allowed = max(1, int(round(sn_count * (1 - max_cut_frac))))
    target = int(round(sn_count * ratio))
    target = max(min_allowed, target)
    target = min(target, total_pt)
    sorted_vals = np.sort(pt_vals)[::-1]
    idx = max(0, min(total_pt - 1, target - 1))
    threshold = max(sorted_vals[idx], min_prob)
    final_mask = (pt_data >= threshold)
    sn_volume_mask = np.zeros_like(pt_data, dtype=bool)
    sn_volume_mask[:, :, slice_mask] = True
    final_mask &= sn_volume_mask
    final_count = int(final_mask.sum())
    if final_count < min_allowed:
        idx = max(0, min(total_pt - 1, min_allowed - 1))
        threshold = max(sorted_vals[idx], min_prob)
        final_mask = (pt_data >= threshold) & sn_volume_mask
        final_count = int(final_mask.sum())
    if final_count < min_allowed and min_prob > min_prob_floor:
        threshold = min_prob_floor
        final_mask = (pt_data >= threshold) & sn_volume_mask
        final_count = int(final_mask.sum())
        status = 'minprob_floor'
    else:
        status = 'ok'
    trimmed = np.where(final_mask, pt_data, 0.0)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"ptctrl_ratio_mask_in_{sub}.nii.gz"
    nib.save(nib.Nifti1Image(trimmed, pt_img.affine, pt_img.header), str(out_path))
    return {
        'subject': sub,
        'status': status,
        'sn_voxels': sn_count,
        'pt_voxels_initial': total_pt,
        'pt_voxels_final': final_count,
        'threshold': float(threshold)
    }


def parse_args():
    p = argparse.ArgumentParser(description="Ratio-trim PT-control masks")
    p.add_argument('--mask-dir', default='MRI_data/TSE/probabilistic_masks_v6',
                   type=Path, help='Directory with sn/ptctrl masks')
    p.add_argument('--output-dir', default='MRI_data/TSE/probabilistic_masks_v6_ratio',
                   type=Path, help='Destination for trimmed PT masks')
    p.add_argument('--subjects', nargs='+', required=True,
                   help='Subjects to process (e.g., 142 503 ...)')
    p.add_argument('--sn-prob', type=float, default=0.15,
                   help='SN probability cutoff used for voxel counts')
    p.add_argument('--ratio', type=float, default=1.0,
                   help='Desired PT/SN voxel ratio (>=0)')
    p.add_argument('--max-cut-frac', type=float, default=0.15,
                   help='Maximum fraction of SN voxels the PT mask may drop below')
    p.add_argument('--min-prob', type=float, default=0.01,
                   help='Minimum PT probability threshold when matching ratio')
    p.add_argument('--min-prob-floor', type=float, default=0.005,
                   help='Absolute floor to fall back on when ratio cannot be met')
    p.add_argument('--summary', type=Path,
                   default=Path('MRI_data/analysis/tests/v6/pt_ratio_summary.csv'),
                   help='CSV file for summary stats')
    return p.parse_args()


def main():
    args = parse_args()
    results = []
    for sub in args.subjects:
        try:
            res = ratio_trim(sub, args.mask_dir, args.output_dir, args.sn_prob,
                             args.ratio, args.min_prob, args.max_cut_frac,
                             args.min_prob_floor)
        except FileNotFoundError as exc:
            res = {'subject': sub, 'status': f'missing:{exc}'}
        results.append(res)
        print(sub, res['status'])
    if args.summary:
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = ['subject', 'status', 'sn_voxels', 'pt_voxels_initial',
                      'pt_voxels_final', 'threshold']
        with args.summary.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow(row)


if __name__ == '__main__':
    main()
