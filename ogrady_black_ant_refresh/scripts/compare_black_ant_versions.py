#!/usr/bin/env python3
"""Compare legacy Black ANT V1 metrics against the Irredeemable rerun."""

import os
from pathlib import Path

import pandas as pd

# Default lookup order mirrors the collaborator's request: keep the original
# V1 autosave, then append both the thresholded and unthresholded outputs from
# the Irredeemable rerun.
DEFAULT_PATHS = {
    'V1': 'MRI_data/analysis/SN_nmdata_autosaved.txt',
    'IRR': 'MRI_data/analysis/SN_nmdata_irredeemable_*_thr*.txt',
}

BASE_COLS = ['mean_raw_sn', 'mean_raw_control', 'std_raw_control', 'CNR', 'intensity']


def resolve_paths():
    resolved = {}
    for label, default in DEFAULT_PATHS.items():
        env_key = f'BLACKANT_{label}_TABLE'
        spec = os.environ.get(env_key, default)
        if '*' in spec:
            matches = sorted(
                Path().glob(spec),
                key=lambda p: p.stat().st_mtime if p.exists() else 0,
                reverse=True,
            )
            resolved[label] = matches[0] if matches else None
        else:
            resolved[label] = Path(spec)
    return resolved


def load_table(label, csv_path):
    if not csv_path.exists():
        raise SystemExit(f"Missing file for {label}: {csv_path}")
    df = pd.read_csv(csv_path)
    cols_needed = [c for c in BASE_COLS if c in df.columns]
    keep_cols = ['sub'] + cols_needed
    renamed = df[keep_cols].rename(columns={c: f"{c}_{label}" for c in keep_cols if c != 'sub'})
    return renamed


def load_irredeemable(csv_path):
    if not csv_path.exists():
        raise SystemExit(f"Missing Irredeemable file: {csv_path}")
    df = pd.read_csv(csv_path)
    thr_cols = [c for c in df.columns if c not in {'sub', 'notes'} and not c.endswith('_unthresholded')]
    unthr_cols = [c for c in df.columns if c.endswith('_unthresholded')]

    thr = df[['sub'] + thr_cols].copy()
    thr.rename(columns={c: f"{c}_IRR_thr" for c in thr_cols}, inplace=True)

    unthr = df[['sub'] + unthr_cols].copy()
    unthr.rename(columns={c: f"{c[:-14]}_IRR_unthr" for c in unthr_cols}, inplace=True)

    merged = thr.merge(unthr, on='sub', how='outer')
    return merged


def main():
    paths = resolve_paths()
    merged = None

    if (v1_path := paths.get('V1')) is not None:
        table = load_table('V1', v1_path)
        merged = table if merged is None else merged.merge(table, on='sub', how='outer')

    if (irr_path := paths.get('IRR')) is not None:
        irr_table = load_irredeemable(irr_path)
        merged = irr_table if merged is None else merged.merge(irr_table, on='sub', how='outer')

    if merged is None:
        raise SystemExit('No tables were loaded; check file paths or overrides.')

    output = Path('tmp_black_ant_V1_IRR_compare.csv')
    merged.to_csv(output, index=False)
    print(f"Wrote {output} ({len(merged)} rows)")


if __name__ == '__main__':
    main()
