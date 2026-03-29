#!/usr/bin/env python3
"""Quick sanity checks for SN Irredeemable metrics tables.

Reads one or more CSVs exported by Irredeemable (V1 or V2) and prints for each:
  * Correlations between CNR and scaled intensity (thresholded + unthresholded)
  * Distribution summaries and extreme subjects for CNR/intensity/sn_voxels
  * Note counts and any subjects still missing data
  * Optional cleanup-table comparisons (main run only)
  * Mean left-right and rostral-caudal differences (Irredeemable V2 columns)

Example:
  python3 Scripts/check_sn_metrics.py \
      --main MRI_data/analysis/SN_nmdata_irredeemable_20251210_thr25.txt \
      --cleanup MRI_data/analysis/SN_nmdata_irredeemable_20251210_thr25_cleanup.txt \
      --expect 145 --watch 530
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics as stats
from collections import Counter
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


VOXEL_VOLUME_MM3 = 0.43 * 0.43 * 3  # Expected SN voxel volume used by Irredeemable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main", required=True, help="Path to primary Irredeemable CSV")
    parser.add_argument(
        "--extra",
        action="append",
        default=[],
        help="Additional run(s) to compare in the form label=path or just a path (label defaults to filename)",
    )
    parser.add_argument("--cleanup", help="Optional cleanup-only CSV")
    parser.add_argument("--expect", type=int, default=0, help="Expected row count (sanity check)")
    parser.add_argument(
        "--watch",
        nargs="*",
        default=[],
        help="Subject IDs that previously failed (report whether they now have data)",
    )
    parser.add_argument(
        "--check-overlap",
        action="store_true",
        help="Verify SN vs SN-control masks have zero overlap (requires mask directory)",
    )
    parser.add_argument(
        "--mask-dir",
        default="MRI_data/TSE/probabilistic_masks_v6",
        help="Directory containing sn/snctrl_groupmask_in_<ID>.nii(.gz) (default: %(default)s)",
    )
    parser.add_argument(
        "--brainstem-dir",
        default="MRI_data/ants_processing/upsampled/TSE",
        help="Directory containing DC_mask_inTSE_<ID>.nii(.gz) for brainstem overlap checks",
    )
    parser.add_argument(
        "--analysis-dir",
        default="MRI_data/analysis",
        help="Directory searched when --recent-tests is enabled (default: %(default)s)",
    )
    parser.add_argument(
        "--recent-tests",
        type=int,
        default=0,
        help="Automatically append the newest Irredeemable test outputs from --analysis-dir",
    )
    parser.add_argument(
        "--recent-tests-pattern",
        default="SN_nmdata_irredeemable_*_test*.txt",
        help="Glob (relative to --analysis-dir) matched for --recent-tests",
    )
    parser.add_argument(
        "--sn-overlap-threshold",
        type=float,
        default=0.15,
        help="Minimum SN probability used when testing mask overlaps (default: %(default)s)",
    )
    parser.add_argument(
        "--pt-variant-root",
        help="Optional directory containing Irredeemable PT-variant tables (summarize PT CNR vs intensity correlations)",
    )
    parser.add_argument(
        "--pt-variant-pattern",
        default="SN_nmdata_irredeemable_v8_*prob*.txt",
        help="Glob pattern searched within --pt-variant-root when summarizing PT variants (default: %(default)s)",
    )
    return parser.parse_args()


def load_rows(path: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def finite_pair(rows: Iterable[Dict[str, str]], key_x: str, key_y: str) -> Tuple[List[float], List[float]]:
    xs, ys = [], []
    for row in rows:
        x = to_float(row.get(key_x, ""))
        y = to_float(row.get(key_y, ""))
        if math.isfinite(x) and math.isfinite(y):
            xs.append(x)
            ys.append(y)
    return xs, ys


def corrcoef(xs: Sequence[float], ys: Sequence[float]) -> float:
    if not xs or not ys:
        return float("nan")
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if den_x == 0 or den_y == 0:
        return float("nan")
    return num / (den_x * den_y)


def correlation_entry(rows: Sequence[Dict[str, str]], key_x: str, key_y: str) -> Tuple[float, int]:
    xs, ys = finite_pair(rows, key_x, key_y)
    if not xs:
        return float("nan"), 0
    return corrcoef(xs, ys), len(xs)


def summarize(values: Iterable[float]) -> Dict[str, float]:
    arr = [v for v in values if math.isfinite(v)]
    if not arr:
        return {}
    arr.sort()
    n = len(arr)

    def percentile(p: float) -> float:
        if n == 1:
            return arr[0]
        idx = p * (n - 1)
        lo = math.floor(idx)
        hi = math.ceil(idx)
        if lo == hi:
            return arr[int(idx)]
        frac = idx - lo
        return arr[lo] * (1 - frac) + arr[hi] * frac

    return {
        "n": n,
        "mean": sum(arr) / n,
        "median": stats.median(arr),
        "std": stats.pstdev(arr) if n > 1 else 0.0,
        "min": arr[0],
        "max": arr[-1],
        "p05": percentile(0.05),
        "p95": percentile(0.95),
    }


def print_summary(title: str, rows: List[Dict[str, str]]) -> None:
    print(f"\n=== {title} ===")
    for key in ("CNR", "intensity", "sn_voxels", "total_volume"):
        stats_dict = summarize(to_float(r.get(key, "")) for r in rows)
        if stats_dict:
            print(f"{key}: {stats_dict}")
        else:
            print(f"{key}: no finite values")


def show_extremes(rows: List[Dict[str, str]], field: str, count: int = 5) -> None:
    valid = [r for r in rows if math.isfinite(to_float(r.get(field, "")))]
    sorted_rows = sorted(valid, key=lambda r: to_float(r[field]))
    print(f"Lowest {count} {field}:")
    for row in sorted_rows[:count]:
        print(
            f"  {row['sub']} -> {field}={row[field]} intensity={row.get('intensity')} sn_voxels={row.get('sn_voxels')} notes={row.get('notes','')}"
        )
    print(f"Highest {count} {field}:")
    for row in sorted_rows[-count:][::-1]:
        print(
            f"  {row['sub']} -> {field}={row[field]} intensity={row.get('intensity')} sn_voxels={row.get('sn_voxels')}"
        )


def list_missing(rows: List[Dict[str, str]]) -> List[str]:
    missing = []
    for row in rows:
        cnr = to_float(row.get("CNR", ""))
        inten = to_float(row.get("intensity", ""))
        if not math.isfinite(cnr) or not math.isfinite(inten):
            missing.append(row.get("sub", ""))
    return missing


def diff_summary(rows: List[Dict[str, str]], key_a: str, key_b: str) -> Dict[str, float]:
    diffs = []
    for row in rows:
        a = to_float(row.get(key_a, ""))
        b = to_float(row.get(key_b, ""))
        if math.isfinite(a) and math.isfinite(b):
            diffs.append(a - b)
    return summarize(diffs)


def mean_pair(rows: List[Dict[str, str]], key_a: str, key_b: str) -> Tuple[float, float]:
    vals_a, vals_b = [], []
    for row in rows:
        a = to_float(row.get(key_a, ""))
        b = to_float(row.get(key_b, ""))
        if math.isfinite(a):
            vals_a.append(a)
        if math.isfinite(b):
            vals_b.append(b)
    mean_a = sum(vals_a) / len(vals_a) if vals_a else float("nan")
    mean_b = sum(vals_b) / len(vals_b) if vals_b else float("nan")
    return mean_a, mean_b


def discover_pt_thresholds(rows: List[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    if not rows:
        return {}
    pattern = re.compile(r"^CNR_PT_thr(\d{3})(?:_(unthresholded|raw|corrected))?$")
    thresholds: Dict[str, Dict[str, str]] = {}
    for key in rows[0].keys():
        match = pattern.match(key)
        if not match:
            continue
        label = match.group(1)
        variant = match.group(2) or "active"
        thresholds.setdefault(label, {})[variant] = key
    return thresholds


def summarize_field(rows: List[Dict[str, str]], field: Optional[str]) -> str:
    if not field:
        return "n=0"
    stats_dict = summarize(to_float(r.get(field, "")) for r in rows)
    if not stats_dict:
        return "n=0"
    return f"n={stats_dict['n']} mean={stats_dict['mean']:.3f} median={stats_dict['median']:.3f}"


def voxel_volume_report(
    rows: List[Dict[str, str]],
    total_key: str,
    voxel_key: str,
    label: str,
    expected: float = VOXEL_VOLUME_MM3,
) -> None:
    ratios = []
    for row in rows:
        total = to_float(row.get(total_key, ""))
        vox = to_float(row.get(voxel_key, ""))
        if math.isfinite(total) and math.isfinite(vox) and vox > 0:
            ratios.append(total / vox)
    if not ratios:
        print(f"Voxel volume ({label}): no finite values")
        return
    stats_dict = summarize(ratios)
    delta = max(abs(r - expected) for r in ratios)
    status = "OK" if delta < 1e-3 else f"deviation max={delta:.4f}"
    print(f"Voxel volume ({label}): {stats_dict} expected={expected:.4f} -> {status}")


def check_bright_counts(
    rows: List[Dict[str, str]],
    bright_fields: Sequence[str],
    sn_voxel_key: str,
    label: str,
) -> None:
    violations = []
    for row in rows:
        sn_vox = to_float(row.get(sn_voxel_key, ""))
        if not math.isfinite(sn_vox) or sn_vox <= 0:
            continue
        for field in bright_fields:
            val = to_float(row.get(field, ""))
            if math.isfinite(val) and val > sn_vox + 1e-3:
                violations.append((row.get("sub", ""), field, val, sn_vox))
    if violations:
        print(f"[WARN] Bright voxel counts exceed SN voxels for {label}:")
        for sub, field, val, sn_vox in violations:
            print(f"  {sub}: {field}={val} > sn_voxels={sn_vox}")
    else:
        print(f"Bright voxel sanity ({label}): all counts within SN voxel totals")

def normalize_sub_id(raw_id: str) -> str:
    text = (raw_id or "").strip()
    if not text:
        return ""
    if text.startswith("HUMAN_EBR-"):
        return text.split("-", 1)[1]
    try:
        return str(int(float(text)))
    except ValueError:
        return text


def find_mask_path(mask_dir: Path, label: str, sub_id: str) -> Optional[Path]:
    for suffix in (".nii", ".nii.gz"):
        candidate = mask_dir / f"{label}_groupmask_in_{sub_id}{suffix}"
        if candidate.exists():
            return candidate
    return None


def make_mask_finder(mask_dir: Path, label: str):
    mask_dir = Path(mask_dir)

    def finder(sub_id: str) -> Optional[Path]:
        return find_mask_path(mask_dir, label, sub_id)

    return finder


def make_brainstem_finder(brainstem_dir: Path):
    brainstem_dir = Path(brainstem_dir)

    def finder(sub_id: str) -> Optional[Path]:
        for suffix in (".nii", ".nii.gz"):
            candidate = brainstem_dir / f"DC_mask_inTSE_{sub_id}{suffix}"
            if candidate.exists():
                return candidate
        return None

    return finder


def parse_run_specs(main_path: str, extra_args: Sequence[str]) -> List[Tuple[str, str]]:
    runs: List[Tuple[str, str]] = []
    runs.append((Path(main_path).stem or "main", main_path))
    for spec in extra_args:
        if "=" in spec:
            label, path = spec.split("=", 1)
            label = label.strip() or Path(path).stem
            runs.append((label, path.strip()))
        else:
            runs.append((Path(spec).stem, spec))
    return runs


def normalize_watch_list(watch_ids: Sequence[str]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for token in watch_ids:
        norm = normalize_sub_id(token)
        if norm:
            mapping[norm] = token
    return mapping


def watchlist_missing(rows: List[Dict[str, str]], watch_map: Dict[str, str]) -> List[str]:
    if not watch_map:
        return []
    missing_norm = {normalize_sub_id(sub) for sub in list_missing(rows)}
    still_missing = [orig for norm, orig in watch_map.items() if norm in missing_norm]
    return still_missing


def gather_nigrosome_intensity_specs(rows: Sequence[Dict[str, str]]) -> Dict[str, List[Tuple[str, str, str]]]:
    """Return mapping of nigrosome prefixes to (label, cnr_key, intensity_key)."""

    if not rows:
        return {}

    header_keys = list(rows[0].keys())
    header_set = set(header_keys)
    specs: Dict[str, List[Tuple[str, str, str]]] = {}
    for key in header_keys:
        if "_CNR_" not in key or key.startswith("CNR_"):
            continue
        prefix, remainder = key.split("_CNR_", 1)
        if not prefix:
            continue
        cnr_parts = remainder.split("_", 1)
        contrast = cnr_parts[0]
        suffix = ""
        if len(cnr_parts) > 1:
            suffix = f"_{cnr_parts[1]}"
        intensity_key = f"{prefix}_intensity{suffix}"
        if intensity_key not in header_set:
            continue
        suffix_label = suffix.lstrip("_")
        label = contrast if not suffix_label else f"{contrast} {suffix_label.replace('_', ' ')}"
        specs.setdefault(prefix, []).append((label, key, intensity_key))

    for entries in specs.values():
        entries.sort(key=lambda entry: entry[0])
    return dict(sorted(specs.items()))


def report_run(
    label: str,
    path: str,
    rows: List[Dict[str, str]],
    expect: int,
    watch_map: Dict[str, str],
    cleanup_rows: Optional[List[Dict[str, str]]] = None,
) -> None:
    print(f"\n=== Run {label}: {Path(path).name} ===")
    print(f"Loaded {len(rows)} rows from {path}")
    if expect:
        delta = len(rows) - expect
        status = "OK" if delta == 0 else f"Δ={delta}"
        print(f"Expected {expect} rows -> {status}")

    print("\nCNR vs intensity correlations:")
    corr_specs = [
        ("CP", "CNR_CP", "intensity"),
        ("PT", "CNR_PT", "intensity"),
        ("BS", "CNR_BS", "intensity"),
    ]
    for desc, x_key, y_key in corr_specs:
        corr, count = correlation_entry(rows, x_key, y_key)
        print(f"  {desc}: corr={corr:.3f} (n={count})")

    print("CNR vs intensity correlations (unthresholded masks):")
    corr_specs_unthr = [
        ("CP", "CNR_CP_unthresholded", "intensity_unthresholded"),
        ("PT", "CNR_PT_unthresholded", "intensity_unthresholded"),
        ("BS", "CNR_BS_unthresholded", "intensity_unthresholded"),
    ]
    for desc, x_key, y_key in corr_specs_unthr:
        corr, count = correlation_entry(rows, x_key, y_key)
        print(f"  {desc}: corr={corr:.3f} (n={count})")

    print("Global CNR vs intensity (raw vs corrected):")
    canonical_specs = [
        ("raw", "CNR_raw", "intensity_raw"),
        ("corrected", "CNR_corrected", "intensity_corrected"),
    ]
    for desc, x_key, y_key in canonical_specs:
        corr, count = correlation_entry(rows, x_key, y_key)
        print(f"  {desc}: corr={corr:.3f} (n={count})")

    nigrosome_specs = gather_nigrosome_intensity_specs(rows)
    if nigrosome_specs:
        print("Nigrosome CNR vs intensity correlations:")
        for prefix, combos in nigrosome_specs.items():
            print(f"  {prefix}:")
            printed = False
            for label, cnr_key, inten_key in combos:
                corr, count = correlation_entry(rows, cnr_key, inten_key)
                if count == 0:
                    print(f"    {label}: corr=nan (n=0)")
                    continue
                printed = True
                print(f"    {label}: corr={corr:.3f} (n={count})")
            if not printed:
                print("    no finite CNR/intensity pairs found")

    print("Top/bottom voxel correlations:")
    edge_specs = [
        ("CP top10", "CNR_CP", "intensity_maximum"),
        ("CP bottom10", "CNR_CP", "intensity_minimum"),
        ("CP median", "CNR_CP", "intensity_median"),
    ]
    for desc, x_key, y_key in edge_specs:
        corr, count = correlation_entry(rows, x_key, y_key)
        print(f"  {desc}: corr={corr:.3f} (n={count})")
    edge_specs_unthr = [
        ("CP top10 (unthr)", "CNR_CP_unthresholded", "intensity_maximum_unthresholded"),
        ("CP bottom10 (unthr)", "CNR_CP_unthresholded", "intensity_minimum_unthresholded"),
        (
            "CP median (unthr)",
            "CNR_CP_unthresholded",
            "intensity_median_unthresholded",
        ),
    ]
    for desc, x_key, y_key in edge_specs_unthr:
        corr, count = correlation_entry(rows, x_key, y_key)
        print(f"  {desc}: corr={corr:.3f} (n={count})")

    print_summary("Main table", rows)
    show_extremes(rows, "CNR")

    notes_counter = Counter()
    for row in rows:
        note_field = row.get("notes", "") or ""
        for token in note_field.split(","):
            token = token.strip()
            if token:
                notes_counter[token] += 1
    print("Notes counts:", notes_counter)

    missing = list_missing(rows)
    print(f"Rows with NaN metrics ({len(missing)}): {missing}")
    if watch_map:
        watch_missing = watchlist_missing(rows, watch_map)
        print(f"Watch list still missing: {watch_missing if watch_missing else 'none'}")

    if cleanup_rows:
        print_summary("Cleanup table", cleanup_rows)

    for label_name, keys in (
        ("volume", ("volume_left", "volume_right")),
        ("intensity", ("intensity_left", "intensity_right")),
        ("max", ("max_left", "max_right")),
    ):
        mean_l, mean_r = mean_pair(rows, keys[0], keys[1])
        diff_stats = diff_summary(rows, keys[0], keys[1])
        print(
            f"Left vs Right {label_name}: mean_left={mean_l:.3f}, mean_right={mean_r:.3f}, diff stats={diff_stats}"
        )

    for label_name, keys in (
        ("volume", ("volume_rostral", "volume_caudal")),
        ("intensity", ("intensity_rostral", "intensity_caudal")),
        ("max", ("max_rostral", "max_caudal")),
    ):
        mean_r, mean_c = mean_pair(rows, keys[0], keys[1])
        diff_stats = diff_summary(rows, keys[0], keys[1])
        print(
            f"Rostral vs Caudal {label_name}: mean_rostral={mean_r:.3f}, mean_caudal={mean_c:.3f}, diff stats={diff_stats}"
        )

    voxel_volume_report(rows, "total_volume", "sn_voxels", "thresholded")
    voxel_volume_report(rows, "total_volume_unthresholded", "sn_voxels_unthresholded", "unthresholded")

    bright_fields = [
        "bright_voxels_cp",
        "bright_voxels_pt",
        "bright_voxels_bs",
        "bright_voxels_left",
        "bright_voxels_right",
        "bright_voxels_rostral",
        "bright_voxels_caudal",
    ]
    bright_fields_unthr = [f"{field}_unthresholded" for field in bright_fields]
    check_bright_counts(rows, bright_fields, "sn_voxels", "thresholded")
    check_bright_counts(rows, bright_fields_unthr, "sn_voxels_unthresholded", "unthresholded")

    pt_thresholds = discover_pt_thresholds(rows)
    if pt_thresholds:
        print("\nPT mask compare thresholds:")
        for label in sorted(pt_thresholds.keys()):
            spec = pt_thresholds[label]
            try:
                prob = int(label) / 100.0
            except ValueError:
                prob = float("nan")
            prob_desc = f"≥ {prob:.2f}" if math.isfinite(prob) else label
            print(f"  thr{label} ({prob_desc}):")
            for variant, desc, inten_key in (
                ("active", "prob", "intensity"),
                ("unthresholded", "unthr", "intensity_unthresholded"),
                ("raw", "raw", "intensity_raw"),
                ("corrected", "corr", "intensity_corrected"),
            ):
                stats_text = summarize_field(rows, spec.get(variant))
                print(f"    {desc}: {stats_text}")
                field = spec.get(variant)
                if field:
                    corr, count = correlation_entry(rows, field, inten_key)
                    print(f"      corr vs {inten_key}: corr={corr:.3f} (n={count})")


def run_afni_overlap(sn_mask: Path, ctrl_mask: Path, sn_threshold: float) -> Optional[int]:
    with tempfile.TemporaryDirectory() as tmpdir:
        overlap_path = Path(tmpdir) / "overlap.nii.gz"
        expr = "step(a)*step(b)"
        if sn_threshold is not None:
            expr = f"step(a-{sn_threshold})*step(b)"
        calc_cmd = [
            "3dcalc",
            "-a",
            str(sn_mask),
            "-b",
            str(ctrl_mask),
            "-expr",
            expr,
            "-prefix",
            str(overlap_path),
        ]
        try:
            subprocess.run(calc_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except (OSError, subprocess.CalledProcessError) as exc:
            print(f"  !! Failed to run {' '.join(calc_cmd)} ({exc})")
            return None
        stat_cmd = ["3dBrickStat", "-non-zero", "-count", str(overlap_path)]
        try:
            result = subprocess.run(stat_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except (OSError, subprocess.CalledProcessError) as exc:
            print(f"  !! Failed to run {' '.join(stat_cmd)} ({exc})")
            return None
        try:
            return int(float(result.stdout.decode().strip()))
        except ValueError:
            return None


def check_mask_overlap(
    rows: List[Dict[str, str]],
    sn_mask_dir: Path,
    background_finder,
    label: str,
    sn_threshold: float,
) -> None:
    if not sn_mask_dir.is_dir():
        print(f"SN mask directory {sn_mask_dir} not found; skipping {label} overlap check.")
        return

    overlaps = []
    missing = []
    for row in rows:
        sub_raw = row.get("sub", "")
        sub_id = normalize_sub_id(sub_raw)
        if not sub_id:
            continue
        sn_path = find_mask_path(sn_mask_dir, "sn", sub_id)
        bg_path = background_finder(sub_id)
        if sn_path is None or bg_path is None:
            missing.append(sub_raw)
            continue
        overlap_vox = run_afni_overlap(sn_path, bg_path, sn_threshold)
        if overlap_vox is None:
            missing.append(sub_raw)
            continue
        if overlap_vox > 0:
            overlaps.append((sub_raw, overlap_vox))

    prefix = f"SN vs {label}"
    if missing:
        print(f"{prefix}: missing mask files or overlap calc failed for {sorted(set(missing))}")
    if overlaps:
        print(f"{prefix}: overlap detected (voxels)")
        for sub, vox in overlaps:
            print(f"  {sub}: {vox} voxels overlap")
    else:
        print(f"{prefix}: no intersections detected.")


def run_overlap_suite(rows: List[Dict[str, str]], args: argparse.Namespace) -> None:
    sn_dir = Path(args.mask_dir)
    specs = [
        ("SN-control", make_mask_finder(sn_dir, "snctrl")),
        ("PT-control", make_mask_finder(sn_dir, "ptctrl")),
    ]
    brainstem_dir = Path(args.brainstem_dir)
    if brainstem_dir.is_dir():
        specs.append(("brainstem/DC", make_brainstem_finder(brainstem_dir)))
    else:
        print(f"Brainstem directory {brainstem_dir} not found; skipping DC overlap check.")

    for label, finder in specs:
        check_mask_overlap(rows, sn_dir, finder, label, args.sn_overlap_threshold)


def collect_recent_tests(
    analysis_dir: str,
    pattern: str,
    limit: int,
    existing_paths: Sequence[str],
) -> List[Tuple[str, str]]:
    if limit <= 0:
        return []
    base_path = Path(analysis_dir)
    if not base_path.is_dir():
        print(f"Analysis dir {analysis_dir} not found; skipping recent-test autodiscovery.")
        return []

    existing: Set[str] = set()
    for raw_path in existing_paths:
        try:
            existing.add(str(Path(raw_path).resolve()))
        except OSError:
            existing.add(str(Path(raw_path).absolute()))

    candidates: List[Tuple[float, Path]] = []
    for match in base_path.glob(pattern):
        try:
            mtime = match.stat().st_mtime
        except OSError:
            continue
        candidates.append((mtime, match))
    candidates.sort(key=lambda item: item[0], reverse=True)

    selected: List[Tuple[str, str]] = []
    for _, candidate in candidates:
        if len(selected) >= limit:
            break
        resolved = str(candidate.resolve())
        if resolved in existing:
            continue
        selected.append((candidate.stem or candidate.name, str(candidate)))
    return selected


def summarize_pt_variants(root: str, pattern: str) -> None:
    root_path = Path(root)
    if not root_path.is_dir():
        print(f"PT variant root {root} not found; skipping summary.")
        return

    matches = sorted(root_path.rglob(pattern))
    if not matches:
        print(f"PT variant summary: no files matched pattern {pattern} under {root}.")
        return

    summary_rows = []
    for path in matches:
        try:
            rows = load_rows(str(path))
        except OSError as exc:
            print(f"[WARN] Unable to read {path}: {exc}")
            continue
        corr, valid = correlation_entry(rows, "CNR_PT", "intensity")
        total = len(rows)
        try:
            label = str(path.relative_to(root_path))
        except ValueError:
            label = str(path)
        summary_rows.append((label, str(path), corr, valid, total))

    if not summary_rows:
        print("PT variant summary: no readable tables found.")
        return

    summary_rows.sort(
        key=lambda item: (
            math.isnan(item[2]),
            -item[2] if math.isfinite(item[2]) else float("-inf"),
        )
    )
    print("\n=== PT variant CNR_PT vs intensity correlations ===")
    print("label | corr | n_valid | total | path")
    best_entry = None
    for label, path, corr, valid, total in summary_rows:
        corr_text = f"{corr:.3f}" if math.isfinite(corr) else "NaN"
        if best_entry is None and math.isfinite(corr):
            best_entry = (label, corr, valid, total)
        print(f"{label} | {corr_text} | {valid} | {total} | {path}")
    if best_entry:
        label, corr, valid, total = best_entry
        print(
            f"\nBest PT variant: {label} corr={corr:.3f} (n={valid}/{total})"
        )


def main() -> None:
    args = parse_args()
    runs = parse_run_specs(args.main, args.extra)
    recent_runs = collect_recent_tests(
        args.analysis_dir,
        args.recent_tests_pattern,
        args.recent_tests,
        [path for _, path in runs],
    )
    runs.extend(recent_runs)

    cleanup_rows = load_rows(args.cleanup) if args.cleanup else []
    watch_map = normalize_watch_list(args.watch)

    for idx, (label, path) in enumerate(runs):
        rows = load_rows(path)
        expect = args.expect if idx == 0 else 0
        cleanup = cleanup_rows if idx == 0 else None
        report_run(label, path, rows, expect, watch_map, cleanup)
        if idx == 0 and args.check_overlap:
            print("\n=== Mask overlap checks ===")
            run_overlap_suite(rows, args)

    if args.pt_variant_root:
        summarize_pt_variants(args.pt_variant_root, args.pt_variant_pattern)


if __name__ == "__main__":
    main()
