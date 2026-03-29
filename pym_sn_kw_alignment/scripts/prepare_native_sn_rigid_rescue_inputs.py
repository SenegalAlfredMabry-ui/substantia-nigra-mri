#!/usr/bin/env python3
"""Stage starting SN masks for native-space rigid rescue.

Priority:
1. hard16 neg20 wins currently measured as better than neg40
2. hard16 routing selected closer approach
3. non66 partial closer pilot outputs with > min improvement
4. frontier66 best measured closer approach
5. fallback to v10 base mask
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def link_or_copy(src: Path, dst: Path) -> None:
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src)


def find_mask(base: Path, subject: str) -> Path | None:
    for ext in (".nii.gz", ".nii"):
        p = base / f"sn_groupmask_in_{subject}{ext}"
        if p.exists():
            return p
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", type=Path, required=True)
    ap.add_argument("--output-root", type=Path, required=True)
    ap.add_argument(
        "--fallback-root",
        type=Path,
        default=Path("${SNHEART_ROOT:-/path/to/SNHEART}/MRI_data/TSE/probabilistic_masks_v10/base"),
    )
    ap.add_argument("--non66-min-improvement", type=float, default=0.1)
    args = ap.parse_args()

    run_dir = args.run_dir
    out_root = args.output_root
    out_root.mkdir(parents=True, exist_ok=True)
    staged_dir = out_root / "staged_masks"
    staged_dir.mkdir(exist_ok=True)
    manual_dir = out_root / "manual_masks"
    manual_dir.mkdir(exist_ok=True)

    subjects: list[str] = []
    group_for_subject: dict[str, str] = {}
    for rel in ("highhat/highhat_rawkw_subjects.txt", "crash/crash_rawkw_subjects.txt"):
        group = rel.split("/", 1)[0]
        path = run_dir / rel
        with path.open(encoding="utf-8") as f:
            for line in f:
                tok = line.strip().split()
                if tok:
                    sid = tok[0][-3:]
                    if sid not in subjects:
                        subjects.append(sid)
                    group_for_subject[sid] = group

    chosen: dict[str, tuple[str, Path]] = {}

    neg20_cmp = run_dir / "hard16_neg20_vs_neg40_partial_latest.tsv"
    if neg20_cmp.exists():
        for row in read_tsv(neg20_cmp):
            if row["winner_so_far"] == "neg20_better" and row["neg20_direction"] == "closer":
                sub = row["subject"]
                src = find_mask(run_dir / "multiopt_outputs" / "hard16_neg20_targeted" / "probabilistic_masks_v10" / "base", sub)
                if src:
                    chosen[sub] = ("hard16_neg20_targeted", src)

    hard16 = run_dir / "hard16_subject_routing_latest.tsv"
    if hard16.exists():
        for row in read_tsv(hard16):
            sub = row["subject"]
            if sub in chosen:
                continue
            if row.get("selected_direction") != "closer":
                continue
            approach = row["selected_approach"]
            src = find_mask(run_dir / "multiopt_outputs" / approach / "probabilistic_masks_v10" / "base", sub)
            if src:
                chosen[sub] = (approach, src)

    non66 = run_dir / "non66_partial_direction_snapshot_latest.tsv"
    if non66.exists():
        for row in read_tsv(non66):
            sub = row["subject"]
            if sub in chosen:
                continue
            if row.get("direction") != "closer":
                continue
            try:
                imp = float(row["improvement_mm"])
            except Exception:
                continue
            if imp < args.non66_min_improvement:
                continue
            approach = row["approach"]
            src = find_mask(run_dir / "multiopt_outputs" / approach / "probabilistic_masks_v10" / "base", sub)
            if src:
                chosen[sub] = (approach, src)

    frontier = run_dir / "multiopt_subject_strategy_latest.tsv"
    if frontier.exists():
        for row in read_tsv(frontier):
            sub = row["subject"]
            if sub in chosen:
                continue
            if row.get("best_direction") != "closer":
                continue
            approach = row["best_approach"]
            src = find_mask(run_dir / "multiopt_outputs" / approach / "probabilistic_masks_v10" / "base", sub)
            if src:
                chosen[sub] = (approach, src)

    manifest = out_root / "starting_mask_manifest.tsv"
    with manifest.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["subject", "group", "source_type", "source_path", "staged_path", "manual_path"])
        for sub in sorted(subjects, key=int):
            manual_src = None
            group = group_for_subject.get(sub, "")
            if group:
                manual_src = find_mask(run_dir / group / "manual_raw_kw", sub)
            manual_dst_path = ""
            if manual_src is not None:
                manual_suffix = ".nii.gz" if manual_src.name.endswith(".nii.gz") else ".nii"
                manual_dst = manual_dir / f"sn_groupmask_in_{sub}{manual_suffix}"
                link_or_copy(manual_src, manual_dst)
                manual_dst_path = str(manual_dst)
            if sub in chosen:
                source_type, src = chosen[sub]
            else:
                src = find_mask(args.fallback_root, sub)
                source_type = "v10_base"
            if src is None:
                w.writerow([sub, group, "missing", "", "", manual_dst_path])
                continue
            suffix = ".nii.gz" if src.name.endswith(".nii.gz") else ".nii"
            dst = staged_dir / f"sn_groupmask_in_{sub}{suffix}"
            link_or_copy(src, dst)
            w.writerow([sub, group, source_type, str(src), str(dst), manual_dst_path])

    print(f"WROTE {manifest}")


if __name__ == "__main__":
    main()
