#!/usr/bin/env python3
"""Summarize multi-approach progress with participant-level best strategy.

Reads per-approach `*_current_now.tsv` files produced by snapshot script and
compares against group pre tables to compute improvement per subject.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List


def read_tsv(path: Path) -> List[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-dir", type=Path, required=True)
    ap.add_argument(
        "--out-tsv",
        type=Path,
        help="Participant-level best-strategy table (default: <run>/multiopt_subject_strategy_latest.tsv)",
    )
    ap.add_argument(
        "--out-summary",
        type=Path,
        help="Summary text (default: <run>/multiopt_subject_strategy_summary_latest.txt)",
    )
    args = ap.parse_args()

    run_dir: Path = args.run_dir
    multi_dir = run_dir / "multiopt_outputs"
    if not multi_dir.exists():
        raise FileNotFoundError(f"Missing {multi_dir}")

    out_tsv = args.out_tsv or (run_dir / "multiopt_subject_strategy_latest.tsv")
    out_summary = args.out_summary or (run_dir / "multiopt_subject_strategy_summary_latest.txt")

    pre_by_group: Dict[str, Dict[str, float]] = {}
    roster_by_group: Dict[str, List[str]] = {}
    for group in ("highhat", "crash"):
        pre_path = run_dir / group / f"{group}_rawkw_vs_v10_pre.tsv"
        roster_path = run_dir / "multiopt_outputs" / "strict_full" / f"{group}_subjects.txt"
        if not pre_path.exists():
            raise FileNotFoundError(pre_path)
        if not roster_path.exists():
            raise FileNotFoundError(roster_path)
        pre_map: Dict[str, float] = {}
        for r in read_tsv(pre_path):
            s = r.get("subject", "")
            try:
                d = float(r.get("delta_mag_mm", "nan"))
            except Exception:
                continue
            if s and math.isfinite(d):
                pre_map[s] = d
        pre_by_group[group] = pre_map
        roster_by_group[group] = [x.strip() for x in roster_path.read_text(encoding="utf-8").splitlines() if x.strip()]

    approaches = sorted([p.name for p in multi_dir.iterdir() if p.is_dir()])
    if not approaches:
        raise RuntimeError(f"No approach dirs under {multi_dir}")

    # subject -> approach -> improvement
    improvements: Dict[str, Dict[str, float]] = {}
    group_for_subject: Dict[str, str] = {}
    for group, subs in roster_by_group.items():
        for s in subs:
            group_for_subject[s] = group
            improvements.setdefault(s, {})

    for approach in approaches:
        adir = multi_dir / approach
        for group in ("highhat", "crash"):
            cur = adir / f"{group}_current_now.tsv"
            if not cur.exists():
                continue
            pre_map = pre_by_group[group]
            for r in read_tsv(cur):
                s = r.get("subject", "")
                if s not in pre_map:
                    continue
                try:
                    post_d = float(r.get("delta_mag_mm", "nan"))
                except Exception:
                    continue
                if not math.isfinite(post_d):
                    continue
                imp = pre_map[s] - post_d
                improvements.setdefault(s, {})[approach] = imp
                group_for_subject.setdefault(s, group)

    rows: List[dict] = []
    closer_total = farther_total = pending_total = 0
    best_counts: Dict[str, int] = {a: 0 for a in approaches}

    for s in sorted(improvements.keys(), key=lambda x: int(x)):
        group = group_for_subject.get(s, "unknown")
        app_map = improvements.get(s, {})
        if not app_map:
            pending_total += 1
            rows.append(
                {
                    "subject": s,
                    "group": group,
                    "status": "pending",
                    "best_approach": "",
                    "best_improvement_mm": "",
                    "best_direction": "",
                    "measured_approaches": "0",
                    "approach_improvements_json": "{}",
                }
            )
            continue

        best_approach = max(app_map.keys(), key=lambda a: app_map[a])
        best_imp = app_map[best_approach]
        if best_imp > 1e-6:
            direction = "closer"
            closer_total += 1
        elif best_imp < -1e-6:
            direction = "farther"
            farther_total += 1
        else:
            direction = "neutral"
        best_counts[best_approach] = best_counts.get(best_approach, 0) + 1

        rows.append(
            {
                "subject": s,
                "group": group,
                "status": "measured",
                "best_approach": best_approach,
                "best_improvement_mm": f"{best_imp:.6f}",
                "best_direction": direction,
                "measured_approaches": str(len(app_map)),
                "approach_improvements_json": json.dumps(
                    {k: round(v, 6) for k, v in sorted(app_map.items())},
                    sort_keys=True,
                ),
            }
        )

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "subject",
                "group",
                "status",
                "best_approach",
                "best_improvement_mm",
                "best_direction",
                "measured_approaches",
                "approach_improvements_json",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)

    total = len(rows)
    measured = total - pending_total
    with out_summary.open("w", encoding="utf-8") as f:
        f.write(f"run_dir\t{run_dir}\n")
        f.write(f"approaches\t{','.join(approaches)}\n")
        f.write(
            f"subjects_total\t{total}\tmeasured\t{measured}\tpending\t{pending_total}"
            f"\tbest_closer\t{closer_total}\tbest_farther\t{farther_total}\n"
        )
        for a in approaches:
            f.write(f"best_count\t{a}\t{best_counts.get(a, 0)}\n")

    print(f"WROTE {out_tsv}")
    print(f"WROTE {out_summary}")
    print(
        f"subjects={total} measured={measured} pending={pending_total} "
        f"best_closer={closer_total} best_farther={farther_total}"
    )


if __name__ == "__main__":
    main()
