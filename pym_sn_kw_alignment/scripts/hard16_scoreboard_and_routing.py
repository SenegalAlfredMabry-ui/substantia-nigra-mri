#!/usr/bin/env python3
"""Build live hard16 rescue scoreboard + routing artifacts.

Outputs:
  - per-subject wide scoreboard with baseline + each approach delta/improvement
  - concise summary with acceptance-rule counts
  - routing TSV/JSON (best-so-far per subject)
  - selected offset JSONs (ready-only and best-so-far)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple


def read_tsv(path: Path) -> List[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_subject_list(path: Path) -> List[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def f6(val: float | None) -> str:
    if val is None or not math.isfinite(val):
        return ""
    return f"{val:.6f}"


def direction(improvement_mm: float) -> str:
    if improvement_mm > 1e-6:
        return "closer"
    if improvement_mm < -1e-6:
        return "farther"
    return "neutral"


def load_pre_maps(run_dir: Path) -> Tuple[Dict[str, float], Dict[str, str]]:
    pre: Dict[str, float] = {}
    group_by_sub: Dict[str, str] = {}
    for group in ("highhat", "crash"):
        pre_path = run_dir / group / f"{group}_rawkw_vs_v10_pre.tsv"
        if not pre_path.exists():
            raise FileNotFoundError(f"Missing pre table: {pre_path}")
        for row in read_tsv(pre_path):
            sub = row.get("subject", "")
            if not sub:
                continue
            try:
                val = float(row.get("delta_mag_mm", "nan"))
            except Exception:
                continue
            if math.isfinite(val):
                pre[sub] = val
                group_by_sub[sub] = group
    return pre, group_by_sub


def discover_approaches(multi_dir: Path, approach_glob: str, include: List[str]) -> List[str]:
    found = sorted(p.name for p in multi_dir.glob(approach_glob) if p.is_dir())
    if include:
        for name in include:
            p = multi_dir / name
            if p.is_dir() and name not in found:
                found.append(name)
    return sorted(set(found))


def load_post_by_approach(run_dir: Path, approaches: List[str]) -> Dict[str, Dict[str, float]]:
    out: Dict[str, Dict[str, float]] = {}
    multi_dir = run_dir / "multiopt_outputs"
    for approach in approaches:
        vals: Dict[str, float] = {}
        for group in ("highhat", "crash"):
            cur_path = multi_dir / approach / f"{group}_current_now.tsv"
            if not cur_path.exists():
                continue
            for row in read_tsv(cur_path):
                sub = row.get("subject", "")
                if not sub:
                    continue
                try:
                    d = float(row.get("delta_mag_mm", "nan"))
                except Exception:
                    continue
                if math.isfinite(d):
                    vals[sub] = d
        out[approach] = vals
    return out


def load_offset_jsons(run_dir: Path, approaches: List[str]) -> Dict[str, Dict[str, dict]]:
    offset_dir = run_dir / "multiopt_offsets"
    out: Dict[str, Dict[str, dict]] = {}
    for approach in approaches:
        p = offset_dir / f"{approach}.json"
        if not p.exists():
            continue
        try:
            out[approach] = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-dir", type=Path, required=True)
    ap.add_argument(
        "--subjects-file",
        type=Path,
        default=None,
        help="Subjects to score (default: <run>/hard16_best_farther_subjects.txt)",
    )
    ap.add_argument(
        "--approach-glob",
        default="hard16_*",
        help="Approach dir glob under multiopt_outputs (default: hard16_*)",
    )
    ap.add_argument(
        "--include-approaches",
        default="strict_full,strict_gain70,strict_gain40,data_driven_v1",
        help="Comma list of extra approaches to include if present",
    )
    ap.add_argument("--min-improvement-mm", type=float, default=1.0)
    ap.add_argument(
        "--require-closer",
        type=int,
        default=1,
        help="1=require best_direction==closer for pass, 0=threshold only",
    )
    ap.add_argument(
        "--scoreboard-tsv",
        type=Path,
        default=None,
        help="Default: <run>/hard16_rescue_scoreboard_latest.tsv",
    )
    ap.add_argument(
        "--summary-txt",
        type=Path,
        default=None,
        help="Default: <run>/hard16_rescue_summary_latest.txt",
    )
    ap.add_argument(
        "--routing-tsv",
        type=Path,
        default=None,
        help="Default: <run>/hard16_subject_routing_latest.tsv",
    )
    ap.add_argument(
        "--routing-json",
        type=Path,
        default=None,
        help="Default: <run>/hard16_subject_routing_latest.json",
    )
    ap.add_argument(
        "--offsets-ready-json",
        type=Path,
        default=None,
        help="Default: <run>/hard16_selected_offsets_ready.json",
    )
    ap.add_argument(
        "--offsets-best-json",
        type=Path,
        default=None,
        help="Default: <run>/hard16_selected_offsets_best_so_far.json",
    )
    args = ap.parse_args()

    run_dir = args.run_dir
    if not run_dir.exists():
        raise FileNotFoundError(f"Missing run dir: {run_dir}")

    subjects_file = args.subjects_file or (run_dir / "hard16_best_farther_subjects.txt")
    if not subjects_file.exists():
        raise FileNotFoundError(f"Missing subjects file: {subjects_file}")

    scoreboard_tsv = args.scoreboard_tsv or (run_dir / "hard16_rescue_scoreboard_latest.tsv")
    summary_txt = args.summary_txt or (run_dir / "hard16_rescue_summary_latest.txt")
    routing_tsv = args.routing_tsv or (run_dir / "hard16_subject_routing_latest.tsv")
    routing_json = args.routing_json or (run_dir / "hard16_subject_routing_latest.json")
    offsets_ready_json = args.offsets_ready_json or (run_dir / "hard16_selected_offsets_ready.json")
    offsets_best_json = args.offsets_best_json or (run_dir / "hard16_selected_offsets_best_so_far.json")

    pre_map, group_map = load_pre_maps(run_dir)
    subjects = read_subject_list(subjects_file)

    include = [x.strip() for x in args.include_approaches.split(",") if x.strip()]
    approaches = discover_approaches(run_dir / "multiopt_outputs", args.approach_glob, include)
    post_by_app = load_post_by_approach(run_dir, approaches)
    offset_by_app = load_offset_jsons(run_dir, approaches)

    wide_rows: List[dict] = []
    routing_rows: List[dict] = []
    routing_payload: Dict[str, dict] = {}
    offsets_ready: Dict[str, dict] = {}
    offsets_best: Dict[str, dict] = {}

    measured = passed = failed = pending = 0
    require_closer = args.require_closer == 1

    for sub in sorted(subjects, key=lambda x: int(x)):
        baseline = pre_map.get(sub)
        group = group_map.get(sub, "")
        row: Dict[str, str] = {"subject": sub, "group": group, "baseline_delta_mm": f6(baseline)}

        best_app = ""
        best_post = None
        best_imp = None
        measured_apps = 0
        app_improvements: Dict[str, float] = {}

        for app in approaches:
            post = post_by_app.get(app, {}).get(sub)
            imp = None
            dirn = ""
            if post is not None and baseline is not None:
                imp = baseline - post
                dirn = direction(imp)
                measured_apps += 1
                app_improvements[app] = imp
                if best_imp is None or imp > best_imp:
                    best_imp = imp
                    best_post = post
                    best_app = app
            row[f"{app}_post_delta_mm"] = f6(post)
            row[f"{app}_improvement_mm"] = f6(imp)
            row[f"{app}_direction"] = dirn

        if best_imp is None:
            status = "pending"
            pending += 1
            best_dir = ""
            best_pct = None
            pass_flag = ""
            decision = "pending_more_results"
        else:
            status = "measured"
            measured += 1
            best_dir = direction(best_imp)
            best_pct = (best_imp / baseline * 100.0) if (baseline is not None and baseline != 0) else None
            pass_rule = best_imp >= args.min_improvement_mm
            if require_closer:
                pass_rule = pass_rule and best_dir == "closer"
            if pass_rule:
                passed += 1
                pass_flag = "pass"
                decision = "ready"
            else:
                failed += 1
                pass_flag = "fail"
                decision = "needs_more_search"

        row["status"] = status
        row["measured_approaches"] = str(measured_apps)
        row["best_approach"] = best_app
        row["best_post_delta_mm"] = f6(best_post)
        row["best_improvement_mm"] = f6(best_imp)
        row["best_percent_reduction"] = f6(best_pct)
        row["best_direction"] = best_dir
        row["pass_closer_ge_threshold"] = pass_flag
        row["approach_improvements_json"] = json.dumps(
            {k: round(v, 6) for k, v in sorted(app_improvements.items())},
            sort_keys=True,
        )
        wide_rows.append(row)

        r = {
            "subject": sub,
            "group": group,
            "status": status,
            "selected_approach": best_app,
            "selected_improvement_mm": f6(best_imp),
            "selected_direction": best_dir,
            "selected_post_delta_mm": f6(best_post),
            "selected_percent_reduction": f6(best_pct),
            "decision": decision,
            "pass_closer_ge_threshold": pass_flag,
        }
        routing_rows.append(r)
        routing_payload[sub] = r

        if best_app and best_app in offset_by_app and sub in offset_by_app[best_app]:
            offsets_best[sub] = offset_by_app[best_app][sub]
            if pass_flag == "pass":
                offsets_ready[sub] = offset_by_app[best_app][sub]

    fixed_cols = [
        "subject",
        "group",
        "baseline_delta_mm",
        "status",
        "measured_approaches",
        "best_approach",
        "best_post_delta_mm",
        "best_improvement_mm",
        "best_percent_reduction",
        "best_direction",
        "pass_closer_ge_threshold",
    ]
    dynamic_cols: List[str] = []
    for app in approaches:
        dynamic_cols.extend(
            [f"{app}_post_delta_mm", f"{app}_improvement_mm", f"{app}_direction"]
        )
    cols = fixed_cols + dynamic_cols + ["approach_improvements_json"]

    with scoreboard_tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for row in wide_rows:
            w.writerow(row)

    route_cols = [
        "subject",
        "group",
        "status",
        "selected_approach",
        "selected_improvement_mm",
        "selected_direction",
        "selected_post_delta_mm",
        "selected_percent_reduction",
        "pass_closer_ge_threshold",
        "decision",
    ]
    with routing_tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=route_cols, delimiter="\t")
        w.writeheader()
        for row in routing_rows:
            w.writerow(row)

    routing_json.write_text(json.dumps(routing_payload, indent=2, sort_keys=True), encoding="utf-8")
    offsets_ready_json.write_text(json.dumps(offsets_ready, indent=2, sort_keys=True), encoding="utf-8")
    offsets_best_json.write_text(json.dumps(offsets_best, indent=2, sort_keys=True), encoding="utf-8")

    summary_lines = [
        f"run_dir\t{run_dir}",
        f"subjects_file\t{subjects_file}",
        f"approaches\t{','.join(approaches)}",
        f"acceptance_rule\trequire_closer={1 if require_closer else 0};min_improvement_mm={args.min_improvement_mm:.3f}",
        f"subjects_total\t{len(subjects)}\tmeasured\t{measured}\tpending\t{pending}\tpass\t{passed}\tfail\t{failed}",
        f"scoreboard_tsv\t{scoreboard_tsv}",
        f"routing_tsv\t{routing_tsv}",
        f"routing_json\t{routing_json}",
        f"offsets_ready_json\t{offsets_ready_json}\tcount\t{len(offsets_ready)}",
        f"offsets_best_json\t{offsets_best_json}\tcount\t{len(offsets_best)}",
    ]
    summary_txt.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print(f"WROTE {scoreboard_tsv}")
    print(f"WROTE {summary_txt}")
    print(f"WROTE {routing_tsv}")
    print(f"WROTE {routing_json}")
    print(f"WROTE {offsets_ready_json} ({len(offsets_ready)} subjects)")
    print(f"WROTE {offsets_best_json} ({len(offsets_best)} subjects)")
    print(
        f"subjects={len(subjects)} measured={measured} pending={pending} "
        f"pass={passed} fail={failed}"
    )


if __name__ == "__main__":
    main()
