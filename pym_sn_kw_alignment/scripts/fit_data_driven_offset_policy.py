#!/usr/bin/env python3
"""Fit a per-subject offset-gain policy from the latest corrected snapshot.

The script builds a clean training table by joining:
  - pre metrics (raw-KW vs v10 pre)
  - latest corrected post metrics (spot_*_hybrid_current.tsv)
  - source offset vectors (JSON used for that corrected run)

It then fits a lightweight response model to estimate subject-specific
"effective response" alpha and derives per-subject gains (not one global gain).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


def read_tsv(path: Path) -> List[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def parse_vec(row: dict, prefix: str = "delta") -> np.ndarray | None:
    try:
        x = float(row[f"{prefix}_x_mm"])
        y = float(row[f"{prefix}_y_mm"])
        z = float(row[f"{prefix}_z_mm"])
    except Exception:
        return None
    if any(math.isnan(v) for v in (x, y, z)):
        return None
    return np.array([x, y, z], dtype=float)


def clip(v: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, v))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument(
        "--source-offset-json",
        type=Path,
        help="Offset JSON used in the latest corrected run (default: hybrid_offsets_conservative.json)",
    )
    ap.add_argument(
        "--output-json",
        type=Path,
        help="Output data-driven offset JSON (default: <run>/multiopt_offsets/data_driven_policy_v1.json)",
    )
    ap.add_argument(
        "--output-training-tsv",
        type=Path,
        help="Output clean training table TSV (default: <run>/data_driven_training_latest.tsv)",
    )
    ap.add_argument(
        "--output-summary",
        type=Path,
        help="Output summary text file (default: <run>/data_driven_policy_summary.txt)",
    )
    ap.add_argument(
        "--gain-max",
        type=float,
        default=1.25,
        help="Maximum allowed subject gain (default: 1.25)",
    )
    args = ap.parse_args()

    run_dir: Path = args.run_dir
    source_offset_json = args.source_offset_json or (run_dir / "hybrid_offsets_conservative.json")
    out_json = args.output_json or (run_dir / "multiopt_offsets" / "data_driven_policy_v1.json")
    out_tsv = args.output_training_tsv or (run_dir / "data_driven_training_latest.tsv")
    out_summary = args.output_summary or (run_dir / "data_driven_policy_summary.txt")

    pre_highhat = run_dir / "highhat" / "highhat_rawkw_vs_v10_pre.tsv"
    pre_crash = run_dir / "crash" / "crash_rawkw_vs_v10_pre.tsv"
    post_highhat = run_dir / "spot_highhat_hybrid_current.tsv"
    post_crash = run_dir / "spot_crash_hybrid_current.tsv"

    for p in [pre_highhat, pre_crash, post_highhat, post_crash, source_offset_json]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required file: {p}")

    with source_offset_json.open(encoding="utf-8") as f:
        offset_data: Dict[str, dict] = json.load(f)

    pre_rows = {
        "highhat": {r["subject"]: r for r in read_tsv(pre_highhat)},
        "crash": {r["subject"]: r for r in read_tsv(pre_crash)},
    }
    post_rows = {
        "highhat": {r["subject"]: r for r in read_tsv(post_highhat)},
        "crash": {r["subject"]: r for r in read_tsv(post_crash)},
    }

    rows: List[dict] = []
    for group in ("highhat", "crash"):
        for sub, pre in pre_rows[group].items():
            p_vec = parse_vec(pre)
            if p_vec is None:
                continue
            p_mag = float(np.linalg.norm(p_vec))

            off = offset_data.get(sub, {})
            delta = off.get("delta_mm")
            if not isinstance(delta, list) or len(delta) != 3:
                continue
            try:
                s_vec = np.array([float(delta[0]), float(delta[1]), float(delta[2])], dtype=float)
            except Exception:
                continue
            if np.any(np.isnan(s_vec)):
                continue

            s_norm2 = float(np.dot(s_vec, s_vec))
            s_mag = float(np.linalg.norm(s_vec))
            dot_ps = float(np.dot(p_vec, s_vec))
            cos_ps = float(dot_ps / (p_mag * s_mag)) if p_mag > 0 and s_mag > 0 else float("nan")
            g_geom = float(dot_ps / s_norm2) if s_norm2 > 0 else float("nan")

            post = post_rows[group].get(sub)
            has_post = post is not None
            q_vec = parse_vec(post) if has_post else None
            q_mag = float("nan")
            imp_obs = float("nan")
            alpha_obs = float("nan")
            if q_vec is not None and s_norm2 > 0:
                q_mag = float(np.linalg.norm(q_vec))
                imp_obs = p_mag - q_mag
                alpha_obs = float(np.dot(p_vec - q_vec, s_vec) / s_norm2)

            rows.append(
                {
                    "subject": sub,
                    "group": group,
                    "pre_dx": float(p_vec[0]),
                    "pre_dy": float(p_vec[1]),
                    "pre_dz": float(p_vec[2]),
                    "pre_mag": p_mag,
                    "src_dx": float(s_vec[0]),
                    "src_dy": float(s_vec[1]),
                    "src_dz": float(s_vec[2]),
                    "src_mag": s_mag,
                    "dot_ps": dot_ps,
                    "cos_ps": cos_ps,
                    "g_geom": g_geom,
                    "has_post": int(has_post and q_vec is not None),
                    "post_dx": float(q_vec[0]) if q_vec is not None else float("nan"),
                    "post_dy": float(q_vec[1]) if q_vec is not None else float("nan"),
                    "post_dz": float(q_vec[2]) if q_vec is not None else float("nan"),
                    "post_mag": q_mag,
                    "improvement_obs": imp_obs,
                    "alpha_obs": alpha_obs,
                }
            )

    train = [
        r
        for r in rows
        if r["has_post"] == 1
        and math.isfinite(r["alpha_obs"])
        and r["src_mag"] > 1e-8
        and r["pre_mag"] > 1e-8
    ]
    if len(train) < 12:
        raise RuntimeError(f"Too few training rows ({len(train)}). Need at least 12.")

    # Robust clip for observed alpha before fitting.
    y_all = np.array([r["alpha_obs"] for r in train], dtype=float)
    q05, q95 = np.quantile(y_all, [0.05, 0.95])
    y = np.clip(y_all, q05, q95)

    x_cols: List[List[float]] = []
    for r in train:
        x_cols.append(
            [
                1.0,
                1.0 if r["group"] == "crash" else 0.0,
                math.log1p(r["pre_mag"]),
                math.log1p(r["src_mag"]),
                0.0 if not math.isfinite(r["cos_ps"]) else r["cos_ps"],
                clip(r["g_geom"], -3.0, 3.0) if math.isfinite(r["g_geom"]) else 0.0,
            ]
        )
    X = np.array(x_cols, dtype=float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)

    # Predicted alpha bounds from training distribution.
    q10, q90 = np.quantile(y, [0.10, 0.90])
    alpha_min = max(0.15, float(q10))
    alpha_max = max(alpha_min + 0.05, float(q90))

    def predict_alpha(r: dict) -> float:
        vec = np.array(
            [
                1.0,
                1.0 if r["group"] == "crash" else 0.0,
                math.log1p(r["pre_mag"]),
                math.log1p(r["src_mag"]),
                0.0 if not math.isfinite(r["cos_ps"]) else r["cos_ps"],
                clip(r["g_geom"], -3.0, 3.0) if math.isfinite(r["g_geom"]) else 0.0,
            ],
            dtype=float,
        )
        a = float(np.dot(vec, beta))
        return clip(a, alpha_min, alpha_max)

    # Fit a simple closeness-probability model to avoid binary max/zero behavior.
    y_cls = np.array([1.0 if r["improvement_obs"] > 0 else 0.0 for r in train], dtype=float)
    beta_cls, *_ = np.linalg.lstsq(X, y_cls, rcond=None)

    def predict_p_closer(r: dict) -> float:
        vec = np.array(
            [
                1.0,
                1.0 if r["group"] == "crash" else 0.0,
                math.log1p(r["pre_mag"]),
                math.log1p(r["src_mag"]),
                0.0 if not math.isfinite(r["cos_ps"]) else r["cos_ps"],
                clip(r["g_geom"], -3.0, 3.0) if math.isfinite(r["g_geom"]) else 0.0,
            ],
            dtype=float,
        )
        p = float(np.dot(vec, beta_cls))
        return clip(p, 0.0, 1.0)

    # Per-subject geometric gain target + learned alpha compensation.
    for r in rows:
        alpha_pred = predict_alpha(r)
        r["alpha_pred"] = alpha_pred
        r["p_closer_pred"] = predict_p_closer(r)
        if r["src_mag"] <= 1e-8 or r["dot_ps"] <= 0:
            g_star = 0.0
        else:
            g_star = r["dot_ps"] / (alpha_pred * (r["src_mag"] ** 2))
        r["gain_star"] = clip(g_star, 0.0, args.gain_max)

    # Calibrate one global safety factor using observed training response.
    candidates = [round(v, 2) for v in np.arange(0.35, 1.01, 0.05)]
    best = None
    for eta in candidates:
        closer = farther = neutral = 0
        imp_vals: List[float] = []
        for r in train:
            g = clip(eta * r["gain_star"] * r["p_closer_pred"], 0.0, args.gain_max)
            p = np.array([r["pre_dx"], r["pre_dy"], r["pre_dz"]], dtype=float)
            s = np.array([r["src_dx"], r["src_dy"], r["src_dz"]], dtype=float)
            q_est = p - r["alpha_obs"] * g * s
            imp = r["pre_mag"] - float(np.linalg.norm(q_est))
            imp_vals.append(imp)
            if imp > 1e-6:
                closer += 1
            elif imp < -1e-6:
                farther += 1
            else:
                neutral += 1
        mean_imp = float(np.mean(imp_vals))
        key = (closer, -farther, mean_imp)
        if best is None or key > best["key"]:
            best = {
                "eta": eta,
                "closer": closer,
                "farther": farther,
                "neutral": neutral,
                "mean_imp": mean_imp,
                "key": key,
            }
    assert best is not None
    eta = best["eta"]

    for r in rows:
        r["gain_final"] = clip(eta * r["gain_star"] * r["p_closer_pred"], 0.0, args.gain_max)

    # Write output JSON (all subjects with a source offset).
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_map: Dict[str, dict] = {}
    for r in rows:
        sub = r["subject"]
        base = dict(offset_data.get(sub, {}))
        g = r["gain_final"]
        src = np.array([r["src_dx"], r["src_dy"], r["src_dz"]], dtype=float)
        scaled = (src * g).tolist()
        base["delta_mm"] = [float(scaled[0]), float(scaled[1]), float(scaled[2])]
        base["policy"] = "data_driven_v1"
        base["data_driven_gain"] = float(g)
        base["data_driven_eta"] = float(eta)
        base["data_driven_alpha_pred"] = float(r["alpha_pred"])
        base["data_driven_train_has_post"] = int(r["has_post"])
        out_map[sub] = base
    with out_json.open("w", encoding="utf-8") as f:
        json.dump(out_map, f, indent=2, sort_keys=True)

    # Write clean training table.
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "subject",
        "group",
        "has_post",
        "pre_dx",
        "pre_dy",
        "pre_dz",
        "pre_mag",
        "src_dx",
        "src_dy",
        "src_dz",
        "src_mag",
        "dot_ps",
        "cos_ps",
        "g_geom",
        "post_dx",
        "post_dy",
        "post_dz",
        "post_mag",
        "improvement_obs",
        "alpha_obs",
        "alpha_pred",
        "p_closer_pred",
        "gain_star",
        "gain_final",
    ]
    with out_tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["group"], int(x["subject"]))):
            row = {}
            for k in fields:
                v = r.get(k, "")
                if isinstance(v, float):
                    row[k] = f"{v:.6f}" if math.isfinite(v) else ""
                else:
                    row[k] = v
            w.writerow(row)

    # Summary.
    train_n = len(train)
    full_n = len(rows)
    gain_vals = [r["gain_final"] for r in rows]
    nonzero = [g for g in gain_vals if g > 1e-6]
    with out_summary.open("w", encoding="utf-8") as f:
        f.write(f"run_dir\t{run_dir}\n")
        f.write(f"source_offset_json\t{source_offset_json}\n")
        f.write(f"rows_total\t{full_n}\n")
        f.write(f"rows_training_has_post\t{train_n}\n")
        f.write(f"alpha_clip_q05_q95\t{q05:.6f}\t{q95:.6f}\n")
        f.write(f"alpha_pred_bounds\t{alpha_min:.6f}\t{alpha_max:.6f}\n")
        f.write(
            "beta\t"
            + "\t".join(
                f"{v:.6f}" for v in beta
            )
            + "\n"
        )
        f.write(
            "beta_cls\t"
            + "\t".join(
                f"{v:.6f}" for v in beta_cls
            )
            + "\n"
        )
        f.write(
            f"eta_selected\t{eta:.2f}\ttrain_closer={best['closer']}\ttrain_farther={best['farther']}"
            f"\ttrain_neutral={best['neutral']}\ttrain_mean_imp={best['mean_imp']:.6f}\n"
        )
        f.write(
            f"gain_distribution\tmin={min(gain_vals):.6f}\tmedian={float(np.median(gain_vals)):.6f}"
            f"\tmax={max(gain_vals):.6f}\tnonzero={len(nonzero)}\n"
        )
        f.write(f"output_json\t{out_json}\n")
        f.write(f"output_training_tsv\t{out_tsv}\n")

    print(f"WROTE {out_json}")
    print(f"WROTE {out_tsv}")
    print(f"WROTE {out_summary}")
    print(
        f"TRAIN rows={train_n} / {full_n}; eta={eta:.2f}; "
        f"pred_train closer={best['closer']} farther={best['farther']} neutral={best['neutral']}"
    )


if __name__ == "__main__":
    main()
