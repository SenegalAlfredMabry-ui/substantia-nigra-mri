#!/usr/bin/env python3
"""Train a no-KW offset policy (kNN + anchor-aware fallback).

Inputs:
  - training dataset TSV (with target_delta_* columns)
  - deployment dataset TSV (same feature schema, target columns optional/blank)

Outputs:
  - model summary
  - training LOOCV diagnostics
  - deployment predictions
  - selected offset JSON (per-subject delta_mm + confidence + QC)

New in v2:
  - Optional multi-approach mode that evaluates anchor subsets and selects the
    best offset approach per deployment subject.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--train-dataset", required=True, type=Path)
    ap.add_argument("--deploy-dataset", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--k", type=int, default=9, help="k neighbors for kNN regression.")
    ap.add_argument(
        "--zmax-gate",
        type=float,
        default=3.0,
        help="Feature z-score gate for selecting kNN over global fallback.",
    )
    ap.add_argument(
        "--distance-quantile",
        type=float,
        default=0.90,
        help="LOOCV nearest-distance quantile used as kNN selection threshold.",
    )
    ap.add_argument(
        "--mag-limit-scale",
        type=float,
        default=1.5,
        help="Scale factor applied to p95 training target magnitude for deployment gate.",
    )
    ap.add_argument(
        "--offset-source-tag",
        default="no_kw_knn_anchor_v1",
        help="Tag saved in output JSON entries.",
    )
    ap.add_argument(
        "--min-train-rows",
        type=int,
        default=20,
        help="Minimum rows required for supervised fitting (default: 20).",
    )
    ap.add_argument(
        "--approach-mode",
        choices=["combined", "auto"],
        default="combined",
        help=(
            "combined: use all anchor features together (legacy behavior). "
            "auto: evaluate anchor subsets and select best approach per subject."
        ),
    )
    ap.add_argument(
        "--approach-select-metric",
        choices=["loocv_mean_err_knn_mm", "loocv_knn_better_rate"],
        default="loocv_mean_err_knn_mm",
        help="Metric used to rank global best approach in auto mode.",
    )
    ap.add_argument(
        "--gate-feature-mode",
        choices=["all", "anchor_relative"],
        default="all",
        help=(
            "Feature subset used for in-distribution gates (distance/zmax). "
            "all: use all modeling features; anchor_relative: use anchor-relative "
            "features (vec/dist/anchor_present + candidate_lr_imbalance)."
        ),
    )
    ap.add_argument(
        "--force-knn-if-input-present",
        type=int,
        choices=[0, 1],
        default=0,
        help=(
            "If 1, force kNN selection for deploy rows with candidate+anchor inputs "
            "even when gates fail; still records gate failures in reason."
        ),
    )
    ap.add_argument(
        "--force-knn-clamp-to-mag-limit",
        type=int,
        choices=[0, 1],
        default=1,
        help="If forcing kNN, clamp predicted vector magnitude to mag-limit when exceeded.",
    )
    ap.add_argument(
        "--allow-partial-anchor-input",
        type=int,
        choices=[0, 1],
        default=1,
        help=(
            "If 1, allow an anchor approach to run when at least one required anchor is "
            "present (missing required-anchor features are imputed)."
        ),
    )
    ap.add_argument(
        "--include-candidate-only-in-auto",
        type=int,
        choices=[0, 1],
        default=1,
        help=(
            "If 1 and approach-mode=auto, include candidate_only (no anchors) as a "
            "backstop so deploy rows with missing anchors can still receive offsets."
        ),
    )
    ap.add_argument(
        "--deploy-selection-policy",
        choices=["default", "local_expected_gain"],
        default="default",
        help=(
            "default: existing gate+confidence routing. "
            "local_expected_gain: require locally estimated expected gain before applying kNN."
        ),
    )
    ap.add_argument(
        "--min-expected-gain-mm",
        type=float,
        default=0.25,
        help=(
            "Minimum locally estimated expected gain (mm) required when "
            "--deploy-selection-policy=local_expected_gain."
        ),
    )
    ap.add_argument(
        "--max-selected-mag-mm",
        type=float,
        default=float("nan"),
        help=(
            "Optional guardrail: reject kNN selection when predicted shift magnitude "
            "exceeds this value (disabled when not finite)."
        ),
    )
    ap.add_argument(
        "--fallback-mode",
        choices=["global_offset", "no_correction"],
        default="global_offset",
        help=(
            "Fallback when kNN is rejected. "
            "global_offset: use global mean vector; no_correction: use [0,0,0]."
        ),
    )
    ap.add_argument(
        "--require-positive-direction-agreement",
        type=int,
        choices=[0, 1],
        default=0,
        help=(
            "If 1, require local expected directional agreement to exceed "
            "--min-expected-direction-agreement before applying kNN."
        ),
    )
    ap.add_argument(
        "--min-expected-direction-agreement",
        type=float,
        default=0.0,
        help=(
            "Minimum local expected direction agreement (cosine-like score, -1..1) "
            "required when --require-positive-direction-agreement=1."
        ),
    )
    ap.add_argument(
        "--candidate-only-min-expected-gain-mm",
        type=float,
        default=float("nan"),
        help=(
            "Optional stricter expected-gain gate applied only to candidate_only "
            "approach (disabled when not finite)."
        ),
    )
    ap.add_argument(
        "--candidate-only-min-confidence",
        type=float,
        default=float("nan"),
        help=(
            "Optional stricter confidence gate applied only to candidate_only "
            "approach (disabled when not finite)."
        ),
    )
    ap.add_argument(
        "--uncertainty-confidence-threshold",
        type=float,
        default=float("nan"),
        help=(
            "If finite, force no_correction fallback when selected confidence is "
            "below this threshold."
        ),
    )
    ap.add_argument(
        "--soft-rescue-min-confidence",
        type=float,
        default=float("nan"),
        help=(
            "Optional low-risk rescue: if kNN was rejected only by conservative "
            "selection thresholds but remains in-distribution, allow it when "
            "confidence is at least this value (disabled when not finite)."
        ),
    )
    ap.add_argument(
        "--soft-rescue-min-expected-gain-mm",
        type=float,
        default=float("nan"),
        help=(
            "Optional low-risk rescue gain floor paired with "
            "--soft-rescue-min-confidence (disabled when not finite)."
        ),
    )
    ap.add_argument(
        "--soft-rescue-min-direction-agreement",
        type=float,
        default=float("nan"),
        help=(
            "Optional low-risk rescue direction-agreement floor paired with "
            "--soft-rescue-min-confidence (disabled when not finite)."
        ),
    )
    ap.add_argument(
        "--soft-rescue-scale",
        type=float,
        default=1.0,
        help=(
            "Optional scale applied to rescued kNN vectors. Use <1.0 to damp "
            "soft-rescue moves."
        ),
    )
    return ap.parse_args()


def read_tsv(path: Path) -> List[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def to_float(row: dict, key: str) -> float:
    try:
        val = row.get(key, "")
        if val is None or val == "":
            return float("nan")
        return float(val)
    except Exception:
        return float("nan")


def finite(v: float) -> bool:
    return math.isfinite(v)


def is_knn_method(method: str) -> bool:
    return method.startswith("knn_")


def is_fallback_method(method: str) -> bool:
    return method in {"global_fallback", "no_correction_fallback"}


def is_anchor_relative_feature(col: str) -> bool:
    return (
        col.startswith("vec_")
        or col.startswith("dist_")
        or col.startswith("anchor_present__")
        or col == "candidate_lr_imbalance"
    )


def nearest_neighbor_distance(x: np.ndarray, X_ref: np.ndarray) -> float:
    if X_ref.shape[0] == 0:
        return math.nan
    d = np.linalg.norm(X_ref - x[None, :], axis=1)
    if d.size == 0:
        return math.nan
    return float(np.min(d))


def collect_feature_columns(train_rows: Sequence[dict], deploy_rows: Sequence[dict]) -> List[str]:
    if not train_rows:
        return []
    blacklist_prefixes = (
        "target_",
        "manual_",
        "anchor_error_",
        "summary__",
    )
    blacklist_exact = {
        "subject",
        "candidate_mask",
        "manual_mask",
        "candidate_present",
        "manual_present",
    }
    cols: List[str] = []
    for key in train_rows[0].keys():
        if key in blacklist_exact:
            continue
        if key.startswith(blacklist_prefixes):
            continue
        keep = (
            key.startswith("candidate_")
            or key.startswith("vec_")
            or key.startswith("dist_")
            or key.startswith("anchor_present__")
        )
        if keep:
            cols.append(key)

    deploy_keys = set(deploy_rows[0].keys()) if deploy_rows else set()
    return [c for c in cols if c in deploy_keys]


def discover_anchor_names(feature_cols: Sequence[str]) -> List[str]:
    names = {
        c.split("anchor_present__", 1)[1]
        for c in feature_cols
        if c.startswith("anchor_present__")
    }
    return sorted(names)


def feature_columns_for_anchors(feature_cols: Sequence[str], anchors: Sequence[str]) -> List[str]:
    if not anchors:
        return [c for c in feature_cols if "__" not in c]
    suffixes = tuple(f"__{a}" for a in anchors)
    out: List[str] = []
    for c in feature_cols:
        if "__" not in c:
            out.append(c)
            continue
        if c.endswith(suffixes):
            out.append(c)
    return out


def approach_id_for(anchors: Sequence[str], all_anchors: Sequence[str]) -> str:
    if list(anchors) == list(all_anchors):
        return "all_anchors"
    if len(anchors) == 1:
        return f"only__{anchors[0]}"
    return "combo__" + "__".join(anchors)


def build_matrix(rows: Sequence[dict], feature_cols: Sequence[str]) -> np.ndarray:
    X = np.full((len(rows), len(feature_cols)), np.nan, dtype=float)
    for i, row in enumerate(rows):
        for j, col in enumerate(feature_cols):
            X[i, j] = to_float(row, col)
    return X


def impute_and_standardize(
    X_train: np.ndarray, X_other: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    means = np.nanmean(X_train, axis=0)
    means = np.where(np.isfinite(means), means, 0.0)

    X_train_i = np.where(np.isfinite(X_train), X_train, means)
    X_other_i = np.where(np.isfinite(X_other), X_other, means)

    stds = np.nanstd(X_train_i, axis=0)
    stds = np.where(stds > 1e-8, stds, 1.0)

    X_train_z = (X_train_i - means) / stds
    X_other_z = (X_other_i - means) / stds
    return X_train_z, X_other_z, means, stds


def knn_predict(
    x: np.ndarray,
    X_ref: np.ndarray,
    Y_ref: np.ndarray,
    k: int,
) -> Tuple[np.ndarray, float, float, np.ndarray, np.ndarray]:
    if X_ref.shape[0] == 0:
        return (
            np.array([math.nan, math.nan, math.nan]),
            math.nan,
            math.nan,
            np.array([], dtype=int),
            np.array([], dtype=float),
        )
    d = np.linalg.norm(X_ref - x[None, :], axis=1)
    order = np.argsort(d)
    kk = max(1, min(k, X_ref.shape[0]))
    idx = order[:kk]
    d_sel = d[idx]
    y_sel = Y_ref[idx]
    w = 1.0 / (d_sel + 1e-6)
    w = w / w.sum()
    pred = (w[:, None] * y_sel).sum(axis=0)
    mags = np.linalg.norm(y_sel, axis=1)
    spread = float(np.sqrt((w * (mags - np.sum(w * mags)) ** 2).sum()))
    return pred, float(d_sel[0]), spread, idx, w


def zmax(x: np.ndarray) -> float:
    if x.size == 0:
        return math.nan
    return float(np.max(np.abs(x)))


def directional_agreement(pred: np.ndarray, true: np.ndarray) -> float:
    if pred.shape != (3,) or true.shape != (3,):
        return math.nan
    if not (np.all(np.isfinite(pred)) and np.all(np.isfinite(true))):
        return math.nan
    pred_n = float(np.linalg.norm(pred))
    true_n = float(np.linalg.norm(true))
    if pred_n <= 1e-8 or true_n <= 1e-8:
        return math.nan
    return float(np.dot(pred, true) / (pred_n * true_n))


def write_tsv(path: Path, rows: Sequence[dict], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def f6(v: float) -> str:
    return f"{v:.6f}" if finite(v) else ""


def run_approach(
    args: argparse.Namespace,
    approach_name: str,
    approach_anchors: Sequence[str],
    feature_cols: Sequence[str],
    train_rows: Sequence[dict],
    deploy_rows: Sequence[dict],
    Y_train: np.ndarray,
) -> dict:
    local_cols = feature_columns_for_anchors(feature_cols, approach_anchors)
    if not local_cols:
        raise SystemExit(f"Approach {approach_name} has no feature columns.")

    gate_cols = list(local_cols)
    if args.gate_feature_mode == "anchor_relative":
        gate_cols = [c for c in local_cols if is_anchor_relative_feature(c)]
    if not gate_cols:
        gate_cols = list(local_cols)

    X_train_raw = build_matrix(train_rows, local_cols)
    X_deploy_raw = build_matrix(deploy_rows, local_cols)
    X_train, X_deploy, feat_means, feat_stds = impute_and_standardize(X_train_raw, X_deploy_raw)

    X_train_gate_raw = build_matrix(train_rows, gate_cols)
    X_deploy_gate_raw = build_matrix(deploy_rows, gate_cols)
    X_train_gate, X_deploy_gate, gate_means, gate_stds = impute_and_standardize(
        X_train_gate_raw, X_deploy_gate_raw
    )

    y_mag = np.linalg.norm(Y_train, axis=1)
    if args.fallback_mode == "no_correction":
        global_vec = np.array([0.0, 0.0, 0.0], dtype=float)
    else:
        global_vec = Y_train.mean(axis=0)
    p95_mag = float(np.percentile(y_mag, 95))
    mag_limit = args.mag_limit_scale * p95_mag

    loocv_rows: List[dict] = []
    nearest_dists: List[float] = []
    err_knn_all: List[float] = []
    err_global_all: List[float] = []
    dir_knn_all: List[float] = []
    for i, row in enumerate(train_rows):
        keep = np.ones(len(train_rows), dtype=bool)
        keep[i] = False
        pred_knn, _, spread, _idx, _w = knn_predict(X_train[i], X_train[keep], Y_train[keep], args.k)
        d0 = nearest_neighbor_distance(X_train_gate[i], X_train_gate[keep])
        pred_global = Y_train[keep].mean(axis=0)
        true = Y_train[i]
        err_knn = float(np.linalg.norm(pred_knn - true))
        err_global = float(np.linalg.norm(pred_global - true))
        dir_knn = directional_agreement(pred_knn, true)
        nearest_dists.append(d0)
        err_knn_all.append(err_knn)
        err_global_all.append(err_global)
        dir_knn_all.append(dir_knn)
        loocv_rows.append(
            {
                "subject": row["subject"],
                "approach": approach_name,
                "true_dx_mm": f6(true[0]),
                "true_dy_mm": f6(true[1]),
                "true_dz_mm": f6(true[2]),
                "knn_dx_mm": f6(pred_knn[0]),
                "knn_dy_mm": f6(pred_knn[1]),
                "knn_dz_mm": f6(pred_knn[2]),
                "global_dx_mm": f6(pred_global[0]),
                "global_dy_mm": f6(pred_global[1]),
                "global_dz_mm": f6(pred_global[2]),
                "err_knn_mm": f6(err_knn),
                "err_global_mm": f6(err_global),
                "direction_agreement_knn": f6(dir_knn),
                "nearest_dist": f6(d0),
                "neighbor_spread_mm": f6(spread),
                "knn_better": "1" if err_knn < err_global else "0",
            }
        )

    loocv_err_arr = np.array(err_knn_all, dtype=float)
    loocv_dir_arr = np.array(dir_knn_all, dtype=float)
    dist_thr = float(np.quantile(np.array(nearest_dists, dtype=float), args.distance_quantile))

    anchor_cols_required = [f"anchor_present__{a}" for a in approach_anchors]

    pred_rows: List[dict] = []
    deploy_eval: List[dict] = []
    selected_offsets: Dict[str, dict] = {}
    train_ref_mag_p75 = float(np.percentile(y_mag, 75))

    for i, row in enumerate(deploy_rows):
        sid = row["subject"]
        cand_ok = row.get("candidate_present") == "1"
        present_anchors = [a for a in approach_anchors if row.get(f"anchor_present__{a}", "0") == "1"]
        missing_anchors = [a for a in approach_anchors if row.get(f"anchor_present__{a}", "0") != "1"]
        anchor_ok = len(missing_anchors) == 0
        if (
            not anchor_ok
            and args.allow_partial_anchor_input == 1
            and len(approach_anchors) > 0
            and len(present_anchors) > 0
        ):
            anchor_ok = True

        method = "missing_input"
        reason: List[str] = []
        pred_knn = np.array([math.nan, math.nan, math.nan], dtype=float)
        pred_global = global_vec.copy()
        pred_sel = np.array([math.nan, math.nan, math.nan], dtype=float)
        d0 = spread = z = conf = math.nan
        expected_baseline_mag_mm = math.nan
        expected_post_err_mm = math.nan
        expected_gain_mm = math.nan
        expected_direction_agreement = math.nan
        knn_mag = math.nan
        qc_pass = 0
        gate_dist_ok = False
        gate_zmax_ok = False
        gate_mag_ok = False
        gate_dir_ok = (
            args.require_positive_direction_agreement == 0
        )  # default true when feature disabled
        prefer_knn = False
        direction_rejected = False

        if not cand_ok:
            reason.append("missing_candidate_mask")
        if not anchor_ok:
            reason.append("missing_anchor_mask")
        elif len(missing_anchors) > 0:
            reason.append("partial_anchor_input")

        if cand_ok and anchor_ok:
            pred_knn, _, spread, nn_idx, nn_w = knn_predict(X_deploy[i], X_train, Y_train, args.k)
            d0 = nearest_neighbor_distance(X_deploy_gate[i], X_train_gate)
            z = zmax(X_deploy_gate[i])
            knn_mag = float(np.linalg.norm(pred_knn))
            if nn_idx.size > 0 and nn_w.size > 0:
                expected_baseline_mag_mm = float(np.sum(nn_w * y_mag[nn_idx]))
                expected_post_err_mm = float(np.sum(nn_w * loocv_err_arr[nn_idx]))
                expected_gain_mm = expected_baseline_mag_mm - expected_post_err_mm
                dir_vals = loocv_dir_arr[nn_idx]
                finite_dir = np.isfinite(dir_vals)
                if np.any(finite_dir):
                    w_eff = nn_w[finite_dir]
                    w_sum = float(np.sum(w_eff))
                    if w_sum > 1e-8:
                        w_eff = w_eff / w_sum
                        expected_direction_agreement = float(np.sum(w_eff * dir_vals[finite_dir]))

            gate_dist_ok = finite(d0) and d0 <= dist_thr
            gate_zmax_ok = finite(z) and z <= args.zmax_gate
            gate_mag_ok = finite(knn_mag) and knn_mag <= mag_limit
            if finite(d0):
                c_dist = max(0.0, 1.0 - min(1.0, d0 / (dist_thr + 1e-6)))
            else:
                c_dist = 0.0
            if finite(spread):
                denom = max(train_ref_mag_p75, 1e-6)
                c_spread = max(0.0, 1.0 - min(1.0, spread / denom))
            else:
                c_spread = 0.0
            conf = 0.5 * (c_dist + c_spread)
            gate_dir_ok = (
                args.require_positive_direction_agreement == 0
                or (
                    finite(expected_direction_agreement)
                    and expected_direction_agreement > args.min_expected_direction_agreement
                )
            )
            prefer_knn = gate_dist_ok and gate_zmax_ok and gate_mag_ok
            if (
                prefer_knn
                and args.deploy_selection_policy == "local_expected_gain"
                and not (
                    finite(expected_gain_mm)
                    and expected_gain_mm >= args.min_expected_gain_mm
                )
            ):
                prefer_knn = False
                reason.append("expected_gain_low")
            if prefer_knn and not gate_dir_ok:
                prefer_knn = False
                direction_rejected = True
                reason.append("direction_disagreement")
            if prefer_knn and approach_name == "candidate_only":
                if finite(args.candidate_only_min_expected_gain_mm) and not (
                    finite(expected_gain_mm)
                    and expected_gain_mm >= args.candidate_only_min_expected_gain_mm
                ):
                    prefer_knn = False
                    reason.append("candidate_only_expected_gain_low")
                if (
                    prefer_knn
                    and finite(args.candidate_only_min_confidence)
                    and not (finite(conf) and conf >= args.candidate_only_min_confidence)
                ):
                    prefer_knn = False
                    reason.append("candidate_only_low_confidence")
            if (
                prefer_knn
                and finite(args.max_selected_mag_mm)
                and finite(knn_mag)
                and knn_mag > args.max_selected_mag_mm
            ):
                prefer_knn = False
                reason.append("selected_magnitude_guardrail")

            soft_rescue_ok = False
            if not prefer_knn and finite(knn_mag):
                soft_rescue_ok = (
                    finite(args.soft_rescue_min_confidence)
                    and finite(args.soft_rescue_min_expected_gain_mm)
                    and finite(args.soft_rescue_min_direction_agreement)
                    and gate_dist_ok
                    and gate_zmax_ok
                    and gate_mag_ok
                    and gate_dir_ok
                    and finite(conf)
                    and conf >= args.soft_rescue_min_confidence
                    and finite(expected_gain_mm)
                    and expected_gain_mm >= args.soft_rescue_min_expected_gain_mm
                    and finite(expected_direction_agreement)
                    and expected_direction_agreement
                    >= args.soft_rescue_min_direction_agreement
                )

            if prefer_knn:
                method = "knn_anchor"
                pred_sel = pred_knn
                qc_pass = 1
            elif soft_rescue_ok:
                method = "knn_soft_rescue"
                pred_sel = pred_knn.copy()
                qc_pass = 1
                reason.append("soft_rescue_override")
                if finite(args.soft_rescue_scale) and args.soft_rescue_scale != 1.0:
                    pred_sel = pred_sel * args.soft_rescue_scale
                    reason.append("soft_rescue_scaled")
            else:
                if not (finite(d0) and d0 <= dist_thr):
                    reason.append("neighbor_distance_high")
                if not (finite(z) and z <= args.zmax_gate):
                    reason.append("feature_out_of_distribution")
                if not (finite(knn_mag) and knn_mag <= mag_limit):
                    reason.append("knn_magnitude_high")
                if not gate_dir_ok:
                    reason.append("direction_gate_fail")

                if args.force_knn_if_input_present == 1 and finite(knn_mag):
                    method = "knn_forced"
                    pred_sel = pred_knn
                    qc_pass = 1
                    reason.append("forced_knn_override")
                    if (
                        args.force_knn_clamp_to_mag_limit == 1
                        and finite(mag_limit)
                        and mag_limit > 0.0
                        and knn_mag > mag_limit
                    ):
                        pred_sel = pred_knn * (mag_limit / knn_mag)
                        method = "knn_forced_clamped"
                        reason.append("knn_scaled_to_mag_limit")
                else:
                    if direction_rejected:
                        method = "no_correction_fallback"
                        pred_sel = np.array([0.0, 0.0, 0.0], dtype=float)
                        reason.append("direction_gate_no_correction_fallback")
                    else:
                        method = (
                            "no_correction_fallback"
                            if args.fallback_mode == "no_correction"
                            else "global_fallback"
                        )
                        pred_sel = pred_global
                    qc_pass = 1
                    if args.force_knn_if_input_present == 1 and not finite(knn_mag):
                        reason.append("forced_knn_unavailable_nonfinite_prediction")

            if (
                finite(args.uncertainty_confidence_threshold)
                and method != "missing_input"
                and (not finite(conf) or conf < args.uncertainty_confidence_threshold)
            ):
                method = "no_correction_fallback"
                pred_sel = np.array([0.0, 0.0, 0.0], dtype=float)
                qc_pass = 1
                reason.append("uncertainty_no_correction_fallback")

            selected_offsets[sid] = {
                "delta_mm": [float(pred_sel[0]), float(pred_sel[1]), float(pred_sel[2])],
                "source": args.offset_source_tag,
                "approach": approach_name,
                "method": method,
                "confidence": float(conf),
                "qc": {
                    "passed": int(qc_pass),
                    "nearest_feature_distance": float(d0) if finite(d0) else None,
                    "feature_zmax": float(z) if finite(z) else None,
                    "neighbor_target_spread_mm": float(spread) if finite(spread) else None,
                    "expected_baseline_mag_mm": float(expected_baseline_mag_mm)
                    if finite(expected_baseline_mag_mm)
                    else None,
                    "expected_post_err_mm": float(expected_post_err_mm)
                    if finite(expected_post_err_mm)
                    else None,
                    "expected_gain_mm": float(expected_gain_mm)
                    if finite(expected_gain_mm)
                    else None,
                    "expected_direction_agreement": float(expected_direction_agreement)
                    if finite(expected_direction_agreement)
                    else None,
                    "distance_threshold": float(dist_thr),
                    "zmax_gate": float(args.zmax_gate),
                    "predicted_mag_limit_mm": float(mag_limit),
                    "min_expected_direction_agreement": float(args.min_expected_direction_agreement),
                    "reason": ";".join(reason),
                },
            }

        deploy_eval.append(
            {
                "subject": sid,
                "approach": approach_name,
                "anchors_required": list(approach_anchors),
                "anchors_present": list(present_anchors),
                "anchors_missing": list(missing_anchors),
                "method": method,
                "qc_pass": int(qc_pass),
                "reason": ";".join(reason),
                "pred_knn": pred_knn,
                "pred_global": pred_global,
                "pred_sel": pred_sel,
                "d0": d0,
                "d0_threshold": float(dist_thr),
                "z": z,
                "zmax_gate": float(args.zmax_gate),
                "spread": spread,
                "expected_baseline_mag_mm": expected_baseline_mag_mm,
                "expected_post_err_mm": expected_post_err_mm,
                "expected_gain_mm": expected_gain_mm,
                "expected_direction_agreement": expected_direction_agreement,
                "conf": conf,
                "knn_mag": knn_mag,
                "knn_mag_limit": float(mag_limit),
                "cand_ok": cand_ok,
                "anchor_ok": anchor_ok,
                "gate_dist_ok": gate_dist_ok,
                "gate_zmax_ok": gate_zmax_ok,
                "gate_mag_ok": gate_mag_ok,
                "gate_dir_ok": gate_dir_ok,
                "prefer_knn": prefer_knn,
            }
        )

        pred_rows.append(
            {
                "subject": sid,
                "approach": approach_name,
                "method": method,
                "qc_pass": str(qc_pass),
                "reason": ";".join(reason),
                "knn_dx_mm": f6(pred_knn[0]),
                "knn_dy_mm": f6(pred_knn[1]),
                "knn_dz_mm": f6(pred_knn[2]),
                "global_dx_mm": f6(pred_global[0]),
                "global_dy_mm": f6(pred_global[1]),
                "global_dz_mm": f6(pred_global[2]),
                "selected_dx_mm": f6(pred_sel[0]),
                "selected_dy_mm": f6(pred_sel[1]),
                "selected_dz_mm": f6(pred_sel[2]),
                "selected_mag_mm": f6(
                    float(np.linalg.norm(pred_sel)) if np.all(np.isfinite(pred_sel)) else math.nan
                ),
                "nearest_feature_distance": f6(d0),
                "feature_zmax": f6(z),
                "neighbor_target_spread_mm": f6(spread),
                "expected_baseline_mag_mm": f6(expected_baseline_mag_mm),
                "expected_post_err_mm": f6(expected_post_err_mm),
                "expected_gain_mm": f6(expected_gain_mm),
                "expected_direction_agreement": f6(expected_direction_agreement),
                "confidence": f6(conf),
            }
        )

    mean_err_knn = float(np.mean(err_knn_all))
    mean_err_global = float(np.mean(err_global_all))
    knn_better_rate = float(np.mean(np.array(err_knn_all) < np.array(err_global_all)))
    finite_dir = np.array([v for v in dir_knn_all if finite(v)], dtype=float)
    mean_dir_knn = float(np.mean(finite_dir)) if finite_dir.size else math.nan
    deploy_knn = sum(1 for r in pred_rows if is_knn_method(r["method"]))
    deploy_knn_anchor = sum(1 for r in pred_rows if r["method"] == "knn_anchor")
    deploy_knn_forced = sum(1 for r in pred_rows if r["method"].startswith("knn_forced"))
    deploy_knn_soft_rescue = sum(
        1 for r in pred_rows if r["method"] == "knn_soft_rescue"
    )
    deploy_global = sum(1 for r in pred_rows if is_fallback_method(r["method"]))
    deploy_missing = sum(1 for r in pred_rows if r["method"] == "missing_input")

    return {
        "approach": approach_name,
        "anchors": list(approach_anchors),
        "feature_cols": list(local_cols),
        "feature_means": {c: float(m) for c, m in zip(local_cols, feat_means)},
        "feature_stds": {c: float(s) for c, s in zip(local_cols, feat_stds)},
        "gate_feature_mode": args.gate_feature_mode,
        "gate_feature_cols": list(gate_cols),
        "gate_feature_means": {c: float(m) for c, m in zip(gate_cols, gate_means)},
        "gate_feature_stds": {c: float(s) for c, s in zip(gate_cols, gate_stds)},
        "global_vec_mm": [float(global_vec[0]), float(global_vec[1]), float(global_vec[2])],
        "distance_threshold": float(dist_thr),
        "zmax_gate": float(args.zmax_gate),
        "predicted_mag_limit_mm": float(mag_limit),
        "force_knn_if_input_present": int(args.force_knn_if_input_present),
        "force_knn_clamp_to_mag_limit": int(args.force_knn_clamp_to_mag_limit),
        "loocv_mean_err_knn_mm": mean_err_knn,
        "loocv_mean_err_global_mm": mean_err_global,
        "loocv_knn_better_rate": knn_better_rate,
        "loocv_mean_direction_agreement_knn": mean_dir_knn,
        "deploy_method_knn_total": deploy_knn,
        "deploy_method_knn_anchor": deploy_knn_anchor,
        "deploy_method_knn_forced": deploy_knn_forced,
        "deploy_method_knn_soft_rescue": deploy_knn_soft_rescue,
        "deploy_method_global_fallback": deploy_global,
        "deploy_method_missing_input": deploy_missing,
        "loocv_rows": loocv_rows,
        "pred_rows": pred_rows,
        "deploy_eval": deploy_eval,
        "selected_offsets": selected_offsets,
    }


def choose_best_eval_for_subject(
    sid: str,
    evals: Sequence[dict],
    global_best_approach: str,
) -> Tuple[dict, str]:
    # Prefer kNN-passing approaches with highest confidence.
    knn = [e for e in evals if is_knn_method(e["method"]) and e["qc_pass"] == 1]
    if knn:
        knn.sort(
            key=lambda e: (
                e["expected_post_err_mm"]
                if finite(e["expected_post_err_mm"])
                else float("inf"),
                e["knn_mag"] if finite(e["knn_mag"]) else float("inf"),
                -(
                    e["expected_gain_mm"]
                    if finite(e["expected_gain_mm"])
                    else -float("inf")
                ),
                -(e["conf"] if finite(e["conf"]) else -1.0),
                e["d0"] if finite(e["d0"]) else float("inf"),
            )
        )
        return knn[0], "best_knn_by_expected_post_error_then_mag_then_gain_then_confidence"

    # Otherwise choose the most in-distribution fallback.
    fallback = [e for e in evals if is_fallback_method(e["method"]) and e["qc_pass"] == 1]
    if fallback:
        fallback.sort(
            key=lambda e: (
                e["d0"] if finite(e["d0"]) else float("inf"),
                -(e["conf"] if finite(e["conf"]) else -1.0),
            )
        )
        return fallback[0], "best_fallback_by_lowest_distance_then_confidence"

    # Last resort: pick global-best approach row (or first).
    for e in evals:
        if e["approach"] == global_best_approach:
            return e, "last_resort_global_best_approach"
    return evals[0], "last_resort_first_available"


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    train_rows_all = read_tsv(args.train_dataset)
    deploy_rows = read_tsv(args.deploy_dataset)
    if not train_rows_all:
        raise SystemExit("Training dataset is empty.")
    if not deploy_rows:
        raise SystemExit("Deployment dataset is empty.")

    train_rows: List[dict] = []
    for row in train_rows_all:
        dx = to_float(row, "target_delta_x_mm")
        dy = to_float(row, "target_delta_y_mm")
        dz = to_float(row, "target_delta_z_mm")
        if not (finite(dx) and finite(dy) and finite(dz)):
            continue
        if row.get("candidate_present") != "1":
            continue
        train_rows.append(row)

    if len(train_rows) < args.min_train_rows:
        raise SystemExit(
            f"Too few train rows with targets ({len(train_rows)}), "
            f"need >= {args.min_train_rows}."
        )

    feature_cols = collect_feature_columns(train_rows, deploy_rows)
    if not feature_cols:
        raise SystemExit("No common feature columns found.")

    anchor_names = discover_anchor_names(feature_cols)

    Y_train = np.array(
        [
            [
                to_float(r, "target_delta_x_mm"),
                to_float(r, "target_delta_y_mm"),
                to_float(r, "target_delta_z_mm"),
            ]
            for r in train_rows
        ],
        dtype=float,
    )

    approach_defs: List[Tuple[str, Tuple[str, ...]]] = []
    if not anchor_names:
        approach_defs.append(("candidate_only", tuple()))
    elif args.approach_mode == "combined":
        approach_defs.append(("all_anchors", tuple(anchor_names)))
    else:
        if args.include_candidate_only_in_auto == 1:
            approach_defs.append(("candidate_only", tuple()))
        # Auto mode: all single anchors, all anchor pairs, plus all anchors.
        for r in range(1, len(anchor_names) + 1):
            for combo in itertools.combinations(anchor_names, r):
                approach_defs.append((approach_id_for(combo, anchor_names), combo))

    approach_results: Dict[str, dict] = {}
    for name, anchors in approach_defs:
        approach_results[name] = run_approach(
            args=args,
            approach_name=name,
            approach_anchors=anchors,
            feature_cols=feature_cols,
            train_rows=train_rows,
            deploy_rows=deploy_rows,
            Y_train=Y_train,
        )

    # Approach ranking table.
    comparison_rows: List[dict] = []
    for name, res in approach_results.items():
        comparison_rows.append(
            {
                "approach": name,
                "anchors": ";".join(res["anchors"]),
                "n_features": str(len(res["feature_cols"])),
                "n_gate_features": str(len(res["gate_feature_cols"])),
                "gate_feature_mode": str(res["gate_feature_mode"]),
                "loocv_mean_err_knn_mm": f6(res["loocv_mean_err_knn_mm"]),
                "loocv_mean_err_global_mm": f6(res["loocv_mean_err_global_mm"]),
                "loocv_knn_better_rate": f6(res["loocv_knn_better_rate"]),
                "loocv_mean_direction_agreement_knn": f6(
                    res["loocv_mean_direction_agreement_knn"]
                ),
                "distance_threshold": f6(res["distance_threshold"]),
                "predicted_mag_limit_mm": f6(res["predicted_mag_limit_mm"]),
                "force_knn_if_input_present": str(res["force_knn_if_input_present"]),
                "force_knn_clamp_to_mag_limit": str(res["force_knn_clamp_to_mag_limit"]),
                "deploy_selection_policy": args.deploy_selection_policy,
                "min_expected_gain_mm": f6(args.min_expected_gain_mm),
                "max_selected_mag_mm": f6(args.max_selected_mag_mm),
                "require_positive_direction_agreement": str(
                    args.require_positive_direction_agreement
                ),
                "min_expected_direction_agreement": f6(
                    args.min_expected_direction_agreement
                ),
                "candidate_only_min_expected_gain_mm": f6(
                    args.candidate_only_min_expected_gain_mm
                ),
                "candidate_only_min_confidence": f6(args.candidate_only_min_confidence),
                "uncertainty_confidence_threshold": f6(
                    args.uncertainty_confidence_threshold
                ),
                "fallback_mode": args.fallback_mode,
                "deploy_method_knn_total": str(res["deploy_method_knn_total"]),
                "deploy_method_knn_anchor": str(res["deploy_method_knn_anchor"]),
                "deploy_method_knn_forced": str(res["deploy_method_knn_forced"]),
                "deploy_method_global_fallback": str(res["deploy_method_global_fallback"]),
                "deploy_method_missing_input": str(res["deploy_method_missing_input"]),
            }
        )

    if args.approach_select_metric == "loocv_mean_err_knn_mm":
        global_best_approach = min(
            approach_results.keys(),
            key=lambda k: approach_results[k]["loocv_mean_err_knn_mm"],
        )
    else:
        global_best_approach = max(
            approach_results.keys(),
            key=lambda k: approach_results[k]["loocv_knn_better_rate"],
        )

    # Select final deployment predictions.
    final_pred_rows: List[dict] = []
    final_offsets: Dict[str, dict] = {}
    selection_trace_rows: List[dict] = []

    if args.approach_mode == "combined" or len(approach_results) == 1:
        chosen = approach_results[global_best_approach]
        final_pred_rows = [dict(r) for r in chosen["pred_rows"]]
        for r in final_pred_rows:
            r["selected_approach"] = global_best_approach
        final_offsets = dict(chosen["selected_offsets"])
        eval_by_subject = {e["subject"]: e for e in chosen["deploy_eval"]}
        for row in final_pred_rows:
            sid = row["subject"]
            e = eval_by_subject[sid]
            selection_trace_rows.append(
                {
                    "subject": sid,
                    "global_best_approach": global_best_approach,
                    "selected_approach": global_best_approach,
                    "selection_rule": "combined_single_approach",
                    "knn_candidate_count": "1" if is_knn_method(e["method"]) and e["qc_pass"] == 1 else "0",
                    "fallback_candidate_count": "1"
                    if is_fallback_method(e["method"]) and e["qc_pass"] == 1
                    else "0",
                    "knn_candidate_approaches": global_best_approach
                    if is_knn_method(e["method"]) and e["qc_pass"] == 1
                    else "",
                    "fallback_candidate_approaches": global_best_approach
                    if is_fallback_method(e["method"]) and e["qc_pass"] == 1
                    else "",
                    "selected_method": e["method"],
                    "selected_reason": e["reason"],
                    "selected_confidence": f6(e["conf"]),
                    "selected_distance": f6(e["d0"]),
                    "selected_distance_threshold": f6(e["d0_threshold"]),
                    "selected_feature_zmax": f6(e["z"]),
                    "selected_zmax_gate": f6(e["zmax_gate"]),
                    "selected_knn_mag_mm": f6(e["knn_mag"]),
                    "selected_knn_mag_limit_mm": f6(e["knn_mag_limit"]),
                    "selected_expected_baseline_mag_mm": f6(e["expected_baseline_mag_mm"]),
                    "selected_expected_post_err_mm": f6(e["expected_post_err_mm"]),
                    "selected_expected_gain_mm": f6(e["expected_gain_mm"]),
                    "selected_expected_direction_agreement": f6(
                        e["expected_direction_agreement"]
                    ),
                    "selected_anchor_required": ";".join(e["anchors_required"]),
                    "selected_anchor_present": ";".join(e["anchors_present"]),
                    "selected_anchor_missing": ";".join(e["anchors_missing"]),
                }
            )
    else:
        for i, row in enumerate(deploy_rows):
            sid = row["subject"]
            evals = [res["deploy_eval"][i] for res in approach_results.values()]
            best, selection_rule = choose_best_eval_for_subject(sid, evals, global_best_approach)
            knn_approaches = sorted(
                e["approach"] for e in evals if is_knn_method(e["method"]) and e["qc_pass"] == 1
            )
            fallback_approaches = sorted(
                e["approach"] for e in evals if is_fallback_method(e["method"]) and e["qc_pass"] == 1
            )

            final_row = {
                "subject": sid,
                "selected_approach": best["approach"],
                "method": best["method"],
                "qc_pass": str(best["qc_pass"]),
                "reason": best["reason"],
                "knn_dx_mm": f6(best["pred_knn"][0]),
                "knn_dy_mm": f6(best["pred_knn"][1]),
                "knn_dz_mm": f6(best["pred_knn"][2]),
                "global_dx_mm": f6(best["pred_global"][0]),
                "global_dy_mm": f6(best["pred_global"][1]),
                "global_dz_mm": f6(best["pred_global"][2]),
                "selected_dx_mm": f6(best["pred_sel"][0]),
                "selected_dy_mm": f6(best["pred_sel"][1]),
                "selected_dz_mm": f6(best["pred_sel"][2]),
                "selected_mag_mm": f6(
                    float(np.linalg.norm(best["pred_sel"]))
                    if np.all(np.isfinite(best["pred_sel"]))
                    else math.nan
                ),
                "nearest_feature_distance": f6(best["d0"]),
                "feature_zmax": f6(best["z"]),
                "neighbor_target_spread_mm": f6(best["spread"]),
                "expected_baseline_mag_mm": f6(best["expected_baseline_mag_mm"]),
                "expected_post_err_mm": f6(best["expected_post_err_mm"]),
                "expected_gain_mm": f6(best["expected_gain_mm"]),
                "expected_direction_agreement": f6(best["expected_direction_agreement"]),
                "confidence": f6(best["conf"]),
            }
            final_pred_rows.append(final_row)
            selection_trace_rows.append(
                {
                    "subject": sid,
                    "global_best_approach": global_best_approach,
                    "selected_approach": best["approach"],
                    "selection_rule": selection_rule,
                    "knn_candidate_count": str(len(knn_approaches)),
                    "fallback_candidate_count": str(len(fallback_approaches)),
                    "knn_candidate_approaches": ";".join(knn_approaches),
                    "fallback_candidate_approaches": ";".join(fallback_approaches),
                    "selected_method": best["method"],
                    "selected_reason": best["reason"],
                    "selected_confidence": f6(best["conf"]),
                    "selected_distance": f6(best["d0"]),
                    "selected_distance_threshold": f6(best["d0_threshold"]),
                    "selected_feature_zmax": f6(best["z"]),
                    "selected_zmax_gate": f6(best["zmax_gate"]),
                    "selected_knn_mag_mm": f6(best["knn_mag"]),
                    "selected_knn_mag_limit_mm": f6(best["knn_mag_limit"]),
                    "selected_expected_baseline_mag_mm": f6(best["expected_baseline_mag_mm"]),
                    "selected_expected_post_err_mm": f6(best["expected_post_err_mm"]),
                    "selected_expected_gain_mm": f6(best["expected_gain_mm"]),
                    "selected_expected_direction_agreement": f6(
                        best["expected_direction_agreement"]
                    ),
                    "selected_anchor_required": ";".join(best["anchors_required"]),
                    "selected_anchor_present": ";".join(best["anchors_present"]),
                    "selected_anchor_missing": ";".join(best["anchors_missing"]),
                }
            )

            if best["qc_pass"] == 1 and np.all(np.isfinite(best["pred_sel"])):
                final_offsets[sid] = {
                    "delta_mm": [
                        float(best["pred_sel"][0]),
                        float(best["pred_sel"][1]),
                        float(best["pred_sel"][2]),
                    ],
                    "source": args.offset_source_tag,
                    "approach": best["approach"],
                    "method": best["method"],
                    "confidence": float(best["conf"]) if finite(best["conf"]) else None,
                    "qc": {
                        "passed": int(best["qc_pass"]),
                        "nearest_feature_distance": float(best["d0"]) if finite(best["d0"]) else None,
                        "feature_zmax": float(best["z"]) if finite(best["z"]) else None,
                        "neighbor_target_spread_mm": float(best["spread"]) if finite(best["spread"]) else None,
                        "reason": best["reason"],
                    },
                }

    selected_approach_by_subject = {r["subject"]: r["selected_approach"] for r in final_pred_rows}
    approach_audit_rows: List[dict] = []
    for name, res in approach_results.items():
        for e in res["deploy_eval"]:
            approach_audit_rows.append(
                {
                    "subject": e["subject"],
                    "approach": name,
                    "anchors_required": ";".join(e["anchors_required"]),
                    "anchors_present": ";".join(e["anchors_present"]),
                    "anchors_missing": ";".join(e["anchors_missing"]),
                    "candidate_present_ok": "1" if e["cand_ok"] else "0",
                    "anchor_present_ok": "1" if e["anchor_ok"] else "0",
                    "gate_distance_ok": "1" if e["gate_dist_ok"] else "0",
                    "gate_zmax_ok": "1" if e["gate_zmax_ok"] else "0",
                    "gate_knn_mag_ok": "1" if e["gate_mag_ok"] else "0",
                    "gate_direction_ok": "1" if e["gate_dir_ok"] else "0",
                    "prefer_knn": "1" if e["prefer_knn"] else "0",
                    "method": e["method"],
                    "qc_pass": str(e["qc_pass"]),
                    "reason": e["reason"],
                    "nearest_feature_distance": f6(e["d0"]),
                    "distance_threshold": f6(e["d0_threshold"]),
                    "feature_zmax": f6(e["z"]),
                    "zmax_gate": f6(e["zmax_gate"]),
                    "knn_pred_mag_mm": f6(e["knn_mag"]),
                    "knn_mag_limit_mm": f6(e["knn_mag_limit"]),
                    "expected_baseline_mag_mm": f6(e["expected_baseline_mag_mm"]),
                    "expected_post_err_mm": f6(e["expected_post_err_mm"]),
                    "expected_gain_mm": f6(e["expected_gain_mm"]),
                    "expected_direction_agreement": f6(e["expected_direction_agreement"]),
                    "expected_gain_ok": "1"
                    if (
                        args.deploy_selection_policy != "local_expected_gain"
                        or (
                            finite(e["expected_gain_mm"])
                            and e["expected_gain_mm"] >= args.min_expected_gain_mm
                        )
                    )
                    else "0",
                    "expected_direction_ok": "1"
                    if (
                        args.require_positive_direction_agreement == 0
                        or (
                            finite(e["expected_direction_agreement"])
                            and e["expected_direction_agreement"]
                            > args.min_expected_direction_agreement
                        )
                    )
                    else "0",
                    "neighbor_target_spread_mm": f6(e["spread"]),
                    "confidence": f6(e["conf"]),
                    "selected_in_final": "1" if selected_approach_by_subject.get(e["subject"]) == name else "0",
                }
            )

    # Persist per-approach artifacts.
    approaches_root = args.output_dir / "approaches"
    for name, res in approach_results.items():
        adir = approaches_root / name
        adir.mkdir(parents=True, exist_ok=True)
        loocv_path = adir / "training_loocv.tsv"
        pred_path = adir / "deployment_predictions.tsv"
        model_json_path = adir / "model_spec.json"
        if res["loocv_rows"]:
            write_tsv(loocv_path, res["loocv_rows"], list(res["loocv_rows"][0].keys()))
        if res["pred_rows"]:
            write_tsv(pred_path, res["pred_rows"], list(res["pred_rows"][0].keys()))
        model_spec = {
            "source": args.offset_source_tag,
            "approach": name,
            "anchors": res["anchors"],
            "k": int(args.k),
            "feature_cols": res["feature_cols"],
            "feature_means": res["feature_means"],
            "feature_stds": res["feature_stds"],
            "gate_feature_mode": res["gate_feature_mode"],
            "gate_feature_cols": res["gate_feature_cols"],
            "gate_feature_means": res["gate_feature_means"],
            "gate_feature_stds": res["gate_feature_stds"],
            "global_vec_mm": res["global_vec_mm"],
            "distance_threshold": res["distance_threshold"],
            "zmax_gate": res["zmax_gate"],
            "predicted_mag_limit_mm": res["predicted_mag_limit_mm"],
            "deploy_selection_policy": args.deploy_selection_policy,
            "min_expected_gain_mm": args.min_expected_gain_mm,
            "max_selected_mag_mm": args.max_selected_mag_mm,
            "fallback_mode": args.fallback_mode,
            "loocv_mean_err_knn_mm": res["loocv_mean_err_knn_mm"],
            "loocv_mean_err_global_mm": res["loocv_mean_err_global_mm"],
            "loocv_knn_better_rate": res["loocv_knn_better_rate"],
            "loocv_mean_direction_agreement_knn": res[
                "loocv_mean_direction_agreement_knn"
            ],
            "require_positive_direction_agreement": int(
                args.require_positive_direction_agreement
            ),
            "min_expected_direction_agreement": float(
                args.min_expected_direction_agreement
            ),
        }
        model_json_path.write_text(json.dumps(model_spec, indent=2, sort_keys=True), encoding="utf-8")

    # Top-level artifacts.
    comparison_path = args.output_dir / "approach_comparison.tsv"
    gate_audit_path = args.output_dir / "approach_gate_audit.tsv"
    selection_trace_path = args.output_dir / "approach_selection_trace.tsv"
    pred_path = args.output_dir / "deployment_predictions.tsv"
    offsets_path = args.output_dir / "deployment_offsets_selected.json"
    summary_path = args.output_dir / "model_summary.txt"
    model_json_path = args.output_dir / "model_spec.json"

    if comparison_rows:
        write_tsv(comparison_path, comparison_rows, list(comparison_rows[0].keys()))
    if approach_audit_rows:
        write_tsv(gate_audit_path, approach_audit_rows, list(approach_audit_rows[0].keys()))
    if selection_trace_rows:
        write_tsv(selection_trace_path, selection_trace_rows, list(selection_trace_rows[0].keys()))
    if final_pred_rows:
        write_tsv(pred_path, final_pred_rows, list(final_pred_rows[0].keys()))
    offsets_path.write_text(json.dumps(final_offsets, indent=2, sort_keys=True), encoding="utf-8")

    best_res = approach_results[global_best_approach]
    deploy_knn = sum(1 for r in final_pred_rows if is_knn_method(r["method"]))
    deploy_knn_anchor = sum(1 for r in final_pred_rows if r["method"] == "knn_anchor")
    deploy_knn_forced = sum(1 for r in final_pred_rows if r["method"].startswith("knn_forced"))
    deploy_knn_soft_rescue = sum(
        1 for r in final_pred_rows if r["method"] == "knn_soft_rescue"
    )
    deploy_global = sum(1 for r in final_pred_rows if is_fallback_method(r["method"]))
    deploy_missing = sum(1 for r in final_pred_rows if r["method"] == "missing_input")

    model_spec_top = {
        "source": args.offset_source_tag,
        "approach_mode": args.approach_mode,
        "approach_select_metric": args.approach_select_metric,
        "gate_feature_mode": args.gate_feature_mode,
        "allow_partial_anchor_input": int(args.allow_partial_anchor_input),
        "include_candidate_only_in_auto": int(args.include_candidate_only_in_auto),
        "force_knn_if_input_present": int(args.force_knn_if_input_present),
        "force_knn_clamp_to_mag_limit": int(args.force_knn_clamp_to_mag_limit),
        "deploy_selection_policy": args.deploy_selection_policy,
        "min_expected_gain_mm": float(args.min_expected_gain_mm),
        "max_selected_mag_mm": float(args.max_selected_mag_mm)
        if finite(args.max_selected_mag_mm)
        else None,
        "require_positive_direction_agreement": int(
            args.require_positive_direction_agreement
        ),
        "min_expected_direction_agreement": float(
            args.min_expected_direction_agreement
        ),
        "candidate_only_min_expected_gain_mm": float(
            args.candidate_only_min_expected_gain_mm
        )
        if finite(args.candidate_only_min_expected_gain_mm)
        else None,
        "candidate_only_min_confidence": float(args.candidate_only_min_confidence)
        if finite(args.candidate_only_min_confidence)
        else None,
        "uncertainty_confidence_threshold": float(
            args.uncertainty_confidence_threshold
        )
        if finite(args.uncertainty_confidence_threshold)
        else None,
        "soft_rescue_min_confidence": float(args.soft_rescue_min_confidence)
        if finite(args.soft_rescue_min_confidence)
        else None,
        "soft_rescue_min_expected_gain_mm": float(
            args.soft_rescue_min_expected_gain_mm
        )
        if finite(args.soft_rescue_min_expected_gain_mm)
        else None,
        "soft_rescue_min_direction_agreement": float(
            args.soft_rescue_min_direction_agreement
        )
        if finite(args.soft_rescue_min_direction_agreement)
        else None,
        "soft_rescue_scale": float(args.soft_rescue_scale),
        "fallback_mode": args.fallback_mode,
        "global_best_approach": global_best_approach,
        "approaches": {
            name: {
                "anchors": res["anchors"],
                "n_features": len(res["feature_cols"]),
                "n_gate_features": len(res["gate_feature_cols"]),
                "gate_feature_mode": res["gate_feature_mode"],
                "loocv_mean_err_knn_mm": res["loocv_mean_err_knn_mm"],
                "loocv_mean_err_global_mm": res["loocv_mean_err_global_mm"],
                "loocv_knn_better_rate": res["loocv_knn_better_rate"],
                "loocv_mean_direction_agreement_knn": res[
                    "loocv_mean_direction_agreement_knn"
                ],
                "distance_threshold": res["distance_threshold"],
                "predicted_mag_limit_mm": res["predicted_mag_limit_mm"],
                "force_knn_if_input_present": res["force_knn_if_input_present"],
                "force_knn_clamp_to_mag_limit": res["force_knn_clamp_to_mag_limit"],
                "deploy_selection_policy": args.deploy_selection_policy,
                "min_expected_gain_mm": float(args.min_expected_gain_mm),
                "max_selected_mag_mm": float(args.max_selected_mag_mm)
                if finite(args.max_selected_mag_mm)
                else None,
                "require_positive_direction_agreement": int(
                    args.require_positive_direction_agreement
                ),
                "min_expected_direction_agreement": float(
                    args.min_expected_direction_agreement
                ),
                "candidate_only_min_expected_gain_mm": float(
                    args.candidate_only_min_expected_gain_mm
                )
                if finite(args.candidate_only_min_expected_gain_mm)
                else None,
                "candidate_only_min_confidence": float(
                    args.candidate_only_min_confidence
                )
                if finite(args.candidate_only_min_confidence)
                else None,
                "uncertainty_confidence_threshold": float(
                    args.uncertainty_confidence_threshold
                )
                if finite(args.uncertainty_confidence_threshold)
                else None,
                "soft_rescue_min_confidence": float(args.soft_rescue_min_confidence)
                if finite(args.soft_rescue_min_confidence)
                else None,
                "soft_rescue_min_expected_gain_mm": float(
                    args.soft_rescue_min_expected_gain_mm
                )
                if finite(args.soft_rescue_min_expected_gain_mm)
                else None,
                "soft_rescue_min_direction_agreement": float(
                    args.soft_rescue_min_direction_agreement
                )
                if finite(args.soft_rescue_min_direction_agreement)
                else None,
                "soft_rescue_scale": float(args.soft_rescue_scale),
                "fallback_mode": args.fallback_mode,
            }
            for name, res in approach_results.items()
        },
        "selected_deploy_counts": {
            "knn_total": deploy_knn,
            "knn_anchor": deploy_knn_anchor,
            "knn_forced": deploy_knn_forced,
            "knn_soft_rescue": deploy_knn_soft_rescue,
            "global_fallback": deploy_global,
            "missing_input": deploy_missing,
        },
    }
    model_json_path.write_text(json.dumps(model_spec_top, indent=2, sort_keys=True), encoding="utf-8")

    summary_lines = [
        f"train_dataset\t{args.train_dataset}",
        f"deploy_dataset\t{args.deploy_dataset}",
        f"n_train_rows\t{len(train_rows)}",
        f"n_deploy_rows\t{len(deploy_rows)}",
        f"approach_mode\t{args.approach_mode}",
        f"approach_count\t{len(approach_results)}",
        f"approach_select_metric\t{args.approach_select_metric}",
        f"gate_feature_mode\t{args.gate_feature_mode}",
        f"allow_partial_anchor_input\t{args.allow_partial_anchor_input}",
        f"include_candidate_only_in_auto\t{args.include_candidate_only_in_auto}",
        f"force_knn_if_input_present\t{args.force_knn_if_input_present}",
        f"force_knn_clamp_to_mag_limit\t{args.force_knn_clamp_to_mag_limit}",
        f"deploy_selection_policy\t{args.deploy_selection_policy}",
        f"min_expected_gain_mm\t{args.min_expected_gain_mm:.6f}",
        f"max_selected_mag_mm\t{f6(args.max_selected_mag_mm)}",
        "require_positive_direction_agreement\t"
        f"{args.require_positive_direction_agreement}",
        "min_expected_direction_agreement\t"
        f"{args.min_expected_direction_agreement:.6f}",
        "candidate_only_min_expected_gain_mm\t"
        f"{f6(args.candidate_only_min_expected_gain_mm)}",
        "candidate_only_min_confidence\t"
        f"{f6(args.candidate_only_min_confidence)}",
        "uncertainty_confidence_threshold\t"
        f"{f6(args.uncertainty_confidence_threshold)}",
        "soft_rescue_min_confidence\t"
        f"{f6(args.soft_rescue_min_confidence)}",
        "soft_rescue_min_expected_gain_mm\t"
        f"{f6(args.soft_rescue_min_expected_gain_mm)}",
        "soft_rescue_min_direction_agreement\t"
        f"{f6(args.soft_rescue_min_direction_agreement)}",
        f"soft_rescue_scale\t{args.soft_rescue_scale:.6f}",
        f"fallback_mode\t{args.fallback_mode}",
        f"global_best_approach\t{global_best_approach}",
        f"n_features\t{len(best_res['feature_cols'])}",
        f"k\t{args.k}",
        f"distance_threshold\t{best_res['distance_threshold']:.6f}",
        f"zmax_gate\t{args.zmax_gate:.6f}",
        f"predicted_mag_limit_mm\t{best_res['predicted_mag_limit_mm']:.6f}",
        f"loocv_mean_err_knn_mm\t{best_res['loocv_mean_err_knn_mm']:.6f}",
        f"loocv_mean_err_global_mm\t{best_res['loocv_mean_err_global_mm']:.6f}",
        f"loocv_knn_better_rate\t{best_res['loocv_knn_better_rate']:.6f}",
        "loocv_mean_direction_agreement_knn\t"
        f"{best_res['loocv_mean_direction_agreement_knn']:.6f}",
        f"deploy_method_knn_total\t{deploy_knn}",
        f"deploy_method_knn_anchor\t{deploy_knn_anchor}",
        f"deploy_method_knn_forced\t{deploy_knn_forced}",
        f"deploy_method_knn_soft_rescue\t{deploy_knn_soft_rescue}",
        f"deploy_method_global_fallback\t{deploy_global}",
        f"deploy_method_missing_input\t{deploy_missing}",
        f"approach_comparison_tsv\t{comparison_path}",
        f"approach_gate_audit_tsv\t{gate_audit_path}",
        f"approach_selection_trace_tsv\t{selection_trace_path}",
        f"approaches_dir\t{approaches_root}",
        f"deploy_tsv\t{pred_path}",
        f"offsets_json\t{offsets_path}",
        f"model_spec_json\t{model_json_path}",
    ]
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print(f"WROTE {summary_path}")
    print(f"WROTE {comparison_path}")
    print(f"WROTE {gate_audit_path}")
    print(f"WROTE {selection_trace_path}")
    print(f"WROTE {pred_path}")
    print(f"WROTE {offsets_path} ({len(final_offsets)} subjects)")


if __name__ == "__main__":
    main()
