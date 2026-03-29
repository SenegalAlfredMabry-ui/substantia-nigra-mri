#!/usr/bin/env bash
set -euo pipefail

# Build and train the no-KW offset policy scaffold (kNN + anchor fallback).
#
# Defaults are wired to the current 108-subject refinement outputs for training.
# Override via environment variables as needed.

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

timestamp="$(date -u +%Y%m%dT%H%M%SZ)"
out_root="${NO_KW_POLICY_OUT_ROOT:-${repo_root}/tmp/no_kw_offset_anchor_policy_${timestamp}}"
mkdir -p "${out_root}"

train_subjects="${NO_KW_POLICY_TRAIN_SUBJECTS:-${repo_root}/tmp/native_sn_axis_gain_refine_all108_20260303T185058Z/subjects_all108.txt}"
train_candidate_dir="${NO_KW_POLICY_TRAIN_CANDIDATE_DIR:-${repo_root}/tmp/native_sn_axis_gain_refine_all108_20260303T185058Z/current_best_masks}"
default_train_manual_dir="${repo_root}/tmp/native_sn_axis_gain_20260302T144546Z/inputs/manual_masks"
train_manual_dir="${NO_KW_POLICY_TRAIN_MANUAL_DIR:-${default_train_manual_dir}}"
train_summary_tsv="${NO_KW_POLICY_TRAIN_SUMMARY_TSV:-${repo_root}/tmp/native_sn_axis_gain_refine_all108_20260303T185058Z/current_best_summary_all108.tsv}"

deploy_subjects="${NO_KW_POLICY_DEPLOY_SUBJECTS:-${repo_root}/Scripts/MOTIP_potential_roster.txt}"
deploy_candidate_dir="${NO_KW_POLICY_DEPLOY_CANDIDATE_DIR:-${repo_root}/MRI_data/TSE/probabilistic_masks_v10/base}"

anchor_dir="${NO_KW_POLICY_ANCHOR_DIR:-${repo_root}/MRI_data/TSE/brainstem_masks_v9}"
anchor_prefixes="${NO_KW_POLICY_ANCHOR_PREFIXES:-DC_mask_inTSE brainstem_mask_inTSE brainstem_4v_mask_inTSE}"

k_neighbors="${NO_KW_POLICY_K:-9}"
zmax_gate="${NO_KW_POLICY_ZMAX_GATE:-5.0}"
dist_q="${NO_KW_POLICY_DISTANCE_QUANTILE:-0.93}"
mag_scale="${NO_KW_POLICY_MAG_LIMIT_SCALE:-1.5}"
min_train_rows="${NO_KW_POLICY_MIN_TRAIN_ROWS:-20}"
approach_mode="${NO_KW_POLICY_APPROACH_MODE:-auto}"
approach_select_metric="${NO_KW_POLICY_APPROACH_SELECT_METRIC:-loocv_mean_err_knn_mm}"
gate_feature_mode="${NO_KW_POLICY_GATE_FEATURE_MODE:-anchor_relative}"
force_knn_if_input_present="${NO_KW_POLICY_FORCE_KNN_IF_INPUT_PRESENT:-0}"
force_knn_clamp_to_mag_limit="${NO_KW_POLICY_FORCE_KNN_CLAMP_TO_MAG_LIMIT:-1}"
allow_partial_anchor_input="${NO_KW_POLICY_ALLOW_PARTIAL_ANCHOR_INPUT:-1}"
include_candidate_only_in_auto="${NO_KW_POLICY_INCLUDE_CANDIDATE_ONLY_IN_AUTO:-1}"
deploy_selection_policy="${NO_KW_POLICY_DEPLOY_SELECTION_POLICY:-local_expected_gain}"
min_expected_gain_mm="${NO_KW_POLICY_MIN_EXPECTED_GAIN_MM:-0.35}"
max_selected_mag_mm="${NO_KW_POLICY_MAX_SELECTED_MAG_MM:-nan}"
fallback_mode="${NO_KW_POLICY_FALLBACK_MODE:-no_correction}"
require_positive_direction_agreement="${NO_KW_POLICY_REQUIRE_POSITIVE_DIRECTION_AGREEMENT:-1}"
min_expected_direction_agreement="${NO_KW_POLICY_MIN_EXPECTED_DIRECTION_AGREEMENT:-0.0}"
candidate_only_min_expected_gain_mm="${NO_KW_POLICY_CANDIDATE_ONLY_MIN_EXPECTED_GAIN_MM:-0.75}"
candidate_only_min_confidence="${NO_KW_POLICY_CANDIDATE_ONLY_MIN_CONFIDENCE:-0.60}"
uncertainty_confidence_threshold="${NO_KW_POLICY_UNCERTAINTY_CONFIDENCE_THRESHOLD:-0.54}"
soft_rescue_min_confidence="${NO_KW_POLICY_SOFT_RESCUE_MIN_CONFIDENCE:-nan}"
soft_rescue_min_expected_gain_mm="${NO_KW_POLICY_SOFT_RESCUE_MIN_EXPECTED_GAIN_MM:-nan}"
soft_rescue_min_direction_agreement="${NO_KW_POLICY_SOFT_RESCUE_MIN_DIRECTION_AGREEMENT:-nan}"
soft_rescue_scale="${NO_KW_POLICY_SOFT_RESCUE_SCALE:-1.0}"
candidate_threshold="${NO_KW_POLICY_CANDIDATE_THRESHOLD:-0.5}"
manual_threshold="${NO_KW_POLICY_MANUAL_THRESHOLD:-0.5}"
anchor_threshold="${NO_KW_POLICY_ANCHOR_THRESHOLD:-0.5}"
augment_train_with_perturb="${NO_KW_POLICY_AUGMENT_TRAIN_WITH_PERTURB:-0}"
augment_preset="${NO_KW_POLICY_AUGMENT_PRESET:-nonkw_rescue_v1}"
augment_vectors_mm="${NO_KW_POLICY_AUGMENT_VECTORS_MM:-}"
augment_keep_original="${NO_KW_POLICY_AUGMENT_KEEP_ORIGINAL:-1}"

count_manual_hits() {
  local manual_dir="$1"
  local subjects_file="$2"
  local hit=0
  local sid=""
  [[ -d "${manual_dir}" ]] || { echo 0; return; }
  while IFS= read -r line || [[ -n "${line}" ]]; do
    sid="$(awk '{print $1}' <<< "${line}")"
    [[ -n "${sid}" ]] || continue
    if [[ -e "${manual_dir}/sn_groupmask_in_${sid}.nii.gz" || -e "${manual_dir}/sn_groupmask_in_${sid}.nii" ]]; then
      hit=$((hit + 1))
    fi
  done < "${subjects_file}"
  echo "${hit}"
}

manual_hits="$(count_manual_hits "${train_manual_dir}" "${train_subjects}")"
manual_dir_locked=0
if [[ -n "${NO_KW_POLICY_TRAIN_MANUAL_DIR:-}" ]]; then
  manual_dir_locked=1
fi

if (( manual_hits < min_train_rows )) && (( manual_dir_locked == 0 )); then
  for fallback_dir in \
    "${repo_root}/tmp/native_sn_axis_gain_20260302T144546Z/inputs/manual_masks_20260312_fix" \
    "${repo_root}/tmp/native_sn_translation_gain_20260227T084336Z/inputs/manual_masks" \
    "${repo_root}/tmp/native_sn_rigid_rescue_20260227T082102Z/inputs/manual_masks"; do
    fallback_hits="$(count_manual_hits "${fallback_dir}" "${train_subjects}")"
    if (( fallback_hits >= min_train_rows )); then
      echo "[WARN] train_manual_dir has only ${manual_hits} usable masks; switching to ${fallback_dir} (${fallback_hits} usable masks)."
      train_manual_dir="${fallback_dir}"
      manual_hits="${fallback_hits}"
      break
    fi
  done
fi

if (( manual_hits < min_train_rows )); then
  echo "[ERROR] train_manual_dir=${train_manual_dir} provides only ${manual_hits} usable masks for ${train_subjects}; need >= ${min_train_rows}." >&2
  echo "[ERROR] Set NO_KW_POLICY_TRAIN_MANUAL_DIR to a valid directory with sn_groupmask_in_<ID>.nii[.gz] files." >&2
  exit 1
fi

train_tsv_base="${out_root}/train_dataset.tsv"
train_tsv="${train_tsv_base}"
deploy_tsv="${out_root}/deploy_dataset.tsv"
model_out="${out_root}/model"

echo "[INFO] out_root=${out_root}"
echo "[INFO] train_subjects=${train_subjects}"
echo "[INFO] train_candidate_dir=${train_candidate_dir}"
echo "[INFO] train_manual_dir=${train_manual_dir}"
echo "[INFO] train_manual_overlap=${manual_hits}"
echo "[INFO] deploy_subjects=${deploy_subjects}"
echo "[INFO] deploy_candidate_dir=${deploy_candidate_dir}"
echo "[INFO] anchor_dir=${anchor_dir}"

python3 Scripts/build_sn_no_kw_dataset.py \
  --subjects-file "${train_subjects}" \
  --candidate-dir "${train_candidate_dir}" \
  --manual-dir "${train_manual_dir}" \
  --summary-tsv "${train_summary_tsv}" \
  --candidate-threshold "${candidate_threshold}" \
  --manual-threshold "${manual_threshold}" \
  --anchor-threshold "${anchor_threshold}" \
  --anchor-dir "${anchor_dir}" \
  --anchor-prefixes ${anchor_prefixes} \
  --output-tsv "${train_tsv_base}"

python3 Scripts/build_sn_no_kw_dataset.py \
  --subjects-file "${deploy_subjects}" \
  --candidate-dir "${deploy_candidate_dir}" \
  --candidate-threshold "${candidate_threshold}" \
  --anchor-threshold "${anchor_threshold}" \
  --anchor-dir "${anchor_dir}" \
  --anchor-prefixes ${anchor_prefixes} \
  --output-tsv "${deploy_tsv}"

if [[ "${augment_train_with_perturb}" == "1" ]]; then
  train_tsv="${out_root}/train_dataset_augmented.tsv"
  augment_cmd=(
    python3 Scripts/augment_sn_no_kw_dataset_with_perturbations.py
    --input-tsv "${train_tsv_base}"
    --output-tsv "${train_tsv}"
    --preset "${augment_preset}"
    --keep-original "${augment_keep_original}"
  )
  if [[ -n "${augment_vectors_mm}" ]]; then
    augment_cmd+=(--perturb-vectors-mm "${augment_vectors_mm}")
  fi
  "${augment_cmd[@]}"
fi

python3 Scripts/train_sn_no_kw_offset_policy.py \
  --train-dataset "${train_tsv}" \
  --deploy-dataset "${deploy_tsv}" \
  --output-dir "${model_out}" \
  --k "${k_neighbors}" \
  --zmax-gate "${zmax_gate}" \
  --distance-quantile "${dist_q}" \
  --mag-limit-scale "${mag_scale}" \
  --approach-mode "${approach_mode}" \
  --approach-select-metric "${approach_select_metric}" \
  --gate-feature-mode "${gate_feature_mode}" \
  --force-knn-if-input-present "${force_knn_if_input_present}" \
  --force-knn-clamp-to-mag-limit "${force_knn_clamp_to_mag_limit}" \
  --allow-partial-anchor-input "${allow_partial_anchor_input}" \
  --include-candidate-only-in-auto "${include_candidate_only_in_auto}" \
  --deploy-selection-policy "${deploy_selection_policy}" \
  --min-expected-gain-mm "${min_expected_gain_mm}" \
  --max-selected-mag-mm "${max_selected_mag_mm}" \
  --require-positive-direction-agreement "${require_positive_direction_agreement}" \
  --min-expected-direction-agreement "${min_expected_direction_agreement}" \
  --candidate-only-min-expected-gain-mm "${candidate_only_min_expected_gain_mm}" \
  --candidate-only-min-confidence "${candidate_only_min_confidence}" \
  --uncertainty-confidence-threshold "${uncertainty_confidence_threshold}" \
  --soft-rescue-min-confidence "${soft_rescue_min_confidence}" \
  --soft-rescue-min-expected-gain-mm "${soft_rescue_min_expected_gain_mm}" \
  --soft-rescue-min-direction-agreement "${soft_rescue_min_direction_agreement}" \
  --soft-rescue-scale "${soft_rescue_scale}" \
  --fallback-mode "${fallback_mode}" \
  --min-train-rows "${min_train_rows}"

cat > "${out_root}/run_manifest.txt" <<EOF
timestamp_utc	${timestamp}
train_subjects	${train_subjects}
train_candidate_dir	${train_candidate_dir}
train_manual_dir	${train_manual_dir}
train_summary_tsv	${train_summary_tsv}
deploy_subjects	${deploy_subjects}
deploy_candidate_dir	${deploy_candidate_dir}
anchor_dir	${anchor_dir}
anchor_prefixes	${anchor_prefixes}
k_neighbors	${k_neighbors}
zmax_gate	${zmax_gate}
distance_quantile	${dist_q}
mag_limit_scale	${mag_scale}
min_train_rows	${min_train_rows}
approach_mode	${approach_mode}
approach_select_metric	${approach_select_metric}
gate_feature_mode	${gate_feature_mode}
force_knn_if_input_present	${force_knn_if_input_present}
force_knn_clamp_to_mag_limit	${force_knn_clamp_to_mag_limit}
allow_partial_anchor_input	${allow_partial_anchor_input}
include_candidate_only_in_auto	${include_candidate_only_in_auto}
deploy_selection_policy	${deploy_selection_policy}
min_expected_gain_mm	${min_expected_gain_mm}
max_selected_mag_mm	${max_selected_mag_mm}
require_positive_direction_agreement	${require_positive_direction_agreement}
min_expected_direction_agreement	${min_expected_direction_agreement}
candidate_only_min_expected_gain_mm	${candidate_only_min_expected_gain_mm}
candidate_only_min_confidence	${candidate_only_min_confidence}
uncertainty_confidence_threshold	${uncertainty_confidence_threshold}
soft_rescue_min_confidence	${soft_rescue_min_confidence}
soft_rescue_min_expected_gain_mm	${soft_rescue_min_expected_gain_mm}
soft_rescue_min_direction_agreement	${soft_rescue_min_direction_agreement}
soft_rescue_scale	${soft_rescue_scale}
fallback_mode	${fallback_mode}
candidate_threshold	${candidate_threshold}
manual_threshold	${manual_threshold}
anchor_threshold	${anchor_threshold}
augment_train_with_perturb	${augment_train_with_perturb}
augment_preset	${augment_preset}
augment_vectors_mm	${augment_vectors_mm}
augment_keep_original	${augment_keep_original}
train_dataset_base_tsv	${train_tsv_base}
train_dataset_tsv	${train_tsv}
deploy_dataset_tsv	${deploy_tsv}
model_dir	${model_out}
EOF

echo "[INFO] Wrote run manifest: ${out_root}/run_manifest.txt"
echo "[INFO] Done."
