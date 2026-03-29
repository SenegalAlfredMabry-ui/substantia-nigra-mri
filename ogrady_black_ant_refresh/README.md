# Black ANT Refresh Kit

Source guide: `Scripts/PYM_TECH/OGRADY/Irredeemable_Field_Guide.txt`

This bundle is a cleaned public copy of the O'Grady Black ANT workflow. Its
purpose is to regenerate the analysis inputs needed for substantia nigra
neuromelanin work when upstream pieces change or when missing subjects need to
be repaired.

In practice, that means:

- recovering missing prerequisites
- rebuilding atlas-derived masks in subject space
- regenerating corrected TSE outputs and metric tables
- handing those refreshed outputs to downstream analysis/reporting code

## Scope

The workflow is organized around five stages:

1. Cleanup prerequisites for late or missing cases.
2. Atlas, warp, and affine preparation.
3. Subject-space mask rebuilding.
4. Neuromelanin metric generation.
5. QC and reporting.

## How The Workflow Works

The core logic is spatial and dependency-driven.

First, the code makes sure each subject has the required upstream products:
reconstructions, anatomical references, warp fields, and anat-to-TSE affine
matrices. Once those prerequisites exist, probabilistic ROIs can be moved from
template space back into each subject's native TSE space. Those subject-space
masks then drive intensity and CNR extraction from corrected and raw
neuromelanin volumes, producing tables that feed statistical analysis and
report generation.

The overall pattern is:

1. repair missing inputs
2. rebuild mask geometry
3. regenerate neuromelanin metrics
4. compare, QC, and pass outputs to analysis notebooks

## Representative Entry Points

If a reviewer wants to understand the workflow quickly, these are good anchor
files:

- `scripts/build_black_ant_masks.sh`
  baseline subject-space mask generator
- `scripts/build_black_ant_masks_v7.sh`
  later-stage evolved builder with additional safeguards and tuning
- `scripts/supervillan_wallet.sh`
  warp-repair wrapper that restores missing prerequisite transforms
- `scripts/unlikable_superhero.sh`
  anat-to-TSE affine regeneration wrapper
- `scripts/black_ant_SN_ASR_metrics.m`
  neuromelanin metric extraction stage
- `scripts/Analysis_Generator.sh`
  handoff into analysis/report generation

## What This Code Shows

- orchestration of multi-stage imaging dependencies across Bash, MATLAB, and R
- practical recovery tooling for incomplete or partially failed cohorts
- use of explicit file conventions and manifests to keep batch outputs
  interpretable
- iteration across script versions rather than pretending the first pipeline
  design was final
- a research workflow that treats data repair, reproducibility, and analysis
  handoff as part of the same engineering problem

## Included Scripts

The `scripts/` folder contains 42 workspace-local files referenced by the
source guide, including:

- cleanup runners such as `Anne_Marie_Hoag_SN_CLEANUP_1.sh`,
  `Lenny_Ballinger_SN_CLEANUP_2.sh`, `Visioneer_SN_CLEANUP_3.sh`, and
  `run_black_goliath_SN_CLEANUP_4.sh`
- atlas and warp builders such as `S.H.I.E.L.D_SN_Warp_v1.sh`,
  `supervillan_wallet.sh`, `unlikable_superhero.sh`, and `warper.sh`
- mask builders from `build_black_ant_masks.sh` through
  `build_black_ant_masks_v7.sh`
- metric and analysis scripts such as `irredeemable_black_ant_V1.m`,
  `black_ant_SN_ASR_metrics.m`, `Analysis_Generator.sh`, and
  `SN_ASR_DRAFT_analysis.Rmd`
- helper utilities such as `compare_black_ant_versions.py`,
  `check_sn_metrics.py`, `supervillain_chunk_status.py`, and
  `pt_mask_ratio_threshold.py`

## Not Bundled

These references appeared in the source guide but were not copied here:

- `SegmentAAN.sh`: external dependency, not present in this workspace
- `arousalNetworkLabels.v10.m`: referenced output/dependency, not present as a
  local script
- `..._cponly.sh` and `..._vSHOWER.sh`: shorthand placeholders in the guide,
  not literal filenames

## Upload Notes

These copies were sanitized for public release. The fun script names and
commentary were preserved, but lab-specific paths, real email addresses, and
hard-coded cohort defaults were replaced with placeholders or example values.
