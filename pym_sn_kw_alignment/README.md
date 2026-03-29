# SN-to-KW Alignment Kit

Source guide: `Scripts/PYM_TECH/PYM/Astonishing_Field_Guide.txt`

This bundle is a cleaned public copy of the Big House alignment workflow. Its
purpose is to improve the placement of native-space substantia nigra masks by
comparing them against raw KW manual tracings, measuring the mismatch, and then
running targeted correction and rescue strategies.

This is the part of the project that is most explicitly geometric. The code is
less about producing any output at all and more about asking whether a proposed
alignment actually lands the mask closer to anatomically traced truth.

## Scope

The workflow is organized around five stages:

1. Truth-source preparation and geometric scoring.
2. Upstream warp refresh and ROI steering.
3. Native subject-space mask rebuilding.
4. Raw-KW correction and native rescue sweeps.
5. Fallback metric packaging and downstream reporting.

## How The Workflow Works

The code treats manual KW tracings as the geometric reference point. Candidate
native masks are scored against that reference using centroid distances, Dice,
and related summaries. Those measurements are then used to decide what to try
next: refresh upstream warps, shift ROI targets, rebuild masks, or run local
native-space rescue sweeps.

That leads to a repeated loop:

1. stage truth and candidate masks
2. score the geometric gap
3. run a correction strategy
4. rescore and compare against baseline
5. keep only the approaches that actually improve the geometry

## Representative Entry Points

These files give a good first pass through the project:

- `scripts/check_sn_alignment_metrics.py`
  canonical scorer for candidate-versus-manual geometry
- `scripts/submit_antpunk_chunks.sh`
  batch launcher for upstream warp refresh experiments
- `scripts/antpunk_v1.sh`
  warp-refinement and ROI-guided alignment runner
- `scripts/build_black_ant_masks_v10.sh`
  native mask emitter used after upstream warps improve
- `scripts/run_raw_kw_correction_loop.sh`
  main correction loop tying truth staging, scoring, and rebuild jobs together
- `scripts/prepare_native_sn_rigid_rescue_inputs.py`
  handoff from warp-space experiments into native-space rescues
- `scripts/launch_native_sn_rigid_rescue_chunks.sh`
  example of chunked native rescue orchestration

## What This Code Shows

- building a measurement loop around a research question instead of relying on
  visual inspection alone
- using batch infrastructure to test many alignment strategies systematically
- separating scoring, routing, correction, and rescue into distinct script
  roles
- preserving exploratory personality and naming while still keeping the
  workflow structured enough to revisit later
- treating failed or underperforming runs as useful information for the next
  iteration rather than dead ends

## Included Scripts

The `scripts/` folder contains 38 workspace-local files referenced by the
source guide, including:

- scoring and routing tools such as `check_sn_alignment_metrics.py`,
  `measure_kw_alignment.py`, `summarize_multiopt_by_subject.py`, and
  `hard16_scoreboard_and_routing.py`
- warp and ROI steering scripts such as `submit_antpunk_chunks.sh`,
  `antpunk_v1.sh`, `shift_roi_to_target.py`, and
  `highhat_prepare_template_targets.py`
- native mask rebuild and prerequisite repair scripts such as
  `build_black_ant_masks_v10.sh`, `supervillan_wallet.sh`, and
  `unlikable_superhero.sh`
- rescue and correction runners such as `run_raw_kw_correction_loop.sh`,
  `prepare_native_sn_rigid_rescue_inputs.py`,
  `launch_native_sn_rigid_rescue_chunks.sh`,
  `launch_native_sn_translation_gain_chunks.sh`, and the associated chunk
  runners
- fallback reporting files such as
  `launch_irredeemable_black_ant_V12_kwsn_chunks.sh`,
  `run_irredeemable_black_ant_V12_kwsn.sh`, and
  `Thunderbolts_v5_CPonly.Rmd`

## Not Bundled

This source guide explicitly marks `compute_virgo_target.py` as retired and not
part of the active stack. It was therefore not copied into this bundle.

## Upload Notes

These copies were sanitized for public release. The script names and informal
tone were preserved, while lab-specific paths, real email addresses, and a
small number of hard-coded cohort defaults were replaced with placeholders or
example values.
