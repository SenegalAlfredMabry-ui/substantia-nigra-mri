# Substantia Nigra Neuromelanin MRI Pipeline

Neuromelanin MRI lets you visualize the substantia nigra (SN) noninvasively —
no contrast agent, no surgery. The signal comes from the iron-melanin complex
that accumulates in dopaminergic neurons over a lifetime. In Parkinson's
disease, those neurons are the ones dying, which makes neuromelanin MRI a
promising window into dopaminergic integrity across the lifespan. But that
window is only useful if the image processing is tight enough to trust.

This repository contains the pipeline I use to build subject-space SN masks
from neuromelanin-sensitive TSE images, extract quantitative signal metrics,
and iteratively improve mask placement against manual expert tracings. It is
research-pipeline code, built to run reproducibly across a large cohort rather
than as a general-purpose toolbox.

<!-- SCREENSHOT 1 -->
<!-- Take this in Freeview:
     Open a subject's TSE image as the background volume.
     Load the SN group mask as an overlay (any bright color works — red or orange reads well).
     Navigate to an axial slice in the midbrain where the bilateral neuromelanin
     signal is clearly visible as two bright spots. The mask should be sitting
     on top of that signal.
     Screenshot the Freeview window with the slice and overlay visible. -->
![SN mask on neuromelanin TSE](docs/img/sn_mask_on_tse.png)
*SN group mask overlaid on a neuromelanin-sensitive TSE image. The bright signal in the midbrain reflects iron-melanin accumulation in dopaminergic neurons.*

---

## The Problem This Solves

The SN sits in the midbrain, is roughly the size of a small grape, and varies
in position and shape across subjects because of normal differences in brain
geometry, neuromelanin content, and age-related signal change. Getting a
group-level atlas mask to land accurately in each subject's native image space
— without a manual tracing for every single person — is the central challenge
this code addresses.

<!-- SCREENSHOT 2 -->
<!-- Take this in Freeview:
     Open the same TSE image.
     Load two overlays: (1) the KW manual tracing in one color (e.g. green),
     and (2) the group mask in another color (e.g. red).
     Find a subject and slice where there is a visible offset between the two —
     the manual tracing and the group mask are not perfectly overlapping.
     This is the "before improvement" picture. It does not need to look bad,
     just realistic. A small offset is normal and expected.
     Screenshot the Freeview window showing both overlays. -->
![KW tracing vs group mask](docs/img/kw_vs_groupmask.png)
*Manual KW tracing (green) vs. group mask (red) on the same subject. Reducing this offset is the goal of the alignment pipeline.*

---

## What This Repository Contains

The pipeline has two parts that address both sides of the same problem.

**`ogrady_black_ant_refresh/` — Mask building and metric extraction**

Warps a probabilistic SN atlas into each subject's native TSE space, builds
binary and probabilistic masks, and extracts neuromelanin intensity and
contrast-to-noise ratio (CNR) metrics for downstream analysis. This is the
production side of the pipeline — getting outputs to exist reliably.

**`pym_sn_kw_alignment/` — Geometry-guided alignment improvement**

Takes subjects where manual KW tracings exist, scores current mask placement
against those tracings using centroid distance and Dice overlap, then runs a
staged native-space rescue to iteratively reduce the error. This is the
quality-improvement side — making sure the outputs are anatomically correct,
not just present.

---

## Requirements

- **ANTs** (2.4+) — deformable registration
- **AFNI** — affine transforms and image arithmetic
- **FSL** — supplementary mask operations
- **Python 3.8+** — `numpy`, `nibabel`, `scipy`, `pandas`
- **MATLAB R2020a+** — neuromelanin metric extraction
- **R** with `rmarkdown`, `lme4`, `ggplot2` — downstream analysis and reporting
- A neuromelanin-sensitive TSE image and T1-weighted MPRAGE per subject
- A probabilistic SN atlas in MNI152 space (see note below)

**Atlas note:** The internal probabilistic atlas is not redistributed here.
The pipeline expects `sn_prob_atlas_mni152.nii.gz` at the path set by
`BLACKANT_ATLAS_ROOT`. A compatible atlas can be built from published SN
atlases in MNI152_2009 space (e.g., Pauli et al. 2018).

---

## Running the Pipeline

Full cohort runs use SLURM. The core sequence is:

```bash
# Rebuild any missing warp prerequisites
sbatch Scripts/supervillan_wallet.sh

# Generate subject-space SN masks
sbatch Scripts/build_black_ant_masks_v13.sh

# Extract neuromelanin metrics
sbatch Scripts/run_irredeemable_black_ant_V13.sh

# Render downstream analysis report
Rscript -e "rmarkdown::render('Scripts/Thunderbolts_v6.Rmd')"
```

For alignment scoring and rescue, see `pym_sn_kw_alignment/README.md`.

---

## Pipeline Stages

### Stage 1 — Warp prerequisites

Before any masks are built, the pipeline verifies that each subject has a
template-to-TSE warp chain and an anat-to-TSE affine. Missing warps are rebuilt
via `supervillan_wallet.sh`; missing affines via `unlikable_superhero.sh`. This
prevents silent failures downstream where a mask lands in the wrong place
because a stale or missing transform was silently reused.

### Stage 2 — Subject-space mask generation

`build_black_ant_masks_v*.sh` warps the probabilistic SN atlas through the
subject-specific warp chain into native TSE space. It thresholds at a
configurable probability cutoff and writes binary and probabilistic masks
alongside control-region masks — putamen, brainstem, dorsal cochlear nucleus —
needed for CNR estimation.

### Stage 3 — Neuromelanin metric extraction

`irredeemable_black_ant_v*.m` reads the subject-space masks and computes
intensity ratios, CNR, and occupancy metrics for the SN, control regions, and
brainstem. Outputs are per-subject rows that merge into a cohort-level table
for statistical analysis.

<!-- SCREENSHOT 3 -->
<!-- Open a terminal and run check_sn_alignment_metrics.py on a small subset,
     or just open the output TSV in a spreadsheet or display the first ~10 rows
     in the terminal with `column -t`. The columns to show are: subject ID,
     Dice, delta_x_mm, delta_y_mm, delta_z_mm, delta_mag_mm.
     Screenshot the terminal or spreadsheet view with those columns visible. -->
![Alignment scoring output](docs/img/alignment_scores.png)
*Per-subject alignment scores: Dice overlap and centroid distance in mm relative to manual KW tracings.*

### Stage 4 — Alignment scoring and rescue *(active development)*

For subjects where manual KW tracings are available,
`check_sn_alignment_metrics.py` scores current mask placement using Dice
overlap and centroid distance. The rescue pipeline in `pym_sn_kw_alignment/`
then runs staged correction:

1. **Scalar translation sweep** — applies a gain ladder to the KW-vs-pre
   centroid delta and accepts the best stage per subject
2. **Per-axis micro-gain rescue** — tests signed per-axis gains independently
   for subjects that did not improve under the scalar sweep
3. **Residual refinement** — targeted follow-up on remaining high-delta
   subjects using a finer local search

Progress across the 108-subject KW-tracing cohort:

| Stage | Subjects improved |
|---|---|
| Scalar translation sweep | 67 / 108 |
| Per-axis gain rescue | 41 / 41 remaining |
| Residual refinement | in progress |

The next step is extending this approach to subjects without manual tracings
using an anchor-geometry offset model trained on the KW-known cohort.

---

## Downstream Analysis

Once the masks and metrics are staged, `Thunderbolts_v*.Rmd` runs mixed-effects
models of neuromelanin signal against age, cardiovascular, and autonomic
variables.

<!-- SCREENSHOT 4 -->
<!-- Open a rendered Thunderbolts HTML report (any version, e.g.
     Scripts/SN_ASR_DRAFT_analysis.html) in a browser.
     Screenshot the first meaningful plot or table — a regression output,
     a CNR-by-group comparison, or an age curve.
     This shows what the pipeline outputs look like once they feed analysis. -->
![Downstream analysis output](docs/img/thunderbolts_output.png)
*Example downstream output: neuromelanin CNR across the cohort entering mixed-effects analysis.*

---

## Scientific Context

This pipeline supports the SNHEART study, a lifespan neuroimaging study of
cardiac and autonomic function. The SN is one of two target nuclei — the locus
coeruleus (LC) is the other — chosen because both carry neuromelanin signal and
both are implicated in early autonomic and dopaminergic dysfunction in
Parkinson's disease.

---

## Repository Layout

```
ogrady_black_ant_refresh/     Mask building and neuromelanin metric extraction
  scripts/                    Shell, Python, and MATLAB scripts
  README.md                   Full workflow guide

pym_sn_kw_alignment/          Native-space SN-to-KW alignment and rescue
  scripts/                    Alignment scoring and rescue scripts
  README.md                   Full workflow guide
```

---

## Development History

This repository reflects the working pipeline, updated as the research
progresses. The alignment rescue work in `pym_sn_kw_alignment/` is actively
evolving. Commit history tracks real development decisions, not a
post-hoc cleanup.

---

## Reproducibility Notes

Path conventions follow the internal study structure and will need adjustment
for external use. Subject lists and atlas files are not redistributed.
Adapting the pipeline to a new dataset means updating path and environment
configuration; the logic of each stage is documented in the per-module READMEs.
