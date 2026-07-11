# btaf530 full-result verification audit

Date: 2026-07-09

Paper: `/Users/licheng/Documents/文章/btaf530.pdf`

Question: can the current local workspace verify all experimental results in Zhao et al.,
Bioinformatics 2025, btaf530?

## Bottom line

No. The current workspace verifies several local DiSTect workflows, but it does not yet
verify all btaf530 paper results under the paper's original experimental design.

The strongest verified evidence is the local independent re-run in `repro/_verify/REPORT.md`.
That re-run confirms the local verification kit is internally reproducible. It does not cover
the full paper result set because the original full HER2 ADVI analysis, baseline comparisons,
and the STARmap PLUS Alzheimer data/protocol still need to be reproduced against the newly
obtained supplementary targets.

## Source-material status

| Source | Status | Notes |
|---|---|---|
| Main PDF | Available | `/Users/licheng/Documents/文章/btaf530.pdf`, 9 pages. |
| OUP supplementary data | Available via Europe PMC mirror | OUP `academic.oup.com` still returns a Cloudflare challenge from CLI, but Europe PMC `PMC12502917/supplementaryFiles` returned the official zip. Local files: `repro/_verify/supplementary_fetch/btaf530_supplementaryFiles.zip`, extracted PDF `repro/_verify/supplementary_fetch/btaf530_supplementaryFiles/btaf530_supplementary_data.pdf`, extracted text `repro/_verify/supplementary_fetch/btaf530_supplementary_data.txt`. |
| PMC full-text XML | Available | `repro/_verify/supplementary_fetch/epmc_fullTextXML.xml`; confirms `btaf530_supplementary_data.pdf` as the associated supplementary material. |
| Zenodo record 17127211 | Available | Contains `DiSTect-main.zip` only: package source + tutorial. No supplementary result tables. |
| Official package source | Available | Downloaded to `/tmp/distect_zenodo`; functions are `dsgd_single`, `dsgd_multiple`, `missing_imputation`. |
| Local workspace | Available | `/Users/licheng/Documents/DiSTect`; includes additional `repro/` verification and upgraded-method code not in Zenodo package. |

## Supplement-derived target map

The official supplement is now local and text-extracted. Key targets found:

- Section A2.1 / Simulation 1: 30 x 30 lattice, 20 covariates, 2000 Gibbs iterations,
  200 simulation iterations, beta `(1, 2, 3, -4, -5, zero_15)`,
  eta `{0.4, 1.6, 2.8}`, `c1=8`, `b1=5`, `b2=50`, `v0=0.000001`.
  Targets: Fig. B1-B6 and Tables C1-C3.
- Section A2.2 / Simulation 2: six slices, 900 spots per slice, eta `1.6`,
  cross-slice settings `{sigma, rho}` = `{0.1,0.1}`, `{0.1,0.4}`, `{0.4,0.1}`,
  `{0.4,0.4}`, 200 ADVI simulations. Targets: Fig. B7-B8 and Table C4.
- Section A2.3 / Simulation 3: same covariates/responses as A2.1, ignorable missingness
  with 10 or 30 masked spots, nonignorable missingness with gamma `(-6,1,4)` and
  `(-5,1,1.6)`, 200 simulations. Targets: Fig. B9-B10 and Tables C5-C6.
- Section A3 / STARmap PLUS Alzheimer: 13-month disease samples, Replicate 1 training
  and Replicate 2 validation, disease label is within 64 pixels of plaque center, neighbor
  radius 60 pixels. Key targets include Cst7 beta `1.580 [1.551, 1.609]`, Trem2 beta
  `0.992 [0.957, 1.027]`, C1qa beta `0.726 [0.693, 0.759]`, Gfap beta
  `0.577 [0.527, 0.627]`, Hexb interaction beta `0.022 [0.011, 0.033]`, and
  Replicate 2 accuracy `97.1%`. Target: Fig. B12.
- Table C7: HER2 hyperparameters are `b1=5`, `b2=50`, `b3=5`, `b4=50`, `b5=0`,
  `b6=10`, `c1=8`, `v0=0.000001`.

## Verification matrix

| Paper result block | Paper target | Local evidence | Status | Gap |
|---|---|---|---|---|
| Workflow / toy package demo | Basic single-slice, multi-slice, prediction, missing imputation | `repro/reproduce.R`; `_verify/reproduce.log` rc=0 | Partial | Demonstrates package mechanics, not all paper simulation grids. |
| Simulation 1 | Bias, SD, coverage across spatial-correlation settings; Fig. B2-B6, Tables C1-C3 | `repro/paper_reproduce/run_sim1.R` completed `N_REP=200` ADVI pass | Not fully verified | Directional proposed-vs-naive result reproduced, but numeric comparison fails, especially eta bias and coverage. See `repro/paper_reproduce/output/sim1_advi_key_comparison_n200.csv`. |
| Simulation 2 | Multi-slice robustness; Fig. B7-B8, Table C4 | `repro/paper_reproduce/run_sim2.R` completed `N_REP=200` ADVI pass | Not fully verified | Numeric comparison fails; beta biases are lower than target and eta bias is higher. See `repro/paper_reproduce/output/sim2_advi_key_comparison_n200.csv`. |
| Simulation 3 | Missing-data inference under ignorable/non-ignorable missingness; Fig. B9-B10, Tables C5-C6 | `repro/paper_reproduce/run_sim3.R` completed `N_REP=200` missing first-pass ADVI pass | Not fully verified | Numeric comparison fails; missing-mechanism convention remains unresolved for nonignorable settings. See `repro/paper_reproduce/output/sim3_advi_key_comparison_n200.csv`. |
| Eta-bias differential diagnosis | Explain why eta bias is systematically higher than paper targets | GLM pseudolikelihood, generation-grid, random-scan Gibbs, ADVI algorithm, intercept, and centered-neighbor diagnostics completed | Mismatch localized, not resolved | The eta gap is not caused by ADVI alone. It persists in plain GLM and single-site Gibbs. The paper uses neighbor average while the package uses neighbor sum, but the sum convention does not recover all eta targets. See `repro/paper_reproduce/ETA_BIAS_DIAGNOSIS.md`. |
| Original simulation script search | Find authors' source scripts for Simulation 1/2/3 | Searched paper-declared GitHub/Zenodo, official fork, visible repo history, tutorial page, issues/PRs/releases/wiki, author/org public repos, and Zenodo keywords | Public script not found | Public resources contain package source, toy data, and tutorials only. See `repro/paper_reproduce/ORIGINAL_SIM_SCRIPT_SEARCH.md`. |
| HER2 data availability | 8 patients; A-D 6 slices, E-H 3 slices; average 346 spots | Current `/tmp/her2st` has counts for 36 slices in git, but local checkout only has labels for A1/B1/C1/D1/E1/F1/G2/H1 | Blocked/partial | Need disease labels for every slice used by the paper, or confirm paper only used labeled sections despite text. |
| HER2 original ADVI | Fig. 2a-i; 30 top genes; eta=2.936 [2.926, 2.946]; runtime 25h32m | `repro/her2_reproduce.R` is explicitly scaled down to 4 patients, one section each, 40 HVG; `_verify/her2_reproduce.log` gives eta=2.665 and different top genes | Not fully verified | Need original full-slice/full-gene ADVI run. Expected runtime is ~25h on M3/24GB per paper. |
| HER2 all-8 local upgraded method | Full 8 patients, PG-CAVI, NB covariates, FDR, calibrated LOPO | `repro/realdata/run_her2st.R`; `repro/realdata/her2_new.log` | Verified as local extension | Useful extension, but not paper-original result because method/covariates/selection differ. |
| HER2 baselines | Giotto, Celina, MERINGUE, PROST, SpatialDE overlap and prediction comparison | No local runnable baseline scripts or logs found | Not verified | Need baseline implementations, selected genes, and prediction protocol from supplement. |
| HER2 deterministic/probabilistic prediction | Fig. 2f-h, within-patient 3-fold and leave-one-patient-out | Local scaled LOPO A/B/C -> D only; upgraded all-8 calibrated LOPO exists | Partial | Need paper's exact within-patient and cross-patient evaluations for all methods. |
| HER2 missing imputation | Fig. 2i, patient A slice 2 posterior probabilities | Local toy missing-imputation demo only | Not verified | Need A2 disease labels/missing locations and paper setup. |
| Alzheimer/STARmap PLUS | Cst7/Trem2/C1qa/Gfap; Hexb*Trem2; Replicate 2 accuracy 97.1%; Fig. B12 | No local data or result log. `SINGLECELL_DATA.md` notes SCP1375/sign-in required. Zenodo 5842625 files are restricted; SCP1375 API requires auth. | Blocked | Need authenticated STARmap PLUS Replicate 1/2 expression data and plaque coordinates/labels. |

## Already verified local evidence

`repro/_verify/progress.txt` shows one complete local verification pass:

- `verify_her2`: rc=0
- `singlecell`: rc=0
- `nuts_compare`: rc=0
- `reproduce`: rc=0
- `her2_reproduce`: rc=0
- `benchmark`: rc=0
- `plot_scaling`: rc=0
- `sim2_scaling`: rc=0

`repro/_verify/REPORT.md` summarizes deterministic checks and timing checks. These are valid
for the local verification kit, but should not be described as full btaf530 reproduction.

## Required next actions for full verification

1. Parse the obtained supplementary material into machine-checkable targets.
   Required pieces now available locally: Section A2, Section A3, Fig. B2-B12, Tables C1-C7.

2. Build an original-paper simulation runner.
   Use the supplement to reproduce Simulation 1/2/3 exactly, then compare bias, SD,
   coverage, runtime, and sensitivity outputs against Tables C1-C6 and Fig. B2-B10.
   The current eta-bias diagnosis and public script search indicate that further progress
   likely requires the authors' exact simulation code or an additional unpublished
   convention, because the public paper formula and local package implementation imply
   different eta scales and no Simulation 1/2/3 scripts were found in the declared
   GitHub/Zenodo resources.

3. Resolve HER2 labels.
   The `her2st` repository has counts for 36 sections, but pathologist label TSVs only for
   A1/B1/C1/D1/E1/F1/G2/H1/J1. Full paper verification needs the disease labels for every
   section actually used by btaf530.

4. Run the original HER2 ADVI job only after labels/protocol are fixed.
   Expected runtime is paper-scale, about 25h32m. Compare top-30 genes, coefficient
   directions/intervals, eta, interaction network, Giotto/SVG overlaps, prediction accuracies,
   and missing imputation.

5. Obtain STARmap PLUS Alzheimer data.
   The local notes identify Broad SCP1375 / restricted Zenodo as likely source. Full
   verification needs Replicate 1 discovery and Replicate 2 prediction under the paper's
   plaque-proximity label rule.

## Current claim wording

Accurate wording:

> We have verified the local DiSTect demonstration and extension workflows, including a
> scaled HER2 reproduction and upgraded all-8-patient PG-CAVI analysis. We have not yet fully
> reproduced all btaf530 paper experiments. The official supplementary material has now been
> obtained via Europe PMC, but full HER2 labeling/protocol resolution, original full ADVI runtime,
> baseline comparisons, and STARmap PLUS data reproduction are still incomplete.
