# btaf530 stepwise reproduction status

Date: 2026-07-11

## UPDATE (2026-07-11): likely sqrt(200) reporting-scale error identified

The strongest reconstruction of Tables C1-C3 is that `avgSEE/avgSEM` were
divided by `sqrt(200)`, turning estimator-level uncertainty into Monte Carlo
uncertainty for the 200-replicate average. The Table C1 NUTS eta values
`0.021/0.022` become `0.297/0.311`, matching the independent local NUTS/GLM
range `0.24-0.32` and restoring consistency with 95% coverage. This substantially
resolves the first reproduction question as a likely reporting-scale issue.
The exact source formula and the remaining ADVI numerical mismatch cannot be
confirmed without the authors' simulation script. See
`ETA_UNCERTAINTY_FINDING.md` and `output/uncertainty_scaling_audit.csv`.

Formula reverse-search over 81 cells in Tables C1-C6 further shows that the
`sqrt(200)` explanation is systematic but incomplete. Scaled RMSE/posterior SD
is the best family for C1-C3 and C5-C6, while no candidate explains C4. This
separates the problem into a likely reporting-scale issue and an independent
multiple-slice/ADVI implementation mismatch. See
`UNCERTAINTY_FORMULA_SEARCH.md`.

## UPDATE (2026-07-11): targeted 10-replicate ADVI grid completed

The 160-fit targeted grid tested meanfield/fullrank, posterior SD/draw MCSE,
neighbor average/sum, fixed/regenerated X, and all four Table C4 covariance
settings. Neighbor sum plus regenerated X is the only single-slice configuration
that approaches Table C1 (fullrank eta bias `0.048`, empirical SD `0.053` versus
paper `0.087/0.044`), though posterior SD remains `0.087` and fullrank has severe
Pareto-k instability. No multiple-slice covariance setting approaches C4.
Details: `TARGETED_ADVI_GRID.md`.

## UPDATE (2026-07-09): eta discrepancy resolved to an uncertainty-column issue

The eta discrepancy is NOT a local bug, ADVI setting, prior mismatch, covariate
scale, or pseudolikelihood-vs-full-likelihood issue. The eta **point estimates**
reproduce (roughly unbiased; proposed << naive signal-gene bias). What does not
reproduce is the eta **uncertainty** columns (avgSEE/avgSEM) in Tables C1-C3:
they are 2.5x-5x below the ~0.11 Cramer-Rao floor for the stated 30x30 design,
established four independent ways (GLM, ADVI meanfield+fullrank, divergence-free
NUTS, and full-likelihood Geyer-Thompson MCMLE), and the floor is scale-invariant
(does not drop below ~0.11 even with covariates removed). The reported ADVI row
is also internally inconsistent (bias ~2x its SE yet 95% coverage). Most likely a
reporting/definitional issue in the SE columns, not an error in the biology.
Full evidence chain and caveats: `ETA_UNCERTAINTY_FINDING.md`. Author query:
`AUTHOR_SIM_SCRIPT_REQUEST_EMAIL.md`.

Source paper: [Zhao, Deng, and Zhang (2025), Bioinformatics, btaf530](https://doi.org/10.1093/bioinformatics/btaf530)

Supplement: supplementary data distributed with the journal article.

## Status table

| Block | Paper target | Current status | Local outputs |
|---|---|---|---|
| A2.1 Simulation 1 ADVI | Tables C1-C3 ADVI columns; Fig. B2-B6 related designs | Completed first full `N_REP=200` ADVI pass. Directionally reproduces proposed-vs-naive improvement, but does not numerically match Tables C1-C3. Main gaps: naive coverage is much lower than target for signal genes, and proposed eta bias is much higher than target. | `output_sim1_n200_advi.log`; `output/sim1_raw_n200.csv`; `output/sim1_aggregate_n200.csv`; `output/sim1_timings_n200.csv`; `output/sim1_advi_key_comparison_n200.csv`. |
| A2.1 Simulation 1 NUTS | Tables C1-C3 NUTS columns; Fig. B2 runtime | Not started. Paper reports NUTS is approximately 23h for this scenario, so this needs explicit long-run scheduling. | None yet. |
| A2.2 Simulation 2 | Fig. B7-B8, Table C4 | Completed first full `N_REP=200` ADVI pass for all four `{sigma2, rho}` settings. Does not numerically match Table C4: beta biases are much lower than target, while eta bias is much higher than target. | `output_sim2_n200_advi.log`; `output/sim2_raw_n200.csv`; `output/sim2_aggregate_n200.csv`; `output/sim2_timings_n200.csv`; `output/sim2_advi_key_comparison_n200.csv`. |
| A2.3 Simulation 3 | Fig. B9-B10, Tables C5-C6 | Completed first full `N_REP=200` ADVI pass with a missing-neighbor first-pass approximation. Does not numerically match Tables C5-C6: eta bias remains high and coverage is below target. Nonignorable missing counts are lower than the paper's stated averages under the direct `M_i=1` missing interpretation of Equation (5). | `output_sim3_n200_advi.log`; `output/sim3_raw_n200.csv`; `output/sim3_aggregate_n200.csv`; `output/sim3_timings_n200.csv`; `output/sim3_missing_counts_n200.csv`; `output/sim3_advi_key_comparison_n200.csv`. |
| HER2 original ADVI | Main Fig. 2, Table C7 hyperparameters | Blocked for full-paper reproduction. Public HER2ST tree has counts/spotfiles for 36 sections, but pathologist labels for A-H exist only for A1/B1/C1/D1/E1/F1/G2/H1. Existing scaled 4-section original-ADVI run and all-8 labeled-section upgraded run do not match paper Fig. 2 targets. | `output/her2_section_availability.csv`; existing logs `repro/_verify/her2_reproduce.log`, `repro/realdata/her2_new.log`; existing figures `repro/fig_her2_coef.png`, `repro/fig_her2_network.png`. |
| STARmap PLUS Alzheimer | Supplement A3, Fig. B12 | Blocked by data access. Local repo/Zenodo DiSTect archive does not contain STARmap PLUS data. Zenodo record `5842625` is public metadata but restricted files; Zenodo API returns 403. Broad Single Cell Portal `SCP1375` API returns 401 without authentication. GitHub searches found no public mirror. | Source checks: `repro/realdata/SINGLECELL_DATA.md`; attempted endpoints `https://zenodo.org/records/5842625`, `https://zenodo.org/api/records/5842625`, `https://singlecell.broadinstitute.org/single_cell/api/v1/studies/SCP1375`. |
| Eta-bias differential diagnosis | Explain systematic eta overestimation in simulations | Completed first diagnostic pass. High eta bias persists in plain GLM conditional pseudolikelihood, random-scan single-site Gibbs, and ADVI meanfield/fullrank pilots. Intercept and centered-neighbor variants do not fix it. RESOLVED (see UPDATE above and `ETA_UNCERTAINTY_FINDING.md`): the reproducible discrepancy is the eta uncertainty columns (avgSEE/avgSEM ~2.5-5x below the ~0.11 Cramer-Rao floor, scale-invariant, confirmed by GLM/ADVI/NUTS/MCMLE), plus an internal SEM-vs-coverage inconsistency; eta point estimates reproduce. | `ETA_UNCERTAINTY_FINDING.md`; `ETA_BIAS_DIAGNOSIS.md`; `output/eta_bias_glm_aggregate_n200.csv`; `output/eta_generation_grid_aggregate_n50.csv`; `output/eta_single_site_glm_raw_n20.csv`; `output/eta_advi_algorithm_raw_n5.csv`; `output/eta_glm_intercept_raw_n200.csv`; `output/eta_glm_centered_raw_n200.csv`. |
| Original simulation script search | Find authors' Simulation 1/2/3 source scripts or historical code | Public search completed across paper-declared GitHub `StaGill/DiSTect`, Zenodo `17127211`, official fork `AnjiDeng/DiSTect`, repository history, tutorial page, issues/PRs/releases/wiki, author/org public repos, and Zenodo keyword search. No original paper simulation scripts were found. | `ORIGINAL_SIM_SCRIPT_SEARCH.md`. |
| Author clarification request | Ask for original scripts / eta convention | Drafted a targeted email to the corresponding author with the exact mismatch, diagnostics already ruled out, and five implementation questions. | `AUTHOR_SIM_SCRIPT_REQUEST_EMAIL.md`. |

## Implementation notes

- The paper's model uses `eta / |N(i)| * sum_j Y_j`, i.e. neighbor average. This differs
  from the current local package implementation in `R/dsgd.R`, which uses an unnormalized
  neighbor sum. Paper reproduction scripts therefore implement the normalized model directly.
- Simulation 1 uses a 30 x 30 lattice, 20 standard normal covariates, beta
  `(1, 2, 3, -4, -5, zero_15)`, eta `{0.4, 1.6, 2.8}`, 2000 Gibbs sweeps, and 200
  simulation replicates for the full ADVI reproduction.
- The default simulator uses checkerboard block Gibbs on the rook-neighbor lattice. This is a
  valid block Gibbs update for this bipartite graph and is much faster than pure R single-site
  random-scan Gibbs.
- Pilot runs produced outputs for `N_REP=1`, `N_REP=3`, and `N_REP=10`. These are smoke tests,
  not final paper verification.

## Current findings

- Simulation 1 ADVI is not fully reproduced yet. The implementation follows the paper-normalized
  formula and 200-replicate design, but key numeric comparisons against Tables C1-C3 fail.
- The strongest qualitative match is that the proposed ADVI model has much smaller signal-gene
  bias than the naive ADVI model.
- The strongest numeric mismatch is eta: target avgBias is `0.087`, `0.083`, `0.095` for eta
  `0.4`, `1.6`, `2.8`, while the current run gives `0.231`, `0.366`, `0.437`.
- ADVI emits repeated Pareto-k warnings, especially for the naive model. These warnings are a
  reproducibility risk and may partly explain unstable coverage.
- Simulation 2 ADVI also does not fully reproduce Table C4. For example, target eta avgBias is
  `0.076`, `0.089`, `0.083`, `0.081` across the four settings, while the current run gives
  `0.373`, `0.374`, `0.333`, `0.375`. Signal beta biases are consistently below the table
  targets, suggesting either a remaining data-generation/modeling mismatch or table values
  generated under a different implementation detail.
- Simulation 3 ADVI also does not fully reproduce Tables C5-C6. Current eta avgBias is
  `0.325`, `0.353`, `0.325`, `0.309` for the four missing settings, compared with targets
  `0.089`, `0.156`, `0.062`, `0.066`. The nonignorable missing generator gives mean missing
  counts `4.17` and `11.90`, while the supplement states approximately `9.9` and `33.3`
  missing spots out of 900. This indicates an unresolved convention or implementation detail
  in Equation (5)'s missing indicator.
- Eta-bias diagnostics now indicate that the systematic eta gap is not caused by ADVI alone.
  A plain GLM conditional pseudolikelihood on the same Simulation 1 datasets gives eta
  average absolute biases `0.250`, `0.291`, and `0.378` for eta `0.4`, `1.6`, and `2.8`.
  Random-scan single-site Gibbs also gives high eta biases `0.270`, `0.236`, and `0.355`.
  The paper formula's neighbor-average convention differs from the local package's
  unnormalized neighbor-sum implementation, but generation-grid diagnostics show that the
  sum convention cannot recover all three paper eta targets and yields implausibly high
  disease rates for larger eta. See `ETA_BIAS_DIAGNOSIS.md`.
- The authors' exact Simulation 1/2/3 scripts were not found in the public code/data
  resources advertised by the paper. The official GitHub, Zenodo record, official fork,
  visible repository history, and tutorial contain package source/toy data only. See
  `ORIGINAL_SIM_SCRIPT_SEARCH.md`.
- Additional reverse-engineering checks also failed to explain the eta gap. Fixed covariates
  across replicates, one-pass parallel generation, `plogis(X beta)` initialization, and
  parallel-vs-checkerboard 2000-sweep generation all keep eta average absolute bias well
  above the paper targets. See `output/eta_generation_extra_grid_raw_n50.csv`.
- HER2 full original reproduction is currently blocked by public label availability. The public
  HER2ST repository exposes 36 count matrices and spotfiles, but only 8 A-H pathologist-label
  files. The paper describes patients A-D with six slices and E-H with three slices, so the
  exact disease-status labels used for the remaining sections are missing from the local/public
  material inspected so far.
- HER2 scaled original-ADVI run on A1/B1/C1/D1 gives eta `2.665 [1.589, 3.999]` and top genes
  `TPT1`, `MMP14`, `SCGB2A2`, etc.; paper target is eta `2.936 [2.926, 2.946]` and top genes
  including `MUCL1`, `NDRG1`, `VEGFA`, `CXCL10`, `SCD`, `ERBB2`. The all-8 local upgraded
  PG-CAVI run is useful but not paper-original.
- STARmap PLUS Alzheimer reproduction is blocked until authenticated data are available. Paper
  targets are: Replicate 1 top genes `Cst7` beta `1.580 [1.551, 1.609]`, `Trem2` beta
  `0.992 [0.957, 1.027]`, `C1qa` beta `0.726 [0.693, 0.759]`, `Gfap` beta
  `0.577 [0.527, 0.627]`; interaction `Hexb` beta `0.022 [0.011, 0.033]`; Replicate 2
  prediction accuracy `97.1%` with disease labels defined as cells within 64 pixels
  (approximately 20 um) of an amyloid-beta plaque center and a 60-pixel neighborhood radius.
