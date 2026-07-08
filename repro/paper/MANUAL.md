# 说明书 — Analysis kit for the scalable/count-aware DiSTect paper

A self-contained research kit on branch `feat/sparse-scalable-m0`. It reproduces the original
DiSTect, implements the upgraded method, and runs every simulation and real-data experiment behind
the paper story (`paper/STORY.md`). Everything is R + rstan; no new package install needed beyond
what the original DiSTect requires (`rstan`, `Matrix`, `ggplot2`; `rjags` only for the original
missing-imputation).

## 0. TL;DR — run everything
```bash
cd /Users/licheng/Documents/DiSTect
bash repro/realdata/fetch_her2st.sh      # ~30s: download labeled HER2 sections
bash repro/run_all.sh                     # runs sims + real-data pipelines, prints results
```
Individual pieces below. Fast scripts (PG-CAVI) finish in seconds; the Stan-based ones compile once
(~30s) then fit.

## 1. Directory map
```
repro/
  m0/          De-risk milestone: sparse O(N) == dense, + scaling benchmark
    neighbors.R  dsgd_sparse.R  verify_her2.R  bench_one.R  run_benchmark.sh
    plot_scaling.R  scaling.csv  fig_scaling_*.png  README.md
  method/      The upgraded method (reusable library)
    neighbors.R          grid-hash neighbor sum + sparse adjacency
    fit_sparse.R         sparse ADVI/NUTS fitters (== model, O(N))
    fit_polyagamma.R     ★ PG-CAVI closed-form inference (the fast/stable engine)
    likelihood_nb.R      ★ NB size-factor normalization (count layer)
    eta_fulllik.R        full-likelihood η de-biasing
    selection_fdr.R      FDR-controlled gene selection (BH + Bayesian)
    predict_calibrated.R posterior-predictive prediction + Brier/ECE
  sim/         Simulation studies (paper Fig 2-5)
    sim1_correctness.R  sim2_scalability.R  sim3_count_regime.R  sim4_pseudolik_bias.R
  realdata/    Real-data pipelines (paper Fig 6-7)
    fetch_her2st.sh  run_her2st.R  run_singlecell.R  SINGLECELL_DATA.md
  paper/       STORY.md (narrative + evidence) + this MANUAL.md
  run_all.sh
  # plus the original-paper reproduction: reproduce.R, her2_reproduce.R, nuts_compare.R
```
★ = the novel pieces.

## 2. The method API (how to call it on your own data)
```r
source("repro/method/neighbors.R")
source("repro/method/fit_polyagamma.R")
source("repro/method/likelihood_nb.R")
source("repro/method/selection_fdr.R")

X   <- scale(nb_pearson(counts))            # counts: spots x genes (raw)  -> NB covariates
fit <- fit_pgcavi_single(y, X, coords = coords, label = patient)  # y: 0/1; coords: x,y
fit$table        # per-gene mean, sd, std_effect (intercept + genes + eta)
fit$eta          # spatial autocorrelation
select_fdr_z(fit$table$std_effect[genes], level = 0.10)   # FDR-controlled selection
```
Key inputs: `coords` must be an integer lattice (Visium/ST already are; snap single-cell coords —
see `run_singlecell.R`). `label` gives per-slice/patient membership (neighbors are label-matched).

## 3. Reproduce, piece by piece
| goal | command | ~time |
|---|---|---|
| Original paper (toy) | `Rscript repro/reproduce.R` | 3 min |
| Original HER2 (4 patients) | `Rscript repro/her2_reproduce.R` | 10 min |
| NUTS vs VI | `Rscript repro/nuts_compare.R` | 3 min |
| **M0** equivalence | `Rscript repro/m0/verify_her2.R` | 3 min |
| **M0** scaling benchmark | `bash repro/m0/run_benchmark.sh && Rscript repro/m0/plot_scaling.R` | ~40 min |
| **Sim1** correctness | `Rscript repro/sim/sim1_correctness.R` | 5 s |
| **Sim2** scaling | `Rscript repro/sim/sim2_scalability.R` | 3 min |
| **Sim3** count regime | `Rscript repro/sim/sim3_count_regime.R` | 20 s |
| **Sim4** η bias | `Rscript repro/sim/sim4_pseudolik_bias.R` | 30 s |
| **HER2** upgraded (all 8) | `Rscript repro/realdata/run_her2st.R` | 15 s |
| **Single-cell** (synthetic) | `Rscript repro/realdata/run_singlecell.R` | 30 s |
| **Single-cell** (real) | `DATA_RDS=x.rds Rscript repro/realdata/run_singlecell.R` | varies |

## 4. Headline numbers already produced (see STORY.md §3)
- Sparse ≡ dense (NUTS): cor(β)=0.9996. Scaling: 1464× at N=3000; PG-CAVI 300k cells in 10 s.
- Count layer: mean-var coupling 0.96→0.4; depth confound 0.95→0.02.
- HER2 all-8 in 3.79 s (vs ~25 h); 26 genes @q<0.10; LOPO acc 0.740 / Brier 0.180.

## 5. Extending toward submission (ordered)
1. Real single-cell dataset → `SINGLECELL_DATA.md`, then `run_singlecell.R` with `DATA_RDS`.
2. Spike-and-slab on PG-CAVI: hook documented at the bottom of `method/fit_polyagamma.R`.
3. Latent-NB measurement model: hook documented in `method/likelihood_nb.R`.
4. Multi-slice random effects in PG-CAVI (batch); currently pooled via `label`.
5. Ablations, external validation, methods/theory write-up. Figure→script map in STORY.md §4.

## 6. Gotchas
- `coords` must be integer lattice for `neighbor_sum`/`build_adjacency`; snap continuous coords.
- Mean-field VI under-estimates variance → selection can be liberal; prefer the spike-and-slab /
  Bayesian-FDR route for final selection.
- `run_her2st.R` needs `bash repro/realdata/fetch_her2st.sh` first (data in `/tmp/her2st`).
- The original `R/dsgd.R` was left unchanged except an ADVI iteration cap (committed to `main`).
