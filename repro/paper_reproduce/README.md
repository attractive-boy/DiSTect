# Paper reproduction

This directory contains a self-contained reconstruction of the simulations in
Zhao, Deng, and Zhang (2025), plus diagnostics for the reported uncertainty of
the spatial parameter `eta`.

## Requirements

- R 4.x
- R packages `Matrix` and `rstan`
- A working C++ toolchain supported by `rstan`

Run every command from the repository root. Scripts locate the repository from
their own file path, so the checkout can live anywhere.

## Main simulation runs

Start with a smoke test:

```bash
N_REP=1 Rscript repro/paper_reproduce/run_sim1.R
N_REP=1 Rscript repro/paper_reproduce/run_sim2.R
N_REP=1 Rscript repro/paper_reproduce/run_sim3.R
```

The paper-scale ADVI comparisons use 200 deterministic replicates:

```bash
N_REP=200 Rscript repro/paper_reproduce/run_sim1.R
N_REP=200 Rscript repro/paper_reproduce/run_sim2.R
N_REP=200 Rscript repro/paper_reproduce/run_sim3.R
Rscript repro/paper_reproduce/compare_sim1_advi.R
Rscript repro/paper_reproduce/compare_sim2_advi.R
Rscript repro/paper_reproduce/compare_sim3_advi.R
```

Outputs are written to `repro/paper_reproduce/output/`. The committed
`*_aggregate_n200.csv` and `*_key_comparison_n200.csv` files are reference
results; raw replicate-level output and logs are intentionally not tracked.

## Eta uncertainty checks

The independent checks used in `ETA_UNCERTAINTY_FINDING.md` are:

```bash
Rscript repro/paper_reproduce/confirm_meanfield.R
Rscript repro/paper_reproduce/confirm_covscale.R
Rscript repro/paper_reproduce/crlb_vs_covscale.R
Rscript repro/paper_reproduce/nuts_sim1_tuned.R
Rscript repro/paper_reproduce/mcmle_autologistic.R
Rscript repro/paper_reproduce/diagnose_uncertainty_scaling.R
Rscript repro/paper_reproduce/reverse_search_uncertainty_formulas.R
N_REP=10 ADVI_ITER=3000 Rscript repro/paper_reproduce/run_targeted_advi_grid.R
Rscript repro/paper_reproduce/make_figures.R
```

These checks support reproduction auditing; they do not claim an exact numeric
match to every table in the paper. See `STATUS.md` for the current boundary and
`ETA_UNCERTAINTY_FINDING.md` for the evidence chain.
