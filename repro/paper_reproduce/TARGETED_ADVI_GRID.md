# Targeted 10-replicate ADVI grid

Date: 2026-07-11

## Scope

This experiment tests the five remaining implementation hypotheses without a
new 200-replicate run:

- meanfield versus fullrank ADVI;
- posterior SD versus Monte Carlo error of the variational draws;
- neighbor average versus neighbor sum;
- fixed versus regenerated covariate matrices;
- the four multiple-slice covariance settings in Table C4.

Each configuration uses 10 deterministic replicates and 3000 ADVI iterations.
The experiment performs 80 single-slice and 80 multiple-slice fits.

## Single-slice result

Paper Table C1 Proposed ADVI eta targets are: bias `0.087`, SEE `0.044`, SEM
`0.044`, coverage `0.95`.

| algorithm | neighbor | X | bias | empirical SD | posterior SD | draw MCSE | coverage |
|---|---|---|---:|---:|---:|---:|---:|
| meanfield | average | fixed | 0.147 | 0.187 | 0.259 | 0.0082 | 1.0 |
| meanfield | average | regenerated | 0.136 | 0.163 | 0.258 | 0.0082 | 1.0 |
| fullrank | average | fixed | 0.601 | 1.444 | 0.448 | 0.0142 | 0.8 |
| fullrank | average | regenerated | 0.411 | 0.771 | 0.388 | 0.0123 | 0.9 |
| meanfield | sum | fixed | 0.106 | 0.115 | 0.072 | 0.0023 | 0.8 |
| meanfield | sum | regenerated | 0.097 | 0.111 | 0.071 | 0.0023 | 0.7 |
| fullrank | sum | fixed | 0.473 | 1.305 | 0.288 | 0.0091 | 0.9 |
| fullrank | sum | regenerated | **0.048** | **0.053** | **0.087** | 0.0028 | 1.0 |

### Interpretation

1. The neighbor-sum convention is the only tested change that moves eta bias
   and empirical SD close to Table C1. This is important because the package
   implementation uses a sum while the paper formula states an average.
2. Fullrank plus neighbor sum plus regenerated X is the closest configuration:
   empirical SD `0.053` versus paper `0.044`. Its posterior SD is still twice
   the paper value and coverage is 100% in this small run.
3. Fixed versus regenerated X has little effect under meanfield. Under fullrank,
   regenerated X avoids several extreme failures, so the apparent effect is
   mainly stability rather than a different uncertainty formula.
4. Fullrank is not robust. Several configurations contain catastrophic eta
   outliers and very large Pareto-k warnings; meanfield is much more stable.

## Posterior SD versus Monte Carlo error

`rstan::vb` does not expose a Stan `se_mean` column. The calculable Monte Carlo
error of the 1000 variational draws is `posterior SD / sqrt(1000)` and ranges
from `0.0023` to `0.0142` in the single-slice grid. These values are generally
below the paper's `0.044`, so confusing posterior SD with variational-draw MCSE
does not explain Table C1.

## Multiple-slice covariance result

Paper Table C4 eta targets across the four settings are bias `0.076-0.089`, SEE
`0.022-0.031`, and SEM `0.022-0.031`.

| algorithm | sigma2 | rho | bias | empirical SD | posterior SD | draw MCSE | coverage |
|---|---:|---:|---:|---:|---:|---:|---:|
| meanfield | 0.1 | 0.1 | 0.171 | 0.191 | 0.108 | 0.0034 | 0.7 |
| meanfield | 0.1 | 0.4 | 0.268 | 0.312 | 0.107 | 0.0034 | 0.4 |
| meanfield | 0.4 | 0.1 | 0.294 | 0.434 | 0.101 | 0.0032 | 0.4 |
| meanfield | 0.4 | 0.4 | 0.255 | 0.321 | 0.103 | 0.0032 | 0.5 |
| fullrank | 0.1 | 0.1 | 0.605 | 1.183 | 0.403 | 0.0128 | 0.8 |
| fullrank | 0.1 | 0.4 | 0.253 | 0.299 | 0.270 | 0.0085 | 1.0 |
| fullrank | 0.4 | 0.1 | 0.153 | 0.195 | 0.224 | 0.0071 | 1.0 |
| fullrank | 0.4 | 0.4 | 0.276 | 0.322 | 0.235 | 0.0074 | 0.7 |

No covariance setting approaches the Table C4 uncertainty. Meanfield posterior
SD remains about four times the paper values; empirical SD is roughly 6-20
times larger. Fullrank is wider and less stable. The duplicated `{0.1, 0.4}`
entry in the supplement text therefore cannot explain the published table.

## Decision

- Do not run a broad 200-replicate ADVI grid.
- The only configuration worth a larger confirmation is single-slice,
  neighbor-sum, regenerated X, comparing meanfield and fullrank with stricter
  convergence diagnostics.
- Table C4 requires investigation of the multiple-slice model/data generator,
  not further tuning of `sigma2`, `rho`, or the displayed uncertainty formula.

## Reproduce

```bash
N_REP=10 ADVI_ITER=3000 Rscript repro/paper_reproduce/run_targeted_advi_grid.R
```

Committed aggregate outputs:

- `output/targeted_single_advi_aggregate_n10.csv`
- `output/targeted_multi_advi_aggregate_n10.csv`

Raw per-fit outputs are intentionally not tracked.
