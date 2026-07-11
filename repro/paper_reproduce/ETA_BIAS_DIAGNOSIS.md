# Eta bias diagnosis

Date: 2026-07-09

Question: why is eta systematically higher than the paper targets in the local
simulation reproductions?

## Short conclusion

The elevated eta bias is not explained by ADVI alone. The same direction and
similar magnitude appear when Stan/ADVI is bypassed and eta is estimated by a
plain conditional logistic GLM on the same generated lattices. Changing the
Gibbs update order, the number of burn-in sweeps, ADVI meanfield/fullrank, adding
an intercept, or centering the neighbor term does not recover the paper's eta
avgBias values.

The current evidence points to an unresolved mismatch between the paper's
reported simulation tables and the simulation/fitting details available from the
paper text plus this repository. The largest concrete inconsistency found so far
is that the paper formula uses the neighbor average,

`eta / |N(i)| * sum_{j in N(i)} Y_j`,

while the package implementation in `R/dsgd.R` uses the unnormalized neighbor
sum. However, switching to the unnormalized sum does not reproduce all eta
settings and gives implausibly high disease rates for larger eta.

## Evidence

### 1. Full ADVI reproduction has high positive eta bias

From `output/sim1_aggregate_n200.csv`, proposed ADVI eta estimates under the
paper-normalized model are:

| eta truth | mean eta_hat | signed bias | avg abs bias | paper ADVI avgBias |
|---:|---:|---:|---:|---:|
| 0.4 | 0.588 | +0.188 | 0.231 | 0.087 |
| 1.6 | 1.907 | +0.307 | 0.366 | 0.083 |
| 2.8 | 3.197 | +0.397 | 0.437 | 0.095 |

The gap is therefore not just an absolute-vs-signed-bias reporting issue.

### 2. GLM conditional pseudolikelihood reproduces the problem

`diagnose_eta_bias.R` fits

`Y ~ X1 + ... + X20 + neighbor_avg - 1`

using base R `glm`, with no Stan and no ADVI. Output:
`output/eta_bias_glm_aggregate_n200.csv`.

| eta truth | GLM mean eta_hat | GLM avg abs bias |
|---:|---:|---:|
| 0.4 | 0.438 | 0.250 |
| 1.6 | 1.773 | 0.291 |
| 2.8 | 3.052 | 0.378 |

Because the high eta error persists outside ADVI, ADVI instability is not the
root cause.

### 3. Gibbs update convention does not fix it

`diagnose_eta_generation_grid.R` varied:

- update mode: parallel vs checkerboard
- sweeps: 20, 100, 500, 2000
- spatial scale: neighbor average vs unnormalized sum

No combination matched the three paper eta targets simultaneously.

The unnormalized neighbor-sum convention can get eta=0.4 close to the target
bias, but for larger eta it produces very high disease rates:

| eta truth | update | sweeps | scale | avg abs bias | disease rate |
|---:|---|---:|---|---:|---:|
| 1.6 | parallel | 500 | sum | 0.140 | 0.711 |
| 2.8 | parallel | 100 | sum | 0.370 | 0.873 |

This makes the repository's unnormalized implementation unlikely to be the full
explanation for the paper tables.

### 4. Single-site Gibbs does not fix it

A slower random-scan single-site Gibbs diagnostic produced:
`output/eta_single_site_glm_raw_n20.csv`.

| eta truth | mean eta_hat | avg abs bias | disease rate |
|---:|---:|---:|---:|
| 0.4 | 0.458 | 0.270 | 0.515 |
| 1.6 | 1.758 | 0.236 | 0.538 |
| 2.8 | 3.085 | 0.355 | 0.587 |

So the checkerboard block update is not the source of the bias.

### 5. ADVI meanfield/fullrank is not enough

A 5-replicate pilot saved to `output/eta_advi_algorithm_raw_n5.csv` compared
meanfield and fullrank ADVI. It did not move eta toward the paper values; both
variants also produced repeated Pareto-k warnings.

### 6. Intercept and centered-neighbor variants were ruled out

Adding an intercept to the GLM worsened eta absolute bias:
`output/eta_glm_intercept_raw_n200.csv`.

Centering `neighbor_avg` also worsened eta absolute bias:
`output/eta_glm_centered_raw_n200.csv`.

### 7. Extra generation conventions were also tested

`output/eta_generation_extra_grid_raw_n50.csv` tested several plausible hidden
implementation choices:

- fixed X reused across replicates
- one parallel update instead of 2000 sweeps
- initialization from `plogis(X beta)` instead of Bernoulli(0.5)
- parallel 2000 updates instead of checkerboard 2000 updates

None recovered the paper eta targets. Representative average absolute biases:

| mode | fixed X | eta=0.4 | eta=1.6 | eta=2.8 |
|---|---:|---:|---:|---:|
| checkerboard 2000 | false | 0.252 | 0.295 | 0.437 |
| one parallel update | false | 0.195 | 0.232 | 0.394 |
| parallel 2000 | false | 0.259 | 0.254 | 0.318 |
| checkerboard 2000 | true | 0.278 | 0.266 | 0.329 |
| parallel 2000 | true | 0.290 | 0.279 | 0.308 |

These values remain well above the paper ADVI eta avgBias targets
`0.087`, `0.083`, and `0.095`.

## Current interpretation

Under the paper-normalized formula and the stated 30 x 30 / 20 covariate / 2000
Gibbs-sweep design, eta is much less tightly recovered than the paper's tables
report. The paper tables' eta standard errors are also much smaller than the
local GLM and ADVI uncertainty scales, suggesting that an unpublished
implementation detail, code path, or reporting convention differs from the
written supplement.

The public search for the authors' exact simulation scripts has now been
completed across the paper-declared GitHub/Zenodo resources, visible repository
history, official fork, tutorial page, and author/org public repositories. No
Simulation 1/2/3 source scripts were found; see `ORIGINAL_SIM_SCRIPT_SEARCH.md`.
The next most useful action is therefore to contact the authors for the exact
scripts or for clarification of the simulation/fitting/reporting convention used
to produce Tables C1-C6.
