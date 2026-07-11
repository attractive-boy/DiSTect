# Reverse search of avgSEE / avgSEM formulas

Date: 2026-07-11

## Question

Can one common aggregation formula reproduce the published uncertainty columns
in Tables C1-C6 from the existing 200-replicate outputs?

## Search design

`reverse_search_uncertainty_formulas.R` evaluates 12 candidate definitions over
81 published cells for each metric, covering Simulation 1-3:

- empirical SD and empirical SD divided by `sqrt(200)`;
- RMSE and RMSE divided by `sqrt(200)`;
- mean absolute error and mean absolute error divided by `sqrt(200)`;
- mean or RMS posterior SD, with and without division by `sqrt(200)`;
- mean posterior variance and `sqrt(mean(variance) / 200)`.

Candidates are scored on the log ratio to the paper target. Log error treats a
twofold overestimate and a twofold underestimate symmetrically.

## Result

No single formula explains all tables.

| scope | metric | best candidate | median candidate / paper | cells within 2x |
|---|---|---|---:|---:|
| all 81 cells | SEE | `RMSE / sqrt(200)` | 0.66 | 44% |
| Tables C1-C3 | SEE | `SD / sqrt(200)` by log RMSE; `RMSE / sqrt(200)` by median | 0.57 / 0.98 | 52% / 45% |
| Table C4 | SEE | unscaled empirical SD | 2.78 | 17% |
| Tables C5-C6 | SEE | `RMSE / sqrt(200)` | 0.75 | 71% |
| all 81 cells | SEM | posterior SD divided by `sqrt(200)` | 0.50 | 47% |
| Tables C1-C3 | SEM | posterior SD divided by `sqrt(200)` | 0.61 | 73% |
| Table C4 | SEM | no satisfactory candidate | - | <=21% |
| Tables C5-C6 | SEM | posterior SD divided by `sqrt(200)` | 0.52 | 46% |

For the 11 eta cells, `RMSE / sqrt(200)` has median ratio `0.985`, but only 36%
of cells fall within 25% of the paper target. This means the near-one median is
meaningful evidence of the reporting scale, but not proof of a universal exact
formula.

## Interpretation

1. `sqrt(200)` is not a coincidence. Every best global SEE/SEM candidate except
   Table C4 contains division by `sqrt(200)`, and the independent NUTS eta check
   remains a strong direct match.
2. It is not the complete explanation. A single scaled formula fits fewer than
   half of all cells within a factor of two.
3. Table C4 is qualitatively different. Neither scaled nor unscaled candidates
   fit it, indicating that the multiple-slice data generation or fitted model is
   mismatched before aggregation is considered.
4. The most plausible combined explanation is an extra `sqrt(200)` reporting
   scale plus separate implementation differences in ADVI, especially for the
   multiple-slice and missing-data models.

## Decision on further ADVI runs

A broad, expensive ADVI grid is not justified yet: changing ADVI settings cannot
repair Table C4's data-generation mismatch or identify the reporting formula.
The next computational experiment should be a small targeted grid (10-20
replicates) restricted to:

- meanfield versus fullrank summaries;
- posterior SD versus Stan-reported Monte Carlo error;
- neighbor average versus neighbor sum;
- fixed versus regenerated covariate matrices;
- the exact multiple-slice covariance settings, including the duplicated
  `{0.1, 0.4}` entry in the supplement text.

Only configurations that improve point-estimate and uncertainty agreement
simultaneously should be promoted to a 200-replicate run.

## Outputs

- `output/uncertainty_formula_search_scores.csv`: formula scores by scope.
- `output/uncertainty_formula_search_best.csv`: best candidate per published
  cell.

Regenerate without model fitting:

```bash
Rscript repro/paper_reproduce/reverse_search_uncertainty_formulas.R
```
