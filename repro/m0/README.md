# M0 — Sparse, O(N) Scalable DiSTect (de-risk prototype)

Milestone **M0** of the scalability direction: prove the published autologistic fit can be made
**O(N)** without changing the model, on real (her2st) and synthetic data — the foundation for the
JAX/GPU port (M1) and the count-likelihood extension (M2).

Branch: `feat/sparse-scalable-m0`. Everything here is self-contained; run from the repo root.

## The insight

`R/dsgd.R` builds a dense `matrix[N,N] prob_neigh` and double-loops all N² pairs on **every**
gradient evaluation to compute `mu[i] = x[i]·β + η·Σ_{j∈N(i)} y[j]`.

But during fitting the neighbor labels `y[j]` are **observed data** (conditional pseudolikelihood),
so the neighbor sum

```
c_i = Σ_{j ∈ N(i)} y_j            (label-matched for multi-slice)
```

is a **constant vector**. The N² loop recomputes a constant every iteration. Precompute `c` once
(O(N·k), k≈4 on a lattice) and the fit becomes `mu = Xβ + η·c (+U[label])` — **O(N) per iteration,
same posterior.**

## Files

| file | role |
|---|---|
| `neighbors.R` | `neighbor_sum()` (grid-hash, O(N·k), dependency-free) + `build_adjacency()` (sparse `A`) |
| `dsgd_sparse.R` | `dsgd_sparse_single/_multiple` (VI or NUTS) + `predict_sparse` (O(N·k·sweep) Gibbs) |
| `verify_her2.R` | equivalence on her2st section B1 |
| `bench_one.R` / `run_benchmark.sh` / `plot_scaling.R` | dense-vs-sparse scaling benchmark |

## Result 1 — equivalence (her2st B1, 269 spots)

- **Neighbor engine exact:** grid-hash `neighbor_sum` == brute-force `dist≤1`, `max|Δc| = 0`.
- **Same posterior (NUTS, exact):** dense vs sparse → `max|Δβ| = 0.062`, `cor(β) = 0.9996`,
  top-10 gene ranking `10/10`, `η` 0.154 vs 0.157. **PASS.**
- **Finding (bonus):** under **VI**, the dense model's N²-repeated uniform priors inflate the ELBO
  magnitude (~1.5e5), so VI's *relative* tolerance triggers **premature convergence** — the released
  dense fit is under-optimized. The sparse ELBO (~1.2e2) optimizes honestly. (Hence equivalence is
  shown via NUTS, which is free of this artifact.)

## Result 2 — scaling (synthetic lattice, VI fit time & peak RSS)

_(filled in by `run_benchmark.sh` → `scaling.csv`; see `fig_scaling_time.png`, `fig_scaling_mem.png`)_

- Per-gradient-eval cost already differed **~80×** at N=269 (0.0062 s dense vs 7.7e-5 s sparse);
  this gap grows ∝ N.
- Dense is intractable beyond a few thousand spots (N² matrix + autodiff graph); sparse reaches
  N=10⁵ on a laptop.

<!-- SCALING_TABLE -->

Measured (synthetic lattice, VI fit time; `scaling.csv`):

| N | dense-ADVI O(N²) | sparse-ADVI O(N) | speedup |
|---|---|---|---|
| 500    | 49.1 s    | 0.36 s | 136× |
| 1000   | 172 s     | 0.48 s | 358× |
| 2000   | 628 s     | 0.69 s | 910× |
| 3000   | 1347 s (22 min) | 0.92 s | **1464×** |
| 10000  | —         | 2.6 s  | dense infeasible |
| 100000 | — (~23 days projected) | 27.6 s | dense infeasible |

Dense grows ∝ N² (≈4× per doubling); sparse ∝ N. The PG-CAVI engine (`method/`) is
faster still — 300,000 spots in ~10 s (see `repro/sim/sim2`).

## Count-likelihood readiness (M2 hook)

The reformulated mean `mu = Xβ + η·c` leaves the covariate layer fully modular: swapping `log(X+1)`
for a negative-binomial / ZINB measurement model with size factors (M2) only changes how `X` enters,
not the spatial term or the O(N) inference. That is the intended path to fold Direction 2 (count
likelihood) into the scalable backbone.

## Reproduce

```bash
Rscript repro/m0/verify_her2.R          # equivalence (neighbor check + NUTS)
bash    repro/m0/run_benchmark.sh        # scaling.csv (dense N≤3000, sparse →1e5)
Rscript repro/m0/plot_scaling.R          # scaling figures
```
