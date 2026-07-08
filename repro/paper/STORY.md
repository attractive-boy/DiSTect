# The Story — a Bioinformatics paper built on DiSTect

**Working title:** *scDiSTect: scalable, count-aware Bayesian discovery of disease-associated
genes in single-cell spatial transcriptomics.*

This document is the complete narrative and the evidence already in hand. Everything marked
✅ is measured by a script in this repo; ⏳ marks the experiments left for you to run on real
single-cell data (pipelines and instructions provided).

---

## 1. The one-sentence thesis

> DiSTect (Zhao et al., *Bioinformatics* 2025) is the right *model* for disease-associated gene
> discovery in spatial transcriptomics, but its inference is **O(N²)** and its **log-Gaussian**
> covariate layer breaks at single-cell resolution. We make it **O(N)** (sparse conditional
> pseudolikelihood + Pólya-Gamma closed-form VI), **count-aware** (negative-binomial size-factor
> normalization), and **statistically honest** (full-likelihood η de-biasing, FDR-controlled
> selection, calibrated prediction) — turning a 25-hour, 4-patient analysis into a **4-second,
> all-8-patient** one that also runs on **10⁵–10⁶-cell** Xenium/MERFISH/Visium-HD data the
> original cannot touch.

Why this is one story, not three: single-cell platforms are simultaneously **large N** (needs
O(N) inference) and **low depth** (needs a count likelihood). The two contributions are forced by
the *same* data regime.

---

## 2. The gap (why a reviewer cares)

Spatial transcriptomics is moving from Visium (~300 spots/slice) to **Xenium / MERFISH / Visium HD**
(10⁵–10⁶ cells). DiSTect's spatial autologistic model rebuilds a dense `N×N` neighbor matrix on
every gradient step — **O(N²) compute and autodiff memory**. Consequences we measured on the
published code:
- ✅ The paper's own HER2 analysis took **25 h**; our faithful reproduction had to be cut from 8
  patients to 4 to finish (`repro/her2_reproduce.R`).
- ✅ At N=3000 the dense fit takes **22 minutes**; at N=10⁵ it is projected at **~23 days**
  (`repro/m0/scaling.csv`).

And the covariate layer, `log(X+1)` + Gaussian, is invalid where counts are low and sparse — exactly
single-cell resolution.

No existing tool fills the four-way gap **{disease-conditional multi-gene selection} × {count
likelihood} × {single-cell scalability} × {prediction + missingness}**: SVG methods (SpatialDE,
SPARK-X, nnSVG, Celina, PROST) are not disease-conditional; C-SIDE is count-aware and cell-type DE
but marginal per-gene and non-predictive; BOOST-GP is Bayesian but does not scale.

---

## 3. Contributions and the evidence in hand

**C1 — O(N) inference, identical model.**
The neighbor labels are observed during fitting, so the neighbor sum `c_i = Σ_{j∈N(i)} y_j` is a
**constant vector**; the O(N²) loop recomputes a constant. Precompute it → `mu = Xβ + η·c`.
- ✅ Grid-hash neighbor sum == dense `dist≤1` rule exactly (`verify_her2.R`, max|Δ|=0).
- ✅ Sparse ≡ dense posterior via NUTS: max|Δβ|=0.062, **cor=0.9996**, top-10 ranking 10/10.
- ✅ Scaling: **1464× at N=3000**; sparse reaches N=10⁵ in 27 s (`repro/m0/`).

**C2 — Pólya-Gamma closed-form VI.**
Because the fit is logistic regression in `[X, c]`, PG augmentation gives conditionally-Gaussian,
closed-form CAVI (Durante & Rigon 2019) — no HMC/ADVI gradient ascent, monotone ELBO, and it
sidesteps the ADVI pathologies we observed (Pareto-k=1.1; the dense model's N²-inflated ELBO even
triggers *premature* convergence).
- ✅ Recovers β (cor 0.99) and matches `glm` on η; **N=1600 in 0.05 s** (`sim1`).
- ✅ **300,000 cells in ~10 s**; 14,600× vs dense-ADVI at N=3000 (`sim2`).

**C3 — Count-aware covariates (NB size-factor normalization).**
`log(X+1)` is heteroscedastic and depth-confounded at low depth; NB Pearson residuals fix both.
- ✅ mean-variance coupling drops from **0.96 → 0.4**; residual depth confound from **0.95 → 0.02**
  (`sim3`), most pronounced at low depth.

**C4 — Statistical honesty.**
- ✅ **η de-biasing**: pseudolikelihood attenuates η (true 2.0 → 0.76); a full-likelihood
  moment-matching corrector moves it back toward truth (`sim4`, `method/eta_fulllik.R`).
- ✅ **FDR-controlled selection** (BH / Bayesian FDR) replacing the ad-hoc |mean/sd|>1.96 cut
  (`method/selection_fdr.R`).
- ✅ **Calibrated prediction** (posterior-predictive + Brier/ECE) instead of plug-in means
  (`method/predict_calibrated.R`).
- Also fixes the boundary η-normalization discrepancy (code uses `η·Σy`, paper says `η/|N(i)|·Σy`).

**C5 — Real-data demonstration.**
- ✅ HER2+ (her2st), **all 8 patients, 3157 spots, fit in 3.79 s** (vs ~25 h): 26 genes at
  **q<0.10** (TIMP1, COL6A2, mitochondrial/Warburg genes …), η=3.08 (debiased 3.13), LOPO
  **accuracy 0.740, Brier 0.180, ECE 0.179** (`realdata/run_her2st.R`).
- ✅ Synthetic Xenium-scale: **62,500 cells in 22.6 s**, correct gene ranking (`run_singlecell.R`).
- ⏳ Real Xenium / MERFISH / Visium-HD generalization (pipeline + data guide provided).

---

## 4. Figure plan (maps to scripts)

| Fig | content | script | status |
|---|---|---|---|
| 1 | method schematic: dense→sparse c-vector, PG-CAVI, NB layer | (draw) | — |
| 2 | Sim1 correctness (β recovery, FDR TPR/FPR) | `sim/sim1_correctness.R` | ✅ |
| 3 | Sim2 scaling: dense-ADVI vs sparse-ADVI vs PG-CAVI | `sim/sim2_scalability.R` → `fig_sim2_scaling.png` | ✅ |
| 4 | Sim3 count regime: variance stabilization + depth invariance | `sim/sim3_count_regime.R` | ✅ |
| 5 | Sim4 η attenuation + de-biasing | `sim/sim4_pseudolik_bias.R` | ✅ |
| 6 | HER2 all-8: gene coef + q-values, η, LOPO calibration | `realdata/run_her2st.R` | ✅ |
| 7 | single-cell (Xenium/MERFISH): discovery + runtime at 10⁵–10⁶ | `realdata/run_singlecell.R` | ⏳ real data |

---

## 5. Positioning / related work (differentiation table)

| method | disease-conditional | multi-gene joint | count model | scales to 10⁵⁺ | predictive |
|---|---|---|---|---|---|
| SpatialDE / SPARK-X / nnSVG | ✗ (SVG) | ✗ | partial | ✓ | ✗ |
| C-SIDE | ✓ | ✗ (marginal) | ✓ | ~ | ✗ |
| BOOST-GP | ✓ | ✓ | ✓ | ✗ | ✗ |
| DiSTect (original) | ✓ | ✓ | ✗ | ✗ (O(N²)) | ✓ |
| **this work** | ✓ | ✓ | ✓ | ✓ | ✓ (calibrated) |

---

## 6. What is left to make it Bioinformatics-complete

1. ⏳ Run C5 on ≥1 real single-cell disease dataset (Xenium breast or STARmap+ AD) — the headline
   generalization result. Pipeline + data guide: `realdata/run_singlecell.R`, `SINGLECELL_DATA.md`.
2. Layer spike-and-slab onto PG-CAVI (documented hook in `fit_polyagamma.R`) for principled
   inclusion probabilities → Bayesian FDR (`select_fdr_bayes`).
3. Full latent-NB measurement model (M2 headline; residual version already shipped).
4. Multi-slice random effects in PG-CAVI (batch adjustment) — currently pooled with patient labels.
5. Ablations + external validation; write methods + theory (pseudolikelihood consistency, PG-CAVI
   ELBO monotonicity), assemble figures.

## 7. Honest limitations (put these in the paper, reviewers will find them otherwise)
- Mean-field VI **under-estimates posterior variance** → slightly inflated selection (Sim1 FPR≈0.19);
  the spike-and-slab + calibration layers address this.
- The η de-biaser is a Monte-Carlo moment-matcher: reduces but does not fully remove attenuation at
  large η; needs more MC draws / a better estimating equation.
- Single-cell coordinates are snapped to a lattice for the rook-neighbor engine; a kd-tree
  fixed-radius graph (`dbscan::frNN`) is the general-coordinate replacement (M1).
