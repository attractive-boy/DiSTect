# Original simulation script search

Date: 2026-07-09

Goal: find the authors' original scripts or historical code used to generate the
btaf530 Simulation 1/2/3 tables, especially to explain the systematic eta-bias
mismatch.

## Bottom line

No original paper Simulation 1/2/3 script was found in the public sources
advertised by the paper or in the visible GitHub history/forks checked below.

The public artifacts contain the DiSTect/BayModDSGD package, toy data, and
tutorial material only. They do not contain the scripts that generated Tables
C1-C6 or Figures B2-B10.

## Checked sources

### Paper-declared GitHub repository

Paper availability statement:

- `https://github.com/StaGill/DiSTect`

GitHub API tree for `StaGill/DiSTect` main contains only:

- `.github/workflows/r.yml`
- `DESCRIPTION`
- `DiSTect.zip`
- `NAMESPACE`
- `R/coef_plot.R`
- `R/dsgd.R`
- `R/missing_imputation.R`
- `R/network_plot.R`
- `R/prediction.R`
- `README.md`
- `data/toy.rda`
- `data/toy2.rda`
- man pages

No `simulation`, `sim`, `paper`, `reproduce`, `HER2`, or `STARmap` source files
are present.

### Paper-declared Zenodo record

Paper availability statement:

- `https://zenodo.org/records/17127211`
- DOI `10.5281/zenodo.17127211`

Zenodo API record contains one file:

- `DiSTect-main.zip`, size `112517`, md5 `7750d2724a88ad8ab6c4cb9bcd2724bd`

Downloaded archive contents:

- `DiSTect-main/BayModDSGD_0.0.0.9000.tar.gz`
- `DiSTect-main/README.md`
- `DiSTect-main/Tutorial/Tutorial.Rmd`
- `DiSTect-main/Tutorial/Tutorial_BayModDSGD.pdf`

This is an early repository/package snapshot, not the simulation source.

### Official fork

GitHub lists one fork of `StaGill/DiSTect`:

- `https://github.com/AnjiDeng/DiSTect`

Its main tree contains only:

- `BayModDSGD_0.0.0.9000.tar.gz`
- `DiSTect.zip`
- `README.md`
- `Tutorial/Tutorial.Rmd`
- `Tutorial/Tutorial_BayModDSGD.pdf`

No simulation scripts are present.

### Repository history

The visible history of the local mirror / official repository includes:

- early `mypkg/BayModDSGD` package skeleton
- `Tutorial/Tutorial.Rmd`, `Tutorial.html`, `Tutorial_BayModDSGD.pdf`
- several uploads/deletions of `BayModDSGD_0.0.0.9000.tar.gz`
- later upload of `DiSTect.zip`
- later package files under `R/`, `data/`, `man/`

Historical file listing did not reveal Simulation 1/2/3 scripts. The package
tarballs were inspected and contain only package source, toy data, and man
pages.

An all-history `git grep` across visible commits was also run for simulation
markers such as `avgBias`, `Simulation 1`, `Simulation 2`, `Simulation 3`,
`Gibbs`, `eta = 0.4`, `sigma2`, `rho`, `gamma`, and `nonignorable`. No paper
simulation runner or result-generation script was found.

The later `DiSTect.zip` archive contains `.Rhistory` files. The top-level
`.Rhistory` was inspected; it contains unrelated MALDI/Seurat clustering work
and package tutorial commands on `toy`/`toy2`, but no paper simulation code.

### Tutorial page

Tutorial URL from README:

- `https://qihuangzhang.github.io/software/DiSTect_tutorial`

The page demonstrates:

- installation
- toy single-slice analysis
- toy multiple-slice analysis
- interaction analysis
- missing-data imputation

It does not include the paper simulation study scripts.

### GitHub issues, PRs, releases, pages, downloads

Checked:

- `StaGill/DiSTect` issues and PRs
- `StaGill/DiSTect` releases/downloads/pages/wiki
- PR `#1` diff

Findings:

- only one closed PR, "Add files via upload", whose diff adds binary
  `DiSTect.zip`
- no releases
- no downloads/pages/wiki content

### Author / organization repositories

Checked public repositories for:

- `qc-zhao`
- `AnjiDeng`
- `StaGill`
- `attractive-boy`

No additional DiSTect simulation repository was found. `qc-zhao` has unrelated
public repositories; `AnjiDeng` has the DiSTect fork and unrelated Galaxy code;
`StaGill` has DiSTect and unrelated repositories.

### Zenodo keyword search

Zenodo searches for `DiSTect`, `BayModDSGD`, and `btaf530` found only the paper
declared `17127211` record for DiSTect. No additional simulation-code record was
found.

The Zenodo metadata for `17127211` reports:

- `conceptrecid`: `17127210`
- `conceptdoi`: `10.5281/zenodo.17127210`
- one version only under the concept DOI
- no `related_identifiers`
- no extra notes or linked records

So there is no separate public Zenodo version carrying simulation scripts.

### External code indexes

Sourcegraph searches with `fork:yes archived:yes` were attempted for:

- `BayModDSGD`
- `prob_neigh[i,j]=eta*y[j]`
- `DiSTect dsgd_single`
- `avgBias eta spatial transcriptomics`

No public indexed source file matching the paper simulation scripts was found.
`searchcode.com` API attempts returned `404 page not found`, so it was not usable.
Jina/Bing and Wayback CDX attempts did not produce useful script hits; Wayback
queries were slow/inconclusive and were stopped rather than treated as evidence
of absence.

## Relevant code discrepancy found in public package

The public package code uses an unnormalized neighbor sum:

```stan
prob_neigh[i,j] = eta * y[j];
mu[i] = dot_product(x[i,1:P], beta) + sum(prob_neigh[i,]);
```

The paper formula uses the neighbor average:

```text
eta / |N(i)| * sum_{j in N(i)} Y_j
```

This discrepancy is real in the public code, but prior diagnostics show that the
unnormalized convention alone does not recover the paper Simulation 1 eta
targets across all eta settings.

## Current conclusion

The authors' exact Simulation 1/2/3 generation and fitting scripts do not appear
to be publicly available in the paper-declared GitHub/Zenodo resources as of this
search, nor in the broader public code-index checks that completed successfully.
To resolve the eta mismatch, the next practical step is to contact the authors
for the exact scripts used to produce Tables C1-C6/Figures B2-B10, or to ask them
to clarify whether the simulations used a different neighbor-scaling, Gibbs
generation, posterior-summary, or bias-reporting convention from the paper text.
