# Single-cell / imaging-based spatial datasets for the generalization experiments

The upgraded method targets single-cell-resolution platforms, where the original
O(N²) log-Gaussian model is infeasible. `run_singlecell.R` runs out-of-the-box on a
synthetic Xenium-like dataset; to use real data, prepare an `.rds` and set `DATA_RDS`.

## Expected `.rds` format
```r
saveRDS(list(
  counts  = <cells x genes integer matrix>,   # raw UMI / molecule counts
  coords  = <cells x 2 matrix>,                # x,y spatial coordinates (any units)
  disease = <length-cells 0/1 vector>          # disease label per cell
), "mydata.rds")
DATA_RDS=mydata.rds Rscript repro/realdata/run_singlecell.R
```

## Recommended public datasets (match the paper's platforms)
| dataset | platform | disease | source |
|---|---|---|---|
| Human Breast Cancer | 10x **Xenium** | breast tumor vs normal | 10x Genomics portal — **✅ scripted & verified**, see below |
| Mouse AD (13-mo) | **STARmap PLUS** | Alzheimer amyloid proximity | Broad Single Cell Portal **SCP1375** (Zeng et al., *Nat Neurosci* 2023); **sign-in required** |
| Mouse/Human brain | **MERFISH** | AD / control | Vizgen / Allen Brain |
| Colorectal / Breast | **Visium HD** (8µm) | tumor region | 10x Genomics datasets portal |

> ⚠️ **Correction:** an earlier version of this table cited **GSE194329** for the STARmap+ AD data.
> That accession is a *different* study (Visium of pediatric DMG/GBM brain tumors) and is **not** the
> Alzheimer's STARmap+ dataset. The STARmap+ AD data is distributed via the Broad Single Cell Portal
> (SCP1375) and a restricted Zenodo record (5842625) — not GEO, and both require access.

## Ready-to-run: 10x Xenium Human Breast Cancer (verified end-to-end)
Single-cell-resolution imaging spatial data (Janesick et al., *Nat Commun* 2023; CC BY 4.0) —
167,780 cells at low depth, the regime the original O(N²) log-Gaussian model cannot fit.
```bash
bash   repro/realdata/fetch_xenium_breast.sh          # ~118 MB (skips the 9.4 GB full bundle)
HVG=50 Rscript repro/realdata/make_xenium_rds.R       # -> /tmp/xenium_breast/..._hvg50.rds
DATA_RDS=/tmp/xenium_breast/xenium_breast_rep1_hvg50.rds Rscript repro/realdata/run_singlecell.R
# -> 167,780 cells x 50 genes fit in ~100s; 42 genes @FDR 0.05; EPCAM/CDH1/KRT8/FASN/TACSTD2/...
```
- `make_xenium_rds.R` keeps the full 313-gene panel by default; set `HVG=K` — the fit is
  **O(N·P²)/iter**, so the full panel hangs at N=1.7e5 (as the README warns: reduce HVGs).
- **Disease labels — three routes** (increasing rigor, all consumed by `make_xenium_rds.R` via `LABEL_CSV`):

  | route | script | circularity | fit |
  |---|---|---|---|
  | marker-score proxy (default) | built into `make_xenium_rds.R` | high (label ≈ covariate threshold) | eta 0.41; top EPCAM/CDH1/KRT8 |
  | unsupervised cluster + markers | `build_xenium_labels.R` | low | eta 0.66; ↑GATA3/FOXA1, ↓PTPRC/TRAC |
  | **gold-standard supervised** | `fetch_xenium_labels.sh` | **none (independent classifier)** | **eta 0.49, AUC 0.995** |

### Gold-standard supervised labels (verified end-to-end)
The Janesick label-transferred Xenium cell types are the paper's Supplementary Data (**MOESM4**),
sheet **`Fig. 3e-j Xenium`** (`Barcode` = integer cell_id, `Cluster` = supervised cell type). Fetch → label → discover:
```bash
bash repro/realdata/fetch_xenium_labels.sh   # downloads MOESM4 xlsx, exports cell_groups.csv (167,780 cells, 38% tumor)
LABEL_CSV=/tmp/xenium_breast/cell_groups.csv HVG=50 \
  OUT_RDS=/tmp/xenium_breast/xenium_supervised.rds Rscript repro/realdata/make_xenium_rds.R
DATA_RDS=/tmp/xenium_breast/xenium_supervised.rds Rscript repro/realdata/xenium_discovery_report.R
```
Result (167,780 cells × 50 HVG, **converged 442 iters, ~5 min**): **eta=0.485, AUC=0.995, 46/50 genes @FDR 0.05.**
Tumor-**enriched** = luminal program (FOXA1, GATA3, MLPH, EPCAM, KRT8, FASN); tumor-**depleted** =
CAF/stromal + immune (MMP2, POSTN, CXCL12, CCDC80, CD4) — a coherent tumor-vs-microenvironment contrast
at single-cell resolution the original O(N²) model cannot produce. `fetch_xenium_labels.sh` needs
`openpyxl` (auto-creates an isolated venv if absent).

## Disease-label construction
- Pathology annotations → cancer vs non-cancer (as in her2st).
- Amyloid-β proximity (STARmap+ AD): label a cell disease=1 if within 20µm of a plaque
  (the paper's rule); set the neighbor radius to ~the cell spacing (≈4 neighbors).
- Coordinates are snapped to a lattice at the neighborhood scale inside the pipeline so
  the rook-neighbor engine applies; adjust the divisor (`/20`) to your cell spacing.
