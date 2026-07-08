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
| Human Breast Cancer | 10x **Xenium** | breast tumor vs normal | 10x Genomics datasets portal |
| Mouse AD (13-mo) | **STARmap PLUS** | Alzheimer amyloid proximity | Zeng et al. 2023 (GSE194329) |
| Mouse/Human brain | **MERFISH** | AD / control | Vizgen / Allen Brain |
| Colorectal / Breast | **Visium HD** (8µm) | tumor region | 10x Genomics datasets portal |

## Disease-label construction
- Pathology annotations → cancer vs non-cancer (as in her2st).
- Amyloid-β proximity (STARmap+ AD): label a cell disease=1 if within 20µm of a plaque
  (the paper's rule); set the neighbor radius to ~the cell spacing (≈4 neighbors).
- Coordinates are snapped to a lattice at the neighborhood scale inside the pipeline so
  the rook-neighbor engine applies; adjust the divisor (`/20`) to your cell spacing.
