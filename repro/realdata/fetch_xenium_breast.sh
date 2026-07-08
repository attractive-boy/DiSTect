#!/bin/bash
# Fetch the 10x Xenium FFPE Human Breast Cancer Rep1 dataset (Janesick et al., Nat Commun 2023),
# CC BY 4.0 -- https://www.10xgenomics.com/datasets. Single-cell-resolution imaging spatial data:
# the regime the upgraded (O(N), count-aware) DiSTect targets and the original O(N^2) model cannot.
#
# Downloads only the 3 small files needed to build the .rds (skips the 9.4 GB full outs.zip):
#   - cell_feature_matrix.tar.gz  (~49 MB)  raw per-cell x gene counts (10x MTX)
#   - cells.csv.gz                (~7.6 MB) per-cell x_centroid,y_centroid (+ QC)
#   - analysis.tar.gz             (~61 MB)  graph-based clusters (optional, for labeling)
set -e
DEST=${1:-/tmp/xenium_breast}
BASE="https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1"
PFX="Xenium_FFPE_Human_Breast_Cancer_Rep1"
mkdir -p "$DEST"; cd "$DEST"

fetch() {  # <filename>
  if [ -s "$1" ]; then echo "  have $1"; return; fi
  echo "  downloading $1 ..."
  curl -fL --retry 3 --connect-timeout 30 -o "$1" "$BASE/$1"
}
fetch "${PFX}_cell_feature_matrix.tar.gz"
fetch "${PFX}_cells.csv.gz"
fetch "${PFX}_analysis.tar.gz"

# unpack matrix (-> cell_feature_matrix/) and analysis (-> analysis/)
[ -d cell_feature_matrix ] || tar xzf "${PFX}_cell_feature_matrix.tar.gz"
[ -d analysis ]            || tar xzf "${PFX}_analysis.tar.gz"

echo "Fetched Xenium Breast Rep1 into $DEST"
echo "Next: Rscript repro/realdata/make_xenium_rds.R   # -> \$DEST/xenium_breast_rep1.rds"
