#!/bin/bash
# Fetch the GOLD-STANDARD supervised cell-type labels for Xenium Human Breast Cancer Rep1
# (Janesick et al., Nat Commun 2023) and export cell_groups.csv (Barcode,celltype) for use as
# LABEL_CSV in make_xenium_rds.R. The labels are the 'Fig. 3e-j Xenium' sheet of the paper's
# Supplementary Data (MOESM4) -- label-transferred from a matched scFFPE/FLEX reference, i.e. an
# INDEPENDENT annotation (removes the circularity of the marker-score / cluster proxies).
set -e
DEST=${1:-/tmp/xenium_breast}
XLSX="$DEST/cell_type_annot.xlsx"
URL="https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-43458-x/MediaObjects/41467_2023_43458_MOESM4_ESM.xlsx"
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$DEST"

# download (resumable) unless a complete, valid xlsx is already present (~37 MB)
if [ -s "$XLSX" ] && python3 -c "import zipfile,sys;zipfile.ZipFile(sys.argv[1]).testzip()" "$XLSX" 2>/dev/null; then
  echo "  have valid $XLSX"
else
  echo "  downloading Supplementary Data (MOESM4, ~37 MB) ..."
  curl -fSL -C - --max-time 1800 -A "Mozilla/5.0" -o "$XLSX" "$URL"
fi

# ensure a Python with openpyxl (use system if present, else an isolated venv -- no system changes)
PYBIN=python3
if ! $PYBIN -c "import openpyxl" 2>/dev/null; then
  if [ ! -x "$DEST/xlsxenv/bin/python" ]; then
    echo "  creating isolated venv for openpyxl ..."
    python3 -m venv "$DEST/xlsxenv"; "$DEST/xlsxenv/bin/pip" install -q openpyxl
  fi
  PYBIN="$DEST/xlsxenv/bin/python"
fi

$PYBIN "$HERE/parse_xenium_labels.py" "$XLSX" "$DEST/cell_groups.csv"
echo "Next: LABEL_CSV=$DEST/cell_groups.csv HVG=50 OUT_RDS=\$DEST/xenium_supervised.rds \\"
echo "        Rscript repro/realdata/make_xenium_rds.R"
