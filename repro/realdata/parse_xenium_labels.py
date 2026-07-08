#!/usr/bin/env python3
"""Extract the Janesick et al. 2023 (Nat Commun) supervised Xenium Rep1 cell types.

Reads the 'Fig. 3e-j Xenium' sheet of the paper's Supplementary Data (MOESM4 xlsx) and writes
cell_groups.csv (Barcode,celltype). Barcode is the integer Xenium cell_id, so it aligns directly
with the cell_feature_matrix. These are label-transferred labels from a matched scFFPE/FLEX
reference -- an INDEPENDENT annotation, not derived from the covariates.

Usage: parse_xenium_labels.py <in.xlsx> <out.csv>
"""
import sys, csv, collections, openpyxl

XLSX = sys.argv[1] if len(sys.argv) > 1 else "/tmp/xenium_breast/cell_type_annot.xlsx"
OUT  = sys.argv[2] if len(sys.argv) > 2 else "/tmp/xenium_breast/cell_groups.csv"
SHEET = "Fig. 3e-j Xenium"

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
if SHEET not in wb.sheetnames:
    sys.exit(f"sheet '{SHEET}' not found; sheets = {wb.sheetnames}")
ws = wb[SHEET]
it = ws.iter_rows(values_only=True)
hdr = list(next(it))
bi, ci = hdr.index("Barcode"), hdr.index("Cluster")

n = 0; cts = collections.Counter()
with open(OUT, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["Barcode", "celltype"])
    for r in it:
        if r[bi] is None:
            continue
        w.writerow([r[bi], r[ci]]); cts[r[ci]] += 1; n += 1

tumor = sum(v for k, v in cts.items() if any(t in str(k) for t in ("Tumor", "DCIS", "Invasive")))
print(f"wrote {OUT}: {n} cells, {len(cts)} cell types")
print(f"tumor cells (Invasive_Tumor/DCIS/Prolif): {tumor} ({100*tumor/n:.0f}%)")
