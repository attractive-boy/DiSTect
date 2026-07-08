#!/bin/bash
# Fetch the labeled HER2+ breast sections from almaan/her2st (Andersson 2021) without
# the large H&E images: sparse, blob-filtered clone + checkout of only the needed files.
# One pathologist-labeled section exists per patient (A1,B1,C1,D1,E1,F1,G2,H1 = A-H).
set -e
DEST=${1:-/tmp/her2st}
if [ -d "$DEST/.git" ]; then echo "her2st already at $DEST"; exit 0; fi
GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 --filter=blob:none --no-checkout \
  https://github.com/almaan/her2st.git "$DEST"
cd "$DEST"
FILES=""
for s in A1 B1 C1 D1 E1 F1 G2 H1; do
  FILES="$FILES data/ST-cnts/$s.tsv.gz data/ST-spotfiles/${s}_selection.tsv \
         data/ST-pat/lbl/${s}_labeled_coordinates.tsv"
done
git checkout HEAD -- $FILES
echo "Fetched labeled sections into $DEST/data (counts + spotfiles + labels)."
