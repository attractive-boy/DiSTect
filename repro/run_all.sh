#!/bin/bash
# Run the full analysis kit end-to-end (simulations + real-data pipelines).
# The M0 scaling benchmark (~40 min) is NOT included here; run it separately:
#   bash repro/m0/run_benchmark.sh && Rscript repro/m0/plot_scaling.R
set -e
cd "$(dirname "$0")/.."          # repo root
R() { echo; echo "########## $1 ##########"; Rscript "$1"; }

# real data (HER2) -- fetch if missing
[ -d /tmp/her2st/data ] || bash repro/realdata/fetch_her2st.sh

echo "===== SIMULATIONS ====="
R repro/sim/sim1_correctness.R
R repro/sim/sim3_count_regime.R
R repro/sim/sim4_pseudolik_bias.R
R repro/sim/sim2_scalability.R          # merges m0/scaling.csv if present

echo; echo "===== REAL DATA ====="
R repro/realdata/run_her2st.R
R repro/realdata/run_singlecell.R       # synthetic fallback unless DATA_RDS set

echo; echo "ALL DONE. See repro/paper/STORY.md for the narrative and repro/paper/MANUAL.md to extend."
