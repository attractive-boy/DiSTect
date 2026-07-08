#!/bin/bash
# Verification re-run: re-execute every previously-skipped script and record outputs
# for comparison against the committed logs/CSVs. Correctness steps first; the timing
# benchmark runs LAST and effectively alone (sequential) so its wall-clock is clean.
set -u
cd /Users/licheng/Documents/DiSTect
V=repro/_verify
mkdir -p "$V"
PROG="$V/progress.txt"
: > "$PROG"

# preserve committed artifacts for side-by-side diff (working tree is about to be overwritten)
cp repro/m0/scaling.csv        "$V/scaling.committed.csv"        2>/dev/null || true
cp repro/sim/sim2_scaling.csv  "$V/sim2_scaling.committed.csv"  2>/dev/null || true

emit(){ echo "$1" | tee -a "$PROG"; }

step(){  # step <name> <cmd...>
  local name=$1; shift
  local t0 rc t1
  t0=$(date +%s)
  emit "STEP_START $name $(date +%H:%M:%S)"
  ( "$@" ) > "$V/$name.log" 2>&1
  rc=$?
  t1=$(date +%s)
  if [ $rc -eq 0 ]; then emit "STEP_DONE $name rc=0 $((t1-t0))s"
  else                   emit "STEP_FAIL $name rc=$rc $((t1-t0))s"; fi
}

emit "VERIFY_BEGIN $(date +%H:%M:%S)"

# ---- correctness / value checks (exact reproduction expected; seeds set) ----
step verify_her2    Rscript repro/m0/verify_her2.R
step singlecell     Rscript repro/realdata/run_singlecell.R
step nuts_compare   Rscript repro/nuts_compare.R
step reproduce      Rscript repro/reproduce.R
step her2_reproduce Rscript repro/her2_reproduce.R

# ---- timing benchmark (MUST be isolated; runs last, nothing else concurrent) ----
step benchmark      bash repro/m0/run_benchmark.sh

# ---- regenerate figures + sim2 merge from the fresh benchmark ----
step plot_scaling   Rscript repro/m0/plot_scaling.R
step sim2_scaling   Rscript repro/sim/sim2_scalability.R

emit "VERIFY_END $(date +%H:%M:%S)"
