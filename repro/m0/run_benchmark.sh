#!/bin/bash
# Dense-vs-sparse scaling benchmark. Each (method,N) runs as a fresh subprocess under
# /usr/bin/time -l so we capture true peak RSS (the dense NxN matrix lives in C++ and is
# invisible to R's gc()). Results append incrementally so partial runs survive.
set -u
cd /Users/licheng/Documents/DiSTect
OUT=repro/m0/scaling.csv
echo "method,N,fit_seconds,peak_rss_mb" > "$OUT"

run() {
  local method=$1 N=$2 tmp
  tmp=$(mktemp)
  /usr/bin/time -l Rscript repro/m0/bench_one.R "$method" "$N" >"$tmp" 2>&1
  local secs rss
  secs=$(grep RESULT "$tmp" | sed -E 's/.*fit_seconds=([0-9.]+).*/\1/')
  rss=$(grep "maximum resident set size" "$tmp" | awk '{printf "%.1f", $1/1048576}')
  echo "$method,$N,${secs:-NA},${rss:-NA}" >> "$OUT"
  echo "[done] $method N=$N  fit=${secs:-NA}s  peakRSS=${rss:-NA}MB"
  rm -f "$tmp"
}

# sparse first: O(N) headline curve, pushes to 1e5 (regime dense cannot reach)
for N in 500 1000 2000 3000 10000 30000 100000; do run sparse "$N"; done
# dense: O(N^2) contrast -- capped at 3000 (N^2 autodiff graph OOMs beyond a few thousand)
for N in 500 1000 2000 3000; do run dense "$N"; done

echo "BENCHMARK DONE -> $OUT"
