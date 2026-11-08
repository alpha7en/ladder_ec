#!/bin/bash
# Run every AutoQ-Para experiment of autoq_models three times; wall time and peak memory by /usr/bin/time (Table 9).
# Usage: bash experiments/autoq_para/run_autoq.sh BUILD_DIR > data/runs/autoq_para/autoq_runs.txt
B=$1
for e in $($B/autoq_models); do
  for r in 1 2 3; do
    /usr/bin/time -f "## wall=%es maxRSS=%MkB" $B/autoq_models $e 2>&1 | tr '\n' ' '; echo
  done
done
