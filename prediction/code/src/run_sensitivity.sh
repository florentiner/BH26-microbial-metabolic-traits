#!/bin/sh
# Sensitivity analyses S1-S4 (docs/deviations.md). usage: sh src/run_sensitivity.sh "nulls A B1 B1r B2" [jobs]
set -e
ARMS="$1"; JOBS="${2:-3}"
cd "$(dirname "$0")"
PY=../.venv/bin/python
for arm in $ARMS; do
  $PY 06_run_arms.py "$arm" --jobs "$JOBS" --label-mode block_nan          # S1
  $PY 06_run_arms.py "$arm" --jobs "$JOBS" --split gtdb                    # S2
  for s in 1 2 3 4 5; do $PY 06_run_arms.py "$arm" --jobs "$JOBS" --split seed$s; done   # S3
done
