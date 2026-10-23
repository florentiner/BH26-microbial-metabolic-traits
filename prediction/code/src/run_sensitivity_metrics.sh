#!/bin/sh
# Metrics for all sensitivity analyses (3 in parallel). usage: sh src/run_sensitivity_metrics.sh
cd "$(dirname "$0")"
PY=../.venv/bin/python
printf '%s\n' "--label-mode block_nan" "--split gtdb" "--split seed1" "--split seed2" "--split seed3" \
  "--split seed4" "--split seed5" "--subset quality_ok" | xargs -P 3 -I{} sh -c "$PY 08_metrics.py {} > /dev/null 2>&1 || echo FAILED {}"
echo metrics-done
