#!/usr/bin/env bash
# Addendum A of docs/preregistration_trees_xai.md. CMH is cheap; FIGS/Boruta share the CPU with the other runners.
set -u
cd "$(dirname "$0")/../src"
PY=../.venv/bin/python
run () { echo "=== $* ($(date +%H:%M))"; "$@" || echo "!!! FAILED: $*"; }

for t in da ox pw; do run $PY 48_alternatives.py cmh --task $t; done
run $PY 47_trees_xai.py cv --task da --model figs --jobs 2
run $PY 48_alternatives.py boruta --task da --jobs 2
run $PY 47_trees_xai.py cv --task ox --model figs --jobs 2
run $PY 48_alternatives.py boruta --task ox --jobs 2
run $PY 47_trees_xai.py cv --task pw --model figs --jobs 2
# consensus and final metrics need the RF XAI outputs from run_trees_xai.sh
until grep -q TREES_XAI_DONE ../logs/run_trees_xai.log 2>/dev/null; do sleep 120; done
for t in da ox pw; do
  run $PY 47_trees_xai.py metrics --task $t          # re-run so FIGS is included
  run $PY 47_trees_xai.py xai --task $t --model figs
  run $PY 48_alternatives.py consensus --task $t
done
echo "ALTERNATIVES_DONE ($(date +%H:%M))"
