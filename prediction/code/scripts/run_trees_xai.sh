#!/usr/bin/env bash
# Tree models + XAI (docs/preregistration_trees_xai.md), cheapest first. Shares the CPU with run_pathways.sh.
set -u
cd "$(dirname "$0")/../src"
PY=../.venv/bin/python
J="${JOBS:-2}"
run () { echo "=== $* ($(date +%H:%M))"; "$@" || echo "!!! FAILED: $*"; }

for m in rf et l1; do run $PY 47_trees_xai.py cv --task da --model $m --jobs $J; done
run $PY 47_trees_xai.py metrics --task da
for m in dt rf l1 et; do run $PY 47_trees_xai.py xai --task da --model $m; done

run $PY 47_trees_xai.py nulls --task ox
for m in dt l1 rf et; do run $PY 47_trees_xai.py cv --task ox --model $m --jobs 2; done
run $PY 47_trees_xai.py metrics --task ox
for m in dt rf l1; do run $PY 47_trees_xai.py xai --task ox --model $m; done

for m in dt l1 rf et; do run $PY 47_trees_xai.py cv --task pw --model $m --jobs $J; done
run $PY 47_trees_xai.py metrics --task pw
for m in dt rf l1; do run $PY 47_trees_xai.py xai --task pw --model $m; done
echo "TREES_XAI_DONE ($(date +%H:%M))"
