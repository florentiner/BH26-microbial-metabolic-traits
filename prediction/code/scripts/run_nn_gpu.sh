#!/usr/bin/env bash
# Addendum C: neural network + Integrated Gradients on the Apple GPU (MPS). One process owns the GPU.
set -u
cd "$(dirname "$0")/../src"
PY=../.venv/bin/python
run () { echo "=== $* ($(date +%H:%M))"; "$@" || echo "!!! FAILED: $*"; }
run $PY 47_trees_xai.py cv --task da --model nn --jobs 1
run $PY 47_trees_xai.py cv --task ox --model nn --jobs 1
# metrics/XAI after the CPU runners, so every model is in the same tables
until grep -q ALTERNATIVES_DONE ../logs/run_alternatives.log 2>/dev/null; do sleep 120; done
for t in da ox; do
  run $PY 47_trees_xai.py metrics --task $t
  run $PY 47_trees_xai.py xai --task $t --model nn
done
for t in da ox; do run $PY 48_alternatives.py consensus --task $t; done   # adds RF-NN agreement
echo "NN_GPU_DONE ($(date +%H:%M))"
