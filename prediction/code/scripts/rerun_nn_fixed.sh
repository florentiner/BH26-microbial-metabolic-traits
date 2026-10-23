#!/usr/bin/env bash
# Re-run the corrected neural network (zero-variance fix) on the Apple GPU, then both attribution variants.
set -u
cd "$(dirname "$0")/../src"
PY=../.venv/bin/python
run () { echo "=== $* ($(date +%H:%M))"; "$@" || echo "!!! FAILED: $*"; }
run $PY 47_trees_xai.py cv --task da --model nn --jobs 1
run $PY 47_trees_xai.py metrics --task da
NN_ATTRIB=ig_zero run $PY 47_trees_xai.py xai --task da --model nn        # pre-registered (Addendum C)
NN_ATTRIB=ig_zero run $PY 48_alternatives.py consensus --task da
run $PY 47_trees_xai.py xai --task da --model nn                          # Expected Gradients (Addendum F)
run $PY 48_alternatives.py consensus --task da
run $PY 47_trees_xai.py cv --task ox --model nn --jobs 1
echo "NN_FIXED_DA_DONE ($(date +%H:%M))"
