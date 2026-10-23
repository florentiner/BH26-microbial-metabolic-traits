#!/usr/bin/env bash
# Addendum M: genome -> frozen Bacformer-large -> MLP head. CPU (tiny network: faster than MPS), traits in parallel.
set -u
cd "$(dirname "$0")/../src"
PY=../.venv/bin/python
export NN_DEVICE=cpu FM_THREADS=2
run () { echo "=== $* ($(date +%H:%M))"; "$@" || echo "!!! FAILED: $*"; }
run $PY 59_fm_mlp.py cv --arm FM_BFL_MLP --jobs 4
run $PY 59_fm_mlp.py external --arm FM_BFL_MLP
run $PY 59_fm_mlp.py summary
echo "FM_MLP_DONE ($(date +%H:%M))"
