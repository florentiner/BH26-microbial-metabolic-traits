#!/usr/bin/env bash
# Full pathway-task pipeline, in the pre-registered order (docs/preregistration_pathways.md).
# Cheap arms and the audit run first, so a problem surfaces before the expensive ones start.
set -u
cd "$(dirname "$0")/../src"
PY=../.venv/bin/python
LOG=../logs
JOBS="${JOBS:-6}"

run () { echo "=== $* ($(date +%H:%M))"; "$@" || echo "!!! FAILED: $*"; }

run $PY 37_pathway_audit.py                      # pre-arm integrity checks (hash, fold blocking, feature purity)
run $PY 32_pathway_arms.py nulls
run $PY 32_pathway_arms.py PA
run $PY 32_pathway_arms.py PB1 --jobs "$JOBS"    # 852-KO energy panel
run $PY 32_pathway_arms.py PB2 --jobs "$JOBS"    # full KO profile (pre-registered primary)
run $PY 33_pathway_metrics.py                    # first read-out: PA vs PB1 vs PB2 vs nulls
run $PY 32_pathway_arms.py PB3 --jobs "$JOBS"    # + module completeness
run $PY 32_pathway_arms.py PB4 --jobs 1          # multi-label MLP (hypothesis P6; torch, single process)
run $PY 33_pathway_metrics.py                    # final metrics over all arms
run $PY 33_pathway_metrics.py --subset no_faprotax
run $PY 33_pathway_metrics.py --subset informative_negatives
run $PY 37_pathway_audit.py                      # re-audit, now including the out-of-fold prediction checks
run $PY 35_pathway_train_final.py                # refit + release gate
run $PY 34_pathway_mag_robustness.py --jobs "$JOBS"
run $PY 36_pathway_crosscheck.py
run $PY 39_task_comparison.py                    # PP5: the two tasks on a prevalence-free footing
run $PY 38_pathway_results_md.py                 # docs/RESULTS_pathways.md, every number read from the CSVs
echo "PATHWAYS_DONE ($(date +%H:%M))"
