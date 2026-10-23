#!/usr/bin/env bash
# Pathway-task rules (Addenda G, I, J/K of docs/preregistration_trees_xai.md), same procedure as DA/OX.
# Needs the pathway RF XAI (run_trees_xai.sh) and consensus (run_alternatives.sh) first.
set -u
cd "$(dirname "$0")/../src"
PY=../.venv/bin/python
run () { echo "=== $* ($(date +%H:%M))"; "$@" || echo "!!! FAILED: $*"; }
until grep -q ALTERNATIVES_DONE ../logs/run_alternatives.log 2>/dev/null; do sleep 120; done
run $PY 51_rulebook.py --task pw
run $PY 53_forest_rules.py --task pw
run $PY 54_rule_interpretation.py --task pw
run $PY 47_trees_xai.py cv --task pw --model dt3 --jobs 2
run $PY 50_trees_xai_report_data.py
run $PY 52_rulebook_report_data.py
run $PY 55_grounded_catalogue_data.py
run $PY 57_universal_data.py
echo "PW_RULES_DONE ($(date +%H:%M))"
