#!/usr/bin/env bash
# Audit-fix rerun (docs/audit_rules_reports.md): novel-only external set, recall fix, strict nesting, direction-aware
# grounding, union pooling of unseen-lineage hits. DA and OX now; PW runs via run_pw_rules.sh when its models finish.
set -u
cd "$(dirname "$0")/../src"
PY=../.venv/bin/python
run () { echo "=== $* ($(date +%H:%M))"; "$@" || echo "!!! FAILED: $*"; }
for t in da ox; do
  run $PY 51_rulebook.py --task $t
  run $PY 53_forest_rules.py --task $t
  run $PY 54_rule_interpretation.py --task $t
done
for t in da ox pw; do
  run $PY 56_universal_rules.py --task $t
  run $PY 56_universal_rules.py --task $t --posthoc
done
run $PY 50_trees_xai_report_data.py
run $PY 52_rulebook_report_data.py
run $PY 55_grounded_catalogue_data.py
run $PY 57_universal_data.py
echo "RULES_AUDIT_RERUN_DONE ($(date +%H:%M))"
