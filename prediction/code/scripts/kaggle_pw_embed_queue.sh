#!/usr/bin/env bash
# Push pathway embedding parts 3 and 4 when a Kaggle GPU slot frees up (max 2 concurrent sessions); log every status
# change. Uses only the approved account (token file passed as $1).
set -u
TOK="$1"; KAG="$2"; DIR="$(dirname "$0")/../kaggle/pathway_embed"
export KAGGLE_API_TOKEN=$(awk '{print $2}' "$TOK")
queue=(3 4); declare -A last
st () { "$KAG" kernels status "${KAGGLE_USER:-your-kaggle-user}/redox-pw-embed-$1" 2>/dev/null | grep -v -i warn | tail -1 | sed 's/.*status //'; }
while true; do
  running=0; done_all=1
  for k in 1 2 3 4; do
    s=$(st $k); [ -z "$s" ] && s="not pushed"
    if [ "${last[$k]:-}" != "$s" ]; then echo "$(date +%H:%M) part$k $s"; last[$k]="$s"; fi
    case "$s" in *RUNNING*|*QUEUED*) running=$((running+1)); done_all=0;; *COMPLETE*|*ERROR*|*CANCEL*) ;; *) done_all=0;; esac
  done
  while [ ${#queue[@]} -gt 0 ] && [ $running -lt 2 ]; do
    k=${queue[0]}
    out=$("$KAG" kernels push -p "$DIR/part$k" 2>&1 | grep -v -i warn | tail -1)
    echo "$(date +%H:%M) push part$k: $out"
    if echo "$out" | grep -q "successfully"; then queue=("${queue[@]:1}"); running=$((running+1)); else break; fi
  done
  if [ ${#queue[@]} -eq 0 ] && [ $done_all -eq 1 ]; then echo "$(date +%H:%M) ALL_PARTS_FINISHED"; break; fi
  sleep 600
done
