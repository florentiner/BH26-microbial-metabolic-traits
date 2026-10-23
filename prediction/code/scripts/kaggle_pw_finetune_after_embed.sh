#!/usr/bin/env bash
# Push the pathway fine-tuning kernel once all four pathway embedding kernels have finished successfully.
set -u
TOK="$1"; KAG="$2"; ROOT="$(dirname "$0")/.."
export KAGGLE_API_TOKEN=$(awk '{print $2}' "$TOK")
until grep -q ALL_PARTS_FINISHED "$ROOT/logs/kaggle_pw_embed.log" 2>/dev/null; do sleep 300; done
ok=1
for k in 1 2 3 4; do
  s=$("$KAG" kernels status "${KAGGLE_USER:-your-kaggle-user}/redox-pw-embed-$k" 2>/dev/null | grep -v -i warn | tail -1)
  echo "$(date +%H:%M) part$k: $s"
  echo "$s" | grep -q COMPLETE || ok=0
done
if [ $ok -eq 1 ]; then
  echo "$(date +%H:%M) push finetune: $("$KAG" kernels push -p "$ROOT/kaggle/pathway_finetune" 2>&1 | grep -v -i warn | tail -1)"
else
  echo "$(date +%H:%M) FINETUNE_NOT_PUSHED: an embedding part did not complete"
fi
prev=""
while true; do
  s=$("$KAG" kernels status ${KAGGLE_USER:-your-kaggle-user}/redox-pw-finetune 2>/dev/null | grep -v -i warn | tail -1)
  [ "$s" != "$prev" ] && echo "$(date +%H:%M) finetune: $s" && prev="$s"
  case "$s" in *COMPLETE*|*ERROR*|*CANCEL*) echo "FINETUNE_DONE"; break;; esac
  [ $ok -eq 0 ] && break
  sleep 600
done
