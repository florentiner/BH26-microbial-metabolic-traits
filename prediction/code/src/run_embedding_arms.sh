#!/bin/sh
# Wait until every labelled-genome protein has an ESM-2 embedding, then build genome features (labelled genomes)
# and run arms C, E, E_MAG + metrics + ledger. usage: sh src/run_embedding_arms.sh
cd "$(dirname "$0")"
PY=../.venv/bin/python
until $PY - <<'PYEOF'
import glob, sys
import numpy as np, pandas as pd
lab = pd.read_csv("../data/processed/labels.tsv", sep="\t")
need = set()
for a in lab.loc[lab.labelled, "genome_id"]:
    need.update(pd.read_parquet(f"../data/processed/proteins/{a}.parquet", columns=["seq_hash"]).seq_hash)
done = set()
for f in glob.glob("../data/processed/esm2_t12_35M/shard_*.npz"):
    done.update(np.load(f)["hashes"].tolist())
print("labelled remaining", len(need - done), flush=True)
sys.exit(0 if need <= done else 1)
PYEOF
do sleep 300; done
$PY 07_genome_embeddings.py --labelled-only || exit 1
for a in C E E_MAG; do $PY 06_run_arms.py $a --jobs 6 || exit 1; done
$PY 15_research_ledger.py --reference B2 --variants C E E_MAG
echo EMB_ARMS_DONE
