"""Experiment F (docs/preregistration_F_full_ko_profile.md): full KO profile (the 10,517 KOs of the released pathway
model EnergyKO-M) for the 358 donor/acceptor training genomes and the new species not already in
validation_new_organisms/pathway100_ko_matrix.parquet, with src/predict.py ko_features_for (same KOfam release and
rules). Resumable, chunks of 20 genomes.
Output: results/additional/full_profile/features/chunk_*.parquet, failed.tsv, run.log
"""
import json
import sys
import time
from pathlib import Path

import joblib
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import predict  # noqa: E402
from redox_common import load_labelled  # noqa: E402

OUT = ROOT / "results/additional/full_profile"
FEAT = OUT / "features"
FEAT.mkdir(parents=True, exist_ok=True)
CHUNK = 20

rel = json.loads((ROOT / "results/models/pathways/RELEASED.json").read_text())
bundle = joblib.load(ROOT / "results/models/pathways" / rel["arm"] / "model.joblib")
KOS = sorted({c for b in bundle.values() for c in b["model"].columns if c.startswith("K")})
(OUT / "ko_list.txt").write_text("\n".join(KOS))
lab, _, _ = load_labelled()
new = set(pd.read_csv(ROOT / "results/additional/x1_calls.csv").accession) | \
      set(pd.read_csv(ROOT / "results/additional/hybrid/holdout_calls.csv").accession)
have = set(pd.read_parquet(ROOT / "validation_new_organisms/pathway100_ko_matrix.parquet").index)
accs = list(lab.genome_id) + sorted(new - have)
done = set()
for f in FEAT.glob("chunk_*.parquet"):
    done |= set(pd.read_parquet(f).index)
failed_f = OUT / "failed.tsv"
failed = set(pd.read_csv(failed_f, sep="\t").genome_id) if failed_f.exists() else set()
todo = [a for a in accs if a not in done and a not in failed]
print(f"KOs={len(KOS)}; {len(accs)} genomes ({len(lab)} training, {len(new - have)} new species), "
      f"{len(done)} done, {len(todo)} to do", flush=True)
t0, n0 = time.time(), 0
for i in range(0, len(todo), CHUNK):
    part = todo[i:i + CHUNK]
    tables, seqs, bad = {}, {}, []
    for acc in part:
        try:
            t, s, _ = predict.prepare_proteins([acc])
            tables.update(t)
            seqs.update(s)
        except (SystemExit, Exception) as e:
            bad.append((acc, str(e)[:200]))
    if bad:
        pd.DataFrame(bad, columns=["genome_id", "error"]).to_csv(failed_f, sep="\t", index=False,
                                                                mode="a", header=not failed_f.exists())
    if tables:
        X = predict.ko_features_for(tables, seqs, KOS).astype("int8")
        X.to_parquet(FEAT / f"chunk_{len(list(FEAT.glob('chunk_*.parquet'))):03d}.parquet")
    n0 += len(part)
    rate = (time.time() - t0) / n0
    print(f"{len(done) + n0}/{len(accs)} genomes, {rate:.1f} s per genome, about {rate * (len(todo) - n0) / 60:.0f} "
          f"min left", flush=True)
print("done", flush=True)
