"""B-clean (docs/preregistration_Bclean_A.md): 852-KO panel features for the added genomes, computed with the
training pipeline (NCBI proteins, pyhmmer, KOfam release of 25 August 2026) through src/predict.py.
Resumable: genomes are processed in chunks, each chunk saved; failed accessions are logged and skipped.
Output: results/additional/bclean/features/chunk_*.parquet, results/additional/bclean/failed.tsv
"""
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import predict  # noqa: E402

OUT = ROOT / "results/additional/bclean"
FEAT = OUT / "features"
FEAT.mkdir(parents=True, exist_ok=True)
CHUNK = 20

accs = pd.read_csv(OUT / "added_genomes.tsv", sep="\t").genome_id.tolist()
done = set()
for f in FEAT.glob("chunk_*.parquet"):
    done |= set(pd.read_parquet(f).index)
failed_f = OUT / "failed.tsv"
failed = set(pd.read_csv(failed_f, sep="\t").genome_id) if failed_f.exists() else set()
todo = [a for a in accs if a not in done and a not in failed]
print(f"{len(accs)} genomes, {len(done)} done, {len(failed)} failed, {len(todo)} to do", flush=True)

t0 = time.time()
for i in range(0, len(todo), CHUNK):
    part = todo[i:i + CHUNK]
    tables, seqs, ok = {}, {}, []
    for acc in part:
        try:
            t, s, _ = predict.prepare_proteins([acc])
            tables.update(t)
            seqs.update(s)
            ok.append(acc)
        except (SystemExit, Exception) as e:  # no NCBI annotation, download error
            with open(failed_f, "a") as fh:
                if fh.tell() == 0:
                    fh.write("genome_id\treason\n")
                fh.write(f"{acc}\t{str(e)[:150]}\n")
    if ok:
        X = predict.ko_features(tables, seqs)
        X.to_parquet(FEAT / f"chunk_{len(list(FEAT.glob('chunk_*.parquet'))):03d}.parquet")
    n = min(i + CHUNK, len(todo))
    rate = (time.time() - t0) / n
    print(f"{n}/{len(todo)} genomes, {rate:.1f} s per genome, about {rate * (len(todo) - n) / 60:.0f} min left",
          flush=True)
print("done", flush=True)
