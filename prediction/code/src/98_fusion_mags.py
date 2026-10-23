"""R1: the TogoMCP donor/acceptor model B2XLE_large (KO panel + Bacformer-large) on the MAG-isolate pairs (X2) and the
habitat MAGs (X4). Bacformer-large vectors come from the Kaggle kernels kaggle/mag_embed (Pyrodigal proteins, verified
identical to src/predict.py --fasta); the KO panel is computed here by src/predict.py from the same FASTA.
Analyses reuse src/88_real_mags.py and src/90_habitat_case.py with outputs in results/additional/fusion/.

usage: python src/98_fusion_mags.py merge | predict | analyze
"""
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
EMB = ROOT / "data/additional/mag_embed"
FNA = ROOT / "data/additional/mags/fna"
PRED = ROOT / "data/additional/mags/pred_fusion"
OUT = ROOT / "results/additional/fusion"


def merge():
    parts = [pd.read_parquet(p) for p in sorted(EMB.glob("k*/genome_features_bacformer_large_mean.parquet"))]
    E = pd.concat(parts)
    E = E[~E.index.duplicated()]
    E.to_parquet(EMB / "bacformer_large_mean_all.parquet")
    print(E.shape, "genomes with a Bacformer-large vector")


def predict(batch=20):
    PRED.mkdir(parents=True, exist_ok=True)
    E = pd.read_parquet(EMB / "bacformer_large_mean_all.parquet")
    accs = sorted(a for a in E.index if (FNA / f"{a}.fna").exists())
    done = {p.stem[2:] for p in PRED.glob("b_*.csv")}
    todo = [a for a in accs if a not in done]
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        out = PRED / f"batch_{i:05d}_{chunk[0]}.csv"
        cmd = [sys.executable, str(ROOT / "src/predict.py"), "--arm", "B2XLE_large", "--pathways", "none",
               "--embeddings", str(EMB / "bacformer_large_mean_all.parquet"), "--out", str(out)]
        for a in chunk:
            cmd += ["--fasta", str(FNA / f"{a}.fna")]
        t0 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print("FAILED", chunk[0], r.stderr[-800:], flush=True)
            continue
        for a in chunk:
            (PRED / f"b_{a}.csv").write_text("done\n")
        print(f"{i + len(chunk)}/{len(todo)} in {time.time() - t0:.0f}s", flush=True)


def analyze():
    OUT.mkdir(parents=True, exist_ok=True)
    m88 = __import__("88_real_mags")
    m88.PRED, m88.OUT = PRED, OUT
    m88.analyze()
    m90 = __import__("90_habitat_case")
    m90.PRED, m90.OUT = PRED, OUT
    m90.analyze()


if __name__ == "__main__":
    {"merge": merge, "predict": predict, "analyze": analyze}[sys.argv[1]]()
