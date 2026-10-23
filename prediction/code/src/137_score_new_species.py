"""KO panel and RedoxKO scores for further new species (pre-registration AT), with the released workflow's functions.

For each accession: NCBI proteins (cached in data/predict_cache), the 852-KO panel with the KofamScan rules of
training, and the final RedoxKO model (results/models/B2) with its stored thresholds. Appends to
validation_new_organisms/ko_panel_auto.parquet and validation_new_organisms/scores_auto_B2.csv.
usage: python src/137_score_new_species.py GCF_... [GCF_ ...]  |  --file accessions.txt
"""
import argparse
import sys
from pathlib import Path

import joblib
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
pr = __import__("predict")
V = ROOT / "validation_new_organisms"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accessions", nargs="*")
    ap.add_argument("--file")
    a = ap.parse_args()
    accs = list(a.accessions) + ([l.strip() for l in open(a.file) if l.strip()] if a.file else [])
    kp, sc = V / "ko_panel_auto.parquet", V / "scores_auto_B2.csv"
    done = set(pd.read_parquet(kp).index) if kp.exists() else set()
    accs = [x for x in dict.fromkeys(accs) if x not in done]
    if not accs:
        print("nothing to do")
        return
    tables, seqs, _ = pr.prepare_proteins(accs)
    X = pr._isolated("ko_features", tables, seqs)
    bundle = joblib.load(ROOT / "results/models/B2/model.joblib")
    rows = []
    for t, b in bundle.items():
        p = b["model"].predict_proba(X)
        for acc, s in zip(X.index, p):
            rows.append(dict(accession=acc, category=t, score=float(s), threshold=float(b["threshold"]),
                             call=int(s >= b["threshold"])))
    pd.concat([pd.read_parquet(kp), X]).to_parquet(kp) if kp.exists() else X.to_parquet(kp)
    new = pd.DataFrame(rows)
    (pd.concat([pd.read_csv(sc), new]) if sc.exists() else new).to_csv(sc, index=False, float_format="%.6f")
    print(f"scored {X.shape[0]} genomes")


if __name__ == "__main__":
    main()
