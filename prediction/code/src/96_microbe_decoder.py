"""Addendum C (X1 comparator): Microbe Decoder machine-learning tool via its public API, from KO content only.

Full KOfam profiles (10 517 KOs) of the 100 primary new species: 92 reused from the pathway test
(validation_new_organisms/pathway100_ko_matrix.parquet, same release and thresholds), the rest scanned here.
Each genome's present KOs are sent as `gene_functions` (no organism name) to POST /compute/ml.

usage: python src/96_microbe_decoder.py kos | query | score
"""
import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
D = ROOT / "data/additional/microbe_decoder"
OUT = ROOT / "results/additional"
API = "https://api.microbe-decoder.org/compute/ml"
MODELS = ["Anaerobe (oxygen tolerance)", "Nitrate (metabolites utilized)", "Methanogenesis (type of metabolism)"]


def kos():
    D.mkdir(parents=True, exist_ok=True)
    calls = pd.read_csv(OUT / "x1_calls.csv")
    accs = sorted(calls.accession.unique())
    M = pd.read_parquet(ROOT / "validation_new_organisms/pathway100_ko_matrix.parquet")
    have = [a for a in accs if a in M.index]
    todo = [a for a in accs if a not in M.index]
    parts = [M.loc[have]]
    if todo:
        predict = __import__("predict")
        tables, seqs, _ = predict.prepare_proteins(todo)
        parts.append(predict.ko_features_for(tables, seqs, list(M.columns)))
    X = pd.concat(parts).reindex(columns=M.columns).fillna(0)
    X.to_parquet(D / "ko_full_100.parquet")
    print(X.shape, "reused", len(have), "scanned", len(todo), "median KOs", int((X > 0).sum(axis=1).median()))


def post(payload, tries=5):
    data = json.dumps(payload).encode()
    for i in range(tries):
        try:
            req = urllib.request.Request(API, data=data, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=600) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(10 * (i + 1))
    raise RuntimeError(err)


def query(only=None):
    X = pd.read_parquet(D / "ko_full_100.parquet")
    out_f = D / "predictions.jsonl"
    done = set()
    if out_f.exists():
        done = {json.loads(l)["accession"] for l in out_f.read_text().splitlines()}
    accs = [a for a in X.index if a not in done][: only or None]
    for n, a in enumerate(accs, 1):
        present = [k for k in X.columns if X.at[a, k] > 0]
        gf = {"gene_oid": list(range(1, len(present) + 1)), "Genome ID": [n] * len(present),
              "Genome Name": [a] * len(present), "KO": [f"KO:{k}" for k in present]}
        res = post({"gene_functions": gf, "model_names": MODELS})
        with open(out_f, "a") as f:
            f.write(json.dumps({"accession": a, "n_kos": len(present), "response": res}) + "\n")
        print(a, len(present), json.dumps(res)[:200], flush=True)
        time.sleep(1)


def score():
    rows = []
    for l in (D / "predictions.jsonl").read_text().splitlines():
        r = json.loads(l)
        p = {d["Model"]: d["Probability"] for d in r["response"]["probabilities"]}
        rows.append({"accession": r["accession"], **p})
    P = pd.DataFrame(rows).set_index("accession")
    calls = pd.read_csv(OUT / "x1_calls.csv")
    s, c = [], []
    for a, cat in zip(calls.accession, calls.category):
        if cat == "acceptor_oxygen":
            v = 1 - P.at[a, MODELS[0]]
        elif cat == "acceptor_nitrate":
            v = P.at[a, MODELS[1]]
        elif cat == "acceptor_CO2_methanogen_acetogen":
            v = P.at[a, MODELS[2]]
        else:
            v = np.nan
        s.append(v)
        c.append(int(v >= 0.5) if pd.notna(v) else 0)
    calls["score_MicrobeDecoder"], calls["call_MicrobeDecoder"] = s, c
    calls[["accession", "category", "score_MicrobeDecoder", "call_MicrobeDecoder"]].to_csv(
        OUT / "x1_microbe_decoder_calls.csv", index=False)
    print(P.describe().round(2).to_string())


if __name__ == "__main__":
    a = sys.argv[1]
    if a == "query":
        query(int(sys.argv[2]) if len(sys.argv) > 2 else None)
    else:
        {"kos": kos, "score": score}[a]()
