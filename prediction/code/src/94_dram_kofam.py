"""Addendum B (X1 comparator): DRAM in KOfam mode via its own distill step.

Step 1 (this env): KOfam annotation of the 100 primary new species for the 3 321 KOs referenced by DRAM's distill
forms, with the same KofamScan rules and thresholds as our models (src/predict.py ko_features_for). Writes a
DRAM-style annotations.tsv (one row per genome x KO present; DRAM's product uses presence only).
Step 2 (DRAM env): `DRAM.py distill` on that table; step 3 (this env): score with the Addendum B mapping.

usage: python src/94_dram_kofam.py annotate | score
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
D = ROOT / "data/additional/dram"
OUT = ROOT / "results/additional"

MAP = {
    "donor_sulfur_oxidation": ["Thiosulfate oxidation by SOX complex, thiosulfate => sulfate"],
    "donor_ammonia_oxidation": ["Bacterial/Archaeal ammonia oxidation"],
    "donor_nitrite_oxidation": ["nitrite=> nitrate"],
    "acceptor_nitrate": ["nitrate => nitrite"],
    "acceptor_sulfate": ["dissimilatory sulfate reduction (and oxidation) sulfate => sulfide"],
    "acceptor_sulfur": ["thiosulfate => sulfite", "tetrathionate => thiosulfate"],
    "acceptor_CO2_methanogen_acetogen": ["Key functional gene"],
}
# M00154, M00155, M00416, M00417 (low affinity) and M00156 (cbb3); bd oxidase (M00153) excluded
OXY_COLUMNS = ["Complex IV Low affinity:", "Complex IV High affinity: Cytochrome c oxidase, cbb3-type"]


def annotate():
    predict = __import__("predict")
    calls = pd.read_csv(OUT / "x1_calls.csv")
    accs = sorted(calls.accession.unique())
    kos = [k.strip() for k in (D / "dram_kos.txt").read_text().split()]
    tables, seqs, _ = predict.prepare_proteins(accs)
    X = predict.ko_features_for(tables, seqs, kos)
    X.to_parquet(D / "ko_presence_dram_kos.parquet")
    rows = []
    for g in X.index:
        present = [k for k in X.columns if X.at[g, k] > 0]
        for i, k in enumerate(present, 1):
            rows.append({"": f"{g}_{i}", "fasta": g, "scaffold": g, "gene_position": i, "ko_id": k,
                         "kegg_hit": k, "rank": "C"})
    pd.DataFrame(rows).to_csv(D / "annotations.tsv", sep="\t", index=False)
    print(X.shape, "KOs present per genome, median:", int((X > 0).sum(axis=1).median()))


def score():
    prod = pd.read_csv(D / "distill_out/product.tsv", sep="\t", index_col=0)
    cols = list(prod.columns)
    def find(name):
        hit = [c for c in cols if c.split(": ")[-1].strip().strip('"') == name.strip() or c.strip() == name.strip()]
        if not hit:
            hit = [c for c in cols if name.strip() in c]
        return hit
    calls = pd.read_csv(OUT / "x1_calls.csv")
    pred, used = {}, {}
    for cat, names in MAP.items():
        cs = sum((find(n) for n in names), [])
        used[cat] = cs
        pred[cat] = prod[cs].astype(str).apply(lambda r: r.str.lower().isin(["true", "1", "1.0"]).any(), axis=1) \
            if cs else pd.Series(False, index=prod.index)
    oxy = [c for c in cols if any(c.startswith(m) for m in OXY_COLUMNS)]
    used["acceptor_oxygen"] = oxy
    pred["acceptor_oxygen"] = (prod[oxy].apply(pd.to_numeric, errors="coerce") >= 0.5).any(axis=1) if oxy \
        else pd.Series(False, index=prod.index)
    P = pd.DataFrame(pred)
    calls["call_DRAM"] = [int(P.at[a, c]) if (a in P.index and c in P.columns) else 0
                          for a, c in zip(calls.accession, calls.category)]
    calls[["accession", "category", "call_DRAM"]].to_csv(OUT / "x1_dram_calls.csv", index=False)
    print({k: v for k, v in used.items()})


if __name__ == "__main__":
    {"annotate": annotate, "score": score}[sys.argv[1]]()
