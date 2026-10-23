"""X1 comparator: BacDive-AI aerobic model (Koblitz et al. 2025, github.com/JKoblitz/bacdive-AI @a1bef3e) on the
primary 100 newly described species.

The tool expects InterProScan Pfam output and uses the set of Pfam accessions with a match E-value <= 1e-20.
InterProScan is not installable here (deviation D2), so matches come from pyhmmer hmmsearch of Pfam-A 38.2 with
gathering (GA) cut-offs, restricted to the 1 586 Pfam families the model uses, keeping a family when any domain
passes GA with independent E-value <= 1e-20. The model is applied exactly as the tool's predict.py does.

usage:
  .venv_sk14/bin/python src/87b_bacdive_ai.py features   # writes the model's Pfam list (scikit-learn 1.4 env)
  .venv/bin/python      src/87b_bacdive_ai.py scan       # pyhmmer scan (main env)
  .venv_sk14/bin/python src/87b_bacdive_ai.py score      # predictions (scikit-learn 1.4 env)
"""
import json
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BD = ROOT / "data/additional/bacdive_ai"
OUT = ROOT / "results/additional"
EVALUE = 1e-20


def features():
    d = pickle.load(open(BD / "models/aerobic_data.p", "rb"))
    (BD / "aerobic_pfams.json").write_text(json.dumps(list(d["categories"])))
    print(len(d["categories"]), "Pfam families")


def scan():
    import pandas as pd
    import pyhmmer
    want = set(json.loads((BD / "aerobic_pfams.json").read_text()))
    hmms = []
    with pyhmmer.plan7.HMMFile(str(BD / "Pfam-A.hmm")) as f:
        for h in f:
            acc = (h.accession.decode() if isinstance(h.accession, bytes) else h.accession).split(".")[0]
            if acc in want:
                hmms.append(h)
    print(len(hmms), "of", len(want), "profiles found in Pfam-A 38.2", flush=True)
    calls = pd.read_csv(ROOT / "results/additional/x1_calls.csv")
    accs = sorted(calls.accession.unique())
    alphabet = pyhmmer.easel.Alphabet.amino()
    rows = []
    for a in accs:
        faa = ROOT / "data/predict_cache/genomes" / a / "protein.faa"
        if not faa.exists():
            faa = ROOT / "data/raw/genomes" / a / "protein.faa"
        with pyhmmer.easel.SequenceFile(str(faa), digital=True, alphabet=alphabet) as sf:
            seqs = sf.read_block()
        found = set()
        for top in pyhmmer.hmmsearch(hmms, seqs, bit_cutoffs="gathering", cpus=8):
            q = top.query.accession
            q = (q.decode() if isinstance(q, bytes) else q).split(".")[0]
            for hit in top:
                if hit.included and any(dom.included and dom.i_evalue <= EVALUE for dom in hit.domains):
                    found.add(q)
                    break
        rows.append({"accession": a, "pfams": ";".join(sorted(found)), "n": len(found)})
        print(a, len(found), flush=True)
    pd.DataFrame(rows).to_csv(BD / "pfam_presence.csv", index=False)


def score():
    import csv
    import numpy as np
    d = pickle.load(open(BD / "models/aerobic_data.p", "rb"))
    clf, cats = d["model"], d["categories"]
    out = []
    with open(BD / "pfam_presence.csv") as f:
        for r in csv.DictReader(f):
            pf = set(r["pfams"].split(";")) if r["pfams"] else set()
            X = [1 if k in pf else 0 for k in cats]
            proba = clf.predict_proba([X])[0]
            out.append((r["accession"], float(proba[list(clf.classes_).index(1)]), int(r["n"])))
    with open(OUT / "x1_bacdive_ai_aerobic.csv", "w") as f:
        f.write("accession,p_aerobic,n_model_pfams_present\n")
        for a, p, n in out:
            f.write(f"{a},{p:.4f},{n}\n")
    print(len(out), "genomes scored; mean p_aerobic", round(float(np.mean([o[1] for o in out])), 3))


if __name__ == "__main__":
    {"features": features, "scan": scan, "score": score}[sys.argv[1]]()
