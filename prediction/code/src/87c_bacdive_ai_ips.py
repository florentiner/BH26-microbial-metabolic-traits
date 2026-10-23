"""BacDive-AI aerobic model on genuine InterProScan 6.0.2.2 Pfam output (removes deviation D2's approximation).

Per genome, the tool's own logic (predict.py of github.com/JKoblitz/bacdive-AI @a1bef3e): Pfam accessions from
column 5 of the InterProScan TSV with E-value (column 9) <= 1e-20, as a set; feature vector over the model's
categories; probability of class 1. Run in the scikit-learn 1.4.0 environment.

usage: data/additional/bacdive_ai/.venv_sk14/bin/python src/87c_bacdive_ai_ips.py
"""
import csv
import pickle
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TSV = ROOT / "data/additional/interproscan/full_out/oxygen49.tsv"
EVALUE = 1e-20

pf = defaultdict(set)
with open(TSV) as f:
    for line in csv.reader(f, delimiter="\t"):
        if line[3] != "Pfam" or float(line[8]) > EVALUE:
            continue
        pf[line[0].split("__")[0]].add(line[4])
d = pickle.load(open(ROOT / "data/additional/bacdive_ai/models/aerobic_data.p", "rb"))
clf, cats = d["model"], d["categories"]
out = ROOT / "results/additional/x1_bacdive_ai_aerobic_interproscan.csv"
with open(out, "w") as f:
    f.write("accession,p_aerobic,n_model_pfams_present\n")
    for acc in sorted(pf):
        X = [1 if k in pf[acc] else 0 for k in cats]
        p = clf.predict_proba([X])[0][list(clf.classes_).index(1)]
        f.write(f"{acc},{p:.4f},{sum(X)}\n")
print(len(pf), "genomes scored ->", out)
