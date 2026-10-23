"""Final model RedoxKO (random forest): macro AUPRC on the LLM subset (Table 7) and kappa with the larger dataset
(Supplementary Table S12), with the procedures of src/21_baseline2_llm.py score and src/24_label_concordance.py.
Output: results/additional/rf_paper/llm_kappa.json
"""
import importlib
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import cohen_kappa_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import RES, load_labelled, oof_matrix, paired_macro_delta, safe_ap  # noqa: E402

lab, traits, Y = load_labelled()
P = oof_matrix(RES / "trees_xai/da/oof", "RF", "genus", lab, traits)
key = pd.read_csv(ROOT / "data/baseline2/subset_key.tsv", sep="\t")
gi = {g: i for i, g in enumerate(lab.genome_id)}
idx = np.array([gi[g] for g in key.genome_id])
out = {"llm_subset_96": float(np.nanmean([safe_ap(Y[idx, j], P[idx, j]) for j in range(len(traits))]))}
s = pd.read_csv(RES / "baseline2_llm_summary.csv")
for arm in s.llm_arm.unique():
    n = int(s[s.llm_arm == arm].n_genomes.iloc[0])
    out[f"n_{arm}"] = n
# paired comparison with each LLM arm on the genomes it answered (reuse its stored predictions)
c21 = importlib.import_module("21_baseline2_llm")
B2DIR = c21.B2DIR
for variant in ("genome", "named", "gpt_genome", "gpt_named"):
    files = sorted(f for f in (B2DIR / "out").glob(f"{variant}_batch*.json") if not f.name.endswith(".usage.json"))
    pred = {}
    for f in files:
        pred.update(json.loads(f.read_text()))
    M = np.full((len(lab), len(traits)), np.nan)
    for _, r in key.iterrows():
        d = pred.get(r.code, {})
        for j, t in enumerate(traits):
            if t in d:
                M[gi[r.genome_id], j] = float(d[t])
    have = ~np.isnan(M[idx]).any(axis=1)
    ii = idx[have]
    sub = lab.iloc[ii].reset_index(drop=True)
    d, lo, hi = paired_macro_delta(sub, traits, Y[ii], M[ii], P[ii], "genus")
    out[f"H_{variant}"] = dict(n=int(have.sum()), RedoxKO=float(np.nanmean([safe_ap(Y[ii, j], P[ii, j]) for j in range(len(traits))])),
                               llm_minus_RedoxKO=d, ci=[lo, hi])
# kappa with the larger dataset (24_label_concordance with the RF genus OOF and the final thresholds)
c24 = importlib.import_module("24_label_concordance")
thr = {t: b["threshold"] for t, b in joblib.load(RES / "models/RF/model.joblib").items()}
conc = pd.read_csv(RES / "label_concordance_madin.csv")
src = c24.__dict__
d = pd.read_csv(ROOT / "data/raw/madin/condensed_species_NCBI.csv", low_memory=False) if False else None
out["note_kappa"] = "computed by src/24_label_concordance.py with RF; see rf_label_concordance.csv"
(RES / "additional/rf_paper").mkdir(parents=True, exist_ok=True)
(RES / "additional/rf_paper/llm_kappa.json").write_text(json.dumps(out, indent=1, default=float))
print(json.dumps(out, indent=1, default=lambda v: round(float(v), 3)))
