"""Experiment H (docs/preregistration_H_hybrid_ranking.md): combined ranking = METABOLIC call (0/1) + released RedoxKO
score, compared with RedoxKO (HH1) and METABOLIC (HH2) by macro AUROC on the primary 100 new species and the 43 withheld
species (HH3). Paired bootstrap over description papers, 2,000 resamples, seed 20261017.
Output: results/additional/hybrid_ranking/results.json, per_category.csv
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/additional/hybrid_ranking"
OUT.mkdir(parents=True, exist_ok=True)
SEED, B = 20261017, 2000
COLS = {"combined": "score_H", "RedoxKO": "score_B2", "METABOLIC": "call_METABOLIC"}


def load(f):
    d = pd.read_csv(f)
    d["call_METABOLIC"] = d.call_METABOLIC.fillna(0).astype(float)
    d["score_H"] = d.call_METABOLIC + d.score_B2
    return d


def macro(d, col):
    v = [roc_auc_score(g.expected, g[col]) for _, g in d.groupby("category") if g.expected.nunique() == 2]
    return float(np.mean(v)) if v else np.nan


def evaluate(d):
    groups = {p: g for p, g in d.groupby("paper")}
    papers = np.array(list(groups))
    rng = np.random.default_rng(SEED)
    stat = lambda s: (macro(s, "score_H") - macro(s, "score_B2"), macro(s, "score_H") - macro(s, "call_METABOLIC"))
    pt = stat(d)
    bs = np.array([stat(pd.concat([groups[p] for p in rng.choice(papers, len(papers))])) for _ in range(B)])
    ci = lambda k: [float(x) for x in np.nanpercentile(bs[:, k], [2.5, 97.5])]
    r = {"n_tests": int(len(d)), "n_papers": int(len(papers)),
         **{f"macro_auroc_{k}": macro(d, c) for k, c in COLS.items()},
         **{f"pooled_auroc_{k}": float(roc_auc_score(d.expected, d[c])) for k, c in COLS.items()},
         "HH1_combined_minus_RedoxKO": pt[0], "HH1_ci": ci(0), "HH1_success": bool(ci(0)[0] > 0),
         "HH2_combined_minus_METABOLIC": pt[1], "HH2_ci": ci(1), "HH2_success": bool(ci(1)[0] > 0)}
    per = [dict(category=t, n=len(g), **{k: roc_auc_score(g.expected, g[c]) for k, c in COLS.items()})
           for t, g in d.groupby("category") if g.expected.nunique() == 2]
    return r, pd.DataFrame(per)


res, pers = {}, []
for name, f in [("primary100", "results/additional/x1_calls.csv"),
                ("withheld43", "results/additional/hybrid/holdout_calls.csv")]:
    r, p = evaluate(load(ROOT / f))
    res[name] = r
    pers.append(p.assign(set=name))
p, w = res["primary100"], res["withheld43"]
res["HH3_confirmation"] = {
    "HH1_same_sign": bool(np.sign(w["HH1_combined_minus_RedoxKO"]) == np.sign(p["HH1_combined_minus_RedoxKO"])),
    "HH2_same_sign": bool(np.sign(w["HH2_combined_minus_METABOLIC"]) == np.sign(p["HH2_combined_minus_METABOLIC"]))}
(OUT / "results.json").write_text(json.dumps(res, indent=2))
pd.concat(pers).to_csv(OUT / "per_category.csv", index=False)
print(json.dumps(res, indent=2, default=lambda x: round(float(x), 4)))
