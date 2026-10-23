"""95% cluster-bootstrap intervals (over description papers) for the macro AUROC of every model in Table 3
(100 newly described species, 255 growth tests). Macro = mean over categories with at least one tested positive and
one tested negative, as in src/82_novel_model_comparison.py. 1,000 resamples, seed 20260917.
Output: results/table3_auroc_cis.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "validation_new_organisms/novel100_model_scores.csv")
# cluster = description paper (PMID), taken from the per-call table; the scores file only keeps the species name
import re
_m = pd.read_csv(ROOT / "validation_new_organisms/mismatch178_B2.csv").drop_duplicates("accession").set_index("accession")
d["paper"] = d.accession.map(_m["paper"].map(lambda s: (re.search(r"PMID (\d+)", str(s)) or re.search(r"(.*)", str(s))).group(1)))
MODELS = {"B2XLE_large": "RedoxKO-ens-BF", "B2": "RedoxKO", "B2XL": "RedoxKO-ens", "E_large": "Bacformer-large alone",
          "C_large": "ESM-C alone"}


def macro_auroc(x, col):
    vals = []
    for _, g in x.groupby("category"):
        y = g.expected.astype(int)
        if 0 < y.sum() < len(y):
            vals.append(roc_auc_score(y, g[col]))
    return float(np.mean(vals))


papers = d.paper.unique()
idx = {p: d.index[d.paper == p].to_numpy() for p in papers}
rng = np.random.default_rng(20260917)
boots = [np.concatenate([idx[p] for p in rng.choice(papers, len(papers))]) for _ in range(1000)]
rows = []
for m, name in MODELS.items():
    col = f"score_{m}"
    b = np.array([macro_auroc(d.loc[i].reset_index(drop=True), col) for i in boots])
    rows.append(dict(model=name, auroc=macro_auroc(d, col), lo=np.percentile(b, 2.5), hi=np.percentile(b, 97.5)))
out = pd.DataFrame(rows)
out.to_csv(ROOT / "results/table3_auroc_cis.csv", index=False, float_format="%.4f")
print(out.round(3).to_string(index=False))
