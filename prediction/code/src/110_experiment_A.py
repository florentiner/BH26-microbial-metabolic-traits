"""Experiment A (docs/preregistration_Bclean_A.md, test HA): on the 43 withheld new species, do the released KO model's
scores rank true above false positives among METABOLIC's positive calls? Paper-cluster bootstrap, 3,000 resamples,
seed 20261017. Success: lower 95% bound above 0.5.
Output: results/additional/bclean/experiment_A.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "results/additional/hybrid/holdout_calls.csv")
m = d[d.call_METABOLIC == 1].copy()
groups = {p: g for p, g in m.groupby("paper")}
papers = np.array(list(groups))

point = roc_auc_score(m.expected, m.score_B2)
rng = np.random.default_rng(20261017)
vals = []
for _ in range(3000):
    s = pd.concat([groups[p] for p in rng.choice(papers, len(papers))])
    if s.expected.nunique() == 2:
        vals.append(roc_auc_score(s.expected, s.score_B2))
lo, hi = np.percentile(vals, [2.5, 97.5])
out = pd.DataFrame([dict(test="HA", n_tests=len(m), n_true=int(m.expected.sum()), n_false=int((m.expected == 0).sum()),
                         n_papers=len(papers), auroc=point, ci_low=lo, ci_high=hi, n_valid_resamples=len(vals),
                         success=bool(lo > 0.5))])
(ROOT / "results/additional/bclean").mkdir(parents=True, exist_ok=True)
out.to_csv(ROOT / "results/additional/bclean/experiment_A.csv", index=False)
print(out.round(3).to_string(index=False))
