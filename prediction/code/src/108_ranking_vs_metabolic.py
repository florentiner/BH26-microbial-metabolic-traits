"""Exploratory (not pre-registered): do KO-model scores add information beyond METABOLIC's yes/no calls on the
100 newly described species? Bootstrap over description papers, as for all new-species results.
Output: results/additional/x6_ranking_vs_metabolic.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "results/additional/x1_calls.csv")
d = d[d.call_METABOLIC.notna()].copy()
groups = {p: g for p, g in d.groupby("paper")}
papers = np.array(list(groups))


def stats(s):
    pos, neg = s[s.call_METABOLIC == 1], s[s.call_METABOLIC == 0]
    q = s.score_B2.quantile(0.75)
    top = s[s.score_B2 >= q].expected.mean()
    met = pos.expected.mean()
    return {
        "auroc_within_metabolic_positive": roc_auc_score(pos.expected, pos.score_B2),
        "auroc_within_metabolic_negative": roc_auc_score(neg.expected, neg.score_B2),
        "precision_ko_top_quartile": top,
        "precision_metabolic_positive_calls": met,
        "precision_difference": top - met,
    }


point = stats(d)
rng = np.random.default_rng(1)
boot = []
for _ in range(3000):
    s = pd.concat([groups[p] for p in rng.choice(papers, len(papers))])
    try:
        boot.append(stats(s))
    except ValueError:  # a resample without both classes in a subset
        continue
b = pd.DataFrame(boot)
out = pd.DataFrame({"estimate": point, "ci_low": b.quantile(0.025), "ci_high": b.quantile(0.975)})
out["n_tests"] = len(d)
out.to_csv(ROOT / "results/additional/x6_ranking_vs_metabolic.csv")
print(out.round(3))
