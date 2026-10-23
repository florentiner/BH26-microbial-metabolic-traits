"""95% cluster-bootstrap intervals for every row of Table 1 (donor/acceptor cross-validation), at genus and family
level, from the saved out-of-fold predictions. Same resampling as the existing intervals: whole genera (or families)
resampled with replacement, 1,000 resamples, seed 20260917, macro AUPRC over the 14 categories.
Output: results/table1_cis.csv
"""
import numpy as np
import pandas as pd

from pathway_common import RES
from redox_common import load_labelled, oof_matrix, safe_ap

ROWS = [("Frequency baseline", RES / "oof", "null_prevalence"), ("Curated marker rules", RES / "oof", "A"),
        ("Taxonomic baseline", RES / "oof", "null_taxNN"), ("ESM-2", RES / "oof", "C"),
        ("Attention-based MIL", RES / "oof", "D_nopf"), ("Bacformer-large", RES / "oof", "E_large"),
        ("RedoxKO", RES / "oof", "B2"), ("RedoxKO-ens", RES / "oof", "B2XL"),
        ("Random forest", RES / "trees_xai/da/oof", "RF"), ("RedoxKO-ens + ESM-2", RES / "oof", "B2XLC"),
        ("RedoxKO-ens-BF", RES / "oof", "B2XLE_large")]

lab, traits, Y = load_labelled()
obs = ~np.isnan(Y)


def macro(P, idx):
    return np.nanmean([safe_ap(Y[idx, j][obs[idx, j]], P[idx, j][obs[idx, j]]) for j in range(len(traits))])


rows = []
for level in ("genus", "family"):
    groups = lab[f"group_{level}"].to_numpy()
    ug = np.unique(groups)
    members = {g: np.where(groups == g)[0] for g in ug}
    rng = np.random.default_rng(20260917)
    boots = [np.concatenate([members[g] for g in rng.choice(ug, len(ug))]) for _ in range(1000)]
    allidx = np.arange(len(lab))
    for name, d, arm in ROWS:
        P = oof_matrix(d, arm, level, lab, traits)
        b = np.array([macro(P, i) for i in boots])
        rows.append(dict(model=name, level=level, auprc=macro(P, allidx), lo=np.percentile(b, 2.5), hi=np.percentile(b, 97.5)))
out = pd.DataFrame(rows)
out.to_csv(RES / "table1_cis.csv", index=False, float_format="%.4f")
print(out.pivot(index="model", columns="level", values=["auprc", "lo", "hi"]).round(3).to_string())
