"""Experiment V (docs/preregistration_V_validity_checks.md): V1 RedoxKO vs METABOLIC on external labels without
FAPROTAX; V2 accuracy on new species against KO-profile distance to the nearest training genome.
Output: results/additional/validity/{v1.json, v2.json, v2_distances.csv}"""
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from metabolic_calls import calls_from_rows  # noqa: E402
from redox_common import load_labelled  # noqa: E402

OUT = ROOT / "results/additional/validity"
OUT.mkdir(parents=True, exist_ok=True)
SEED, B = 20261018, 2000
thr = {t: b["threshold"] for t, b in joblib.load(ROOT / "results/models/B2/model.joblib").items()}
ci = lambda a: [float(x) for x in np.nanpercentile(a, [2.5, 97.5])]


def boot(d, cl, fn):
    g = {k: x for k, x in d.groupby(cl)}
    keys = np.array(list(g))
    rng = np.random.default_rng(SEED)
    return np.array([fn(pd.concat([g[k] for k in rng.choice(keys, len(keys))])) for _ in range(B)])


# ------------------------------------------------------------------ V1
lab = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t").set_index("genome_id")
p = pd.read_csv(ROOT / "results/external_validation_predictions.csv").query("arm == 'B2'")
M = calls_from_rows(pd.read_parquet(ROOT / "data/additional/metabolic_hmm/external.parquet"))
p = p[p.genome_id.isin(M.index) & ~p.genome_id.map(lab.faprotax_source.astype(bool))]
d = pd.DataFrame({"category": p.trait.to_numpy(), "y": p.y.astype(int).to_numpy(), "s": p.score.to_numpy(),
                  "c_M": [int(M.at[g, c]) for g, c in zip(p.genome_id, p.trait)],
                  "cl": lab.loc[p.genome_id, "genus"].to_numpy()})
d["c_R"] = (d.s >= d.category.map(thr)).astype(int)
cnt = d.groupby("category").y.agg(["sum", "size"])
cats = [c for c in cnt.index if cnt.at[c, "sum"] >= 5 and cnt.at[c, "size"] - cnt.at[c, "sum"] >= 5]
d = d[d.category.isin(cats)]


def v1stat(s):
    a, m = [], []
    for _, g in s.groupby("category"):
        if g.y.nunique() < 2:
            continue
        a.append(roc_auc_score(g.y, g.s) - roc_auc_score(g.y, g.c_M))
        m.append(matthews_corrcoef(g.y, g.c_R) - matthews_corrcoef(g.y, g.c_M))
    return [np.mean(a), np.mean(m)]


per = []
for c, g in d.groupby("category"):
    per.append(dict(category=c, n=len(g), pos=int(g.y.sum()), auroc_R=roc_auc_score(g.y, g.s),
                    auroc_M=roc_auc_score(g.y, g.c_M), mcc_R=matthews_corrcoef(g.y, g.c_R),
                    mcc_M=matthews_corrcoef(g.y, g.c_M)))
per = pd.DataFrame(per)
bs = boot(d, "cl", v1stat)
pt = v1stat(d)
y = d.y.to_numpy()
v1 = {"genomes": int(p.genome_id.nunique()), "pairs": int(len(d)), "categories": cats,
      "macro_auroc_RedoxKO": float(per.auroc_R.mean()), "macro_auroc_METABOLIC": float(per.auroc_M.mean()),
      "d_macro_auroc": float(pt[0]), "d_macro_auroc_ci": ci(bs[:, 0]), "success": bool(ci(bs[:, 0])[0] > 0),
      "macro_mcc_RedoxKO": float(per.mcc_R.mean()), "macro_mcc_METABOLIC": float(per.mcc_M.mean()),
      "d_macro_mcc": float(pt[1]), "d_macro_mcc_ci": ci(bs[:, 1]),
      "sens_R": float((d.c_R[y == 1] == 1).mean()), "spec_R": float((d.c_R[y == 0] == 0).mean()),
      "sens_M": float((d.c_M[y == 1] == 1).mean()), "spec_M": float((d.c_M[y == 0] == 0).mean()),
      "per_category": per.round(3).to_dict(orient="records")}
(OUT / "v1.json").write_text(json.dumps(v1, indent=2))
print(json.dumps({k: v for k, v in v1.items() if k != "per_category"}, indent=1, default=float))
print(per.round(3).to_string())

# ------------------------------------------------------------------ V2
F = ROOT / "results/additional/full_profile"
kos = (F / "ko_list.txt").read_text().split()
X = pd.concat([pd.read_parquet(f) for f in sorted((F / "features").glob("chunk_*.parquet"))] +
              [pd.read_parquet(ROOT / "validation_new_organisms/pathway100_ko_matrix.parquet")])
X = (X[~X.index.duplicated()].reindex(columns=kos).fillna(0) > 0)
labt, _, _ = load_labelled()
T = X.loc[labt.genome_id].to_numpy(dtype=bool)
x1 = pd.read_csv(ROOT / "results/additional/x1_calls.csv")
accs = x1.accession.unique()
Q = X.loc[accs].to_numpy(dtype=bool)
inter = Q.astype(np.int32) @ T.T.astype(np.int32)
union = Q.sum(1)[:, None] + T.sum(1)[None, :] - inter
dist = pd.Series(1 - (inter / union).max(1), index=accs, name="distance")
dist.to_csv(OUT / "v2_distances.csv")
tert = pd.qcut(dist.rank(method="first"), 3, labels=["near", "middle", "far"])
n = pd.DataFrame({"acc": x1.accession, "category": x1.category, "y": x1.expected.astype(int), "s": x1.score_B2,
                  "c_R": x1.call_B2.astype(int), "c_M": x1.call_METABOLIC.astype(int), "cl": x1.paper})
n["dist"], n["tert"] = n.acc.map(dist), n.acc.map(tert)


def bacc(y, c):
    return 0.5 * ((c[y == 1] == 1).mean() + (c[y == 0] == 0).mean()) if len(set(y)) == 2 else np.nan


def macro_auc(s):
    v = [roc_auc_score(g.y, g.s) for _, g in s.groupby("category") if g.y.nunique() == 2]
    return float(np.mean(v)) if v else np.nan


def v2stat(s):
    out = []
    for t in ("near", "middle", "far"):
        g = s[s.tert == t]
        out += [bacc(g.y.to_numpy(), g.c_R.to_numpy()), bacc(g.y.to_numpy(), g.c_M.to_numpy()), macro_auc(g)]
    near, far = out[0:3], out[6:9]
    return out + [far[0] - near[0], far[1] - near[1], (far[0] - near[0]) - (far[1] - near[1]), far[2] - near[2]]


pt, bs = v2stat(n), boot(n, "cl", v2stat)
names = [f"{t}_{m}" for t in ("near", "middle", "far") for m in ("bacc_RedoxKO", "bacc_METABOLIC", "macro_auroc_RedoxKO")]
names += ["far_minus_near_RedoxKO", "far_minus_near_METABOLIC", "diff_of_diffs_R_minus_M", "far_minus_near_auroc_RedoxKO"]
v2 = {k: {"value": float(v), "ci": ci(bs[:, i])} for i, (k, v) in enumerate(zip(names, pt))}
v2["tertile_distance_ranges"] = {t: [float(dist[tert == t].min()), float(dist[tert == t].max())] for t in ("near", "middle", "far")}
v2["n_tests_per_tertile"] = n.tert.value_counts().to_dict()
for tool, col in (("RedoxKO", "c_R"), ("METABOLIC", "c_M")):
    r = spearmanr(n.dist, (n[col] == n.y).astype(int))
    v2[f"spearman_correct_vs_distance_{tool}"] = {"rho": float(r.statistic), "p": float(r.pvalue)}
(OUT / "v2.json").write_text(json.dumps(v2, indent=2, default=str))
print(json.dumps(v2, indent=1, default=str))
