"""Feature attribution of the released KO model RedoxKO (file code B2) itself, with the method used for the random
forests (src/47_trees_xai.py): refit RedoxKO in each outer genus-blocked fold with the configuration its own
cross-validation selected for that fold, compute SHAP values on the held-out genomes (exact linear SHAP for logistic
regression, TreeSHAP for XGBoost, averaged over the three seeds), call a KO stable if it ranks in the top ten by mean
|SHAP| in at least 60% of folds, and test enrichment of curated marker KOs among stable KOs (hypergeometric, BH across
traits). Also reports the overlap with the random-forest stable set. Descriptive analysis; no model is changed.
Outputs: results/additional/xai_redoxko/{stable.csv, enrichment.csv, overlap_rf.csv}
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import fit_models, load_features, load_labelled  # noqa: E402

xai = importlib.import_module("47_trees_xai")
m8 = importlib.import_module("08_metrics")
OUT = ROOT / "results/additional/xai_redoxko"
OUT.mkdir(parents=True, exist_ok=True)

lab, traits, Y = load_labelled()
Xdf = load_features("energy_panel", lab.genome_id)
Xdf = Xdf.loc[:, Xdf.std(axis=0) > 0]
X, cols = Xdf.to_numpy(dtype=np.float32), list(Xdf.columns)
fold = lab.fold_genus.to_numpy()
oof = pd.read_parquet(ROOT / "results/oof/B2.parquet")
oof = oof[oof.level == "genus"]
cg = xai.reference_groups()


def shap_values(models, Xte):
    import shap
    vals = []
    for m in models:
        if hasattr(m, "named_steps"):  # StandardScaler + LogisticRegression: exact linear SHAP on the logit scale
            sc, lr = m.named_steps["standardscaler"], m.named_steps["logisticregression"]
            z = (Xte - sc.mean_) / sc.scale_
            vals.append(z * lr.coef_[0])
        else:
            vals.append(shap.TreeExplainer(m).shap_values(Xte))
    return np.mean(vals, axis=0)


def one(j, t):
    dt = oof[oof.trait == t]
    rows = []
    for k in sorted(dt.fold.unique()):
        cfg = json.loads(dt[dt.fold == k].config.iloc[0])
        tr, te = fold != k, fold == k
        sv = shap_values(fit_models(cfg, X[tr], Y[tr, j]), X[te])
        imp = np.abs(sv).mean(axis=0)
        present = X[te] > 0
        direction = np.array([sv[present[:, c], c].mean() - sv[~present[:, c], c].mean()
                              if present[:, c].any() and (~present[:, c]).any() else np.nan for c in range(len(cols))])
        for rank, c in enumerate(np.argsort(-imp)[:10]):
            rows.append(dict(trait=t, fold=int(k), KO=cols[c], rank=rank + 1, mean_abs_shap=float(imp[c]),
                             shap_direction=float(direction[c]), model=cfg.get("model")))
    return rows


rows = sum(Parallel(n_jobs=4)(delayed(one)(j, t) for j, t in enumerate(traits)), [])
top = pd.DataFrame(rows)
nf = top.groupby("trait").fold.nunique()
stab = top.groupby(["trait", "KO"]).agg(n_folds_top10=("fold", "nunique"), mean_rank=("rank", "mean"),
                                         mean_abs_shap=("mean_abs_shap", "mean"),
                                         shap_direction=("shap_direction", "mean")).reset_index()
stab["stability"] = stab.n_folds_top10 / stab.trait.map(nf)
stable = stab[stab.stability >= 0.6].copy()
N = len(cols)
enr = []
for t in traits:
    known = xai.known_markers("da", t, cg) & set(cols)
    s = stable[stable.trait == t]
    stable.loc[stable.trait == t, "known_marker"] = stable.loc[stable.trait == t, "KO"].isin(known)
    x, n, K = int(s.KO.isin(known).sum()), len(s), len(known)
    enr.append(dict(trait=t, n_known_markers_in_features=K, n_stable_top=n, n_stable_top_known=x,
                    p_hypergeom=hypergeom.sf(x - 1, N, K, n) if n and K else np.nan,
                    stable_top_kos=", ".join(s.sort_values("mean_rank").KO)))
enr = pd.DataFrame(enr)
ok = enr.p_hypergeom.notna()
enr.loc[ok, "q_bh"] = m8.bh(enr.loc[ok, "p_hypergeom"].to_numpy())
rf = pd.read_csv(ROOT / "results/trees_xai/da/xai_stable_rf.csv")
ov = []
for t in traits:
    a, b = set(stable[stable.trait == t].KO), set(rf[rf.trait == t].KO)
    ov.append(dict(trait=t, n_redoxko=len(a), n_rf=len(b), n_shared=len(a & b),
                   jaccard=len(a & b) / len(a | b) if a | b else np.nan))
stable.sort_values(["trait", "mean_rank"]).to_csv(OUT / "stable.csv", index=False)
enr.to_csv(OUT / "enrichment.csv", index=False)
pd.DataFrame(ov).to_csv(OUT / "overlap_rf.csv", index=False)
print(enr[["trait", "n_known_markers_in_features", "n_stable_top", "n_stable_top_known", "q_bh"]].to_string(index=False))
print(pd.DataFrame(ov).to_string(index=False))
