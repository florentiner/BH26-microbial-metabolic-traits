"""Experiment D (docs/preregistration_D_pu_bagging.md): bagging positive-unlabeled learning for the KO model, with the
released B2 configuration per category, compared with the released B2 on the newly described species (HD1-HD3) and in
genus-blocked cross-validation (HD4). Released models are not touched.
Outputs: results/models_pu/model.joblib, results/additional/pu/*.
"""
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import (best_mcc_threshold, fit_calibrator, fit_models, inner_folds, load_features,  # noqa: E402
                          load_labelled, safe_ap)

SEED, K = 20261017, 20
OUT = ROOT / "results/additional/pu"
OUT.mkdir(parents=True, exist_ok=True)

lab, traits, Y = load_labelled()
X0 = load_features("energy_panel", lab.genome_id)
X0 = X0.loc[:, X0.std(axis=0) > 0]
cols = list(X0.columns)
groups = lab.group_genus.to_numpy()
b2 = joblib.load(ROOT / "results/models/B2/model.joblib")


def fit_bag(X, y, cfg):
    pos, unl = np.where(y == 1)[0], np.where(y == 0)[0]
    n = min(len(pos), len(unl))
    models = []
    for b in range(K):
        rng = np.random.default_rng(SEED + b)
        idx = np.r_[pos, rng.choice(unl, n, replace=False)]
        models += fit_models(cfg, X[idx], y[idx])
    return models


def predict(models, X):
    return np.mean([m.predict_proba(X)[:, 1] for m in models], axis=0)


def one_category(t):
    j = traits.index(t)
    y = Y[:, j]
    X = X0.to_numpy(dtype=np.float32)
    cfg = b2[t]["config"]
    oof = np.zeros(len(y))
    for itr, iva in inner_folds(y, groups):
        oof[iva] = predict(fit_bag(X[itr], y[itr], cfg), X[iva])
    final = fit_bag(X, y, cfg)
    thr = best_mcc_threshold(y, oof)
    cal = fit_calibrator(y, oof)
    # HD4: outer genus-blocked folds
    outer = np.zeros(len(y))
    fold = lab.fold_genus.to_numpy()
    for k in np.unique(fold):
        te = fold == k
        outer[te] = predict(fit_bag(X[~te], y[~te], cfg), X[te])
    return t, {"models": final, "threshold": thr, "calibrator": cal, "config": cfg, "inner_auprc": safe_ap(y, oof)}, outer


res = Parallel(n_jobs=7)(delayed(one_category)(t) for t in traits)
bundle = {t: b for t, b, _ in res}
outer = {t: o for t, _, o in res}
(ROOT / "results/models_pu").mkdir(parents=True, exist_ok=True)
joblib.dump({t: {k: v for k, v in b.items()} for t, b in bundle.items()}, ROOT / "results/models_pu/model.joblib")

# ---------------------------------------------------------------- new species
val = (pd.read_parquet(ROOT / "validation_new_organisms/ko_panel_validation.parquet") > 0).astype(np.float32)


def score(calls):
    d = calls.copy()
    for t in d.category.unique():
        sel = d.category == t
        Xv = val.reindex(index=d.loc[sel, "accession"], columns=cols).fillna(0.0).to_numpy(dtype=np.float32)
        d.loc[sel, "score_PU"] = predict(bundle[t]["models"], Xv)
        d.loc[sel, "call_PU"] = (d.loc[sel, "score_PU"] >= bundle[t]["threshold"]).astype(int)
    return d


def macro_auroc(df, col):
    v = [roc_auc_score(g.expected, g[col]) for _, g in df.groupby("category") if g.expected.nunique() == 2]
    return float(np.mean(v)) if v else np.nan


def sens_spec(df, col):
    y, c = df.expected.to_numpy(), df[col].to_numpy()
    return (c[y == 1] == 1).mean(), (c[y == 0] == 0).mean()


def compare(d, B=2000):
    groups_p = {p: g for p, g in d.groupby("paper")}
    papers = np.array(list(groups_p))
    rng = np.random.default_rng(SEED)
    stat = lambda s: (macro_auroc(s, "score_PU") - macro_auroc(s, "score_B2"),
                      matthews_corrcoef(s.expected, s.call_PU) - matthews_corrcoef(s.expected, s.call_B2))
    pt = stat(d)
    bs = np.array([stat(pd.concat([groups_p[p] for p in rng.choice(papers, len(papers))])) for _ in range(B)])
    se_pu, sp_pu = sens_spec(d, "call_PU")
    se_b2, sp_b2 = sens_spec(d, "call_B2")
    return {"n_tests": int(len(d)), "macro_auroc_B2": macro_auroc(d, "score_B2"), "macro_auroc_PU": macro_auroc(d, "score_PU"),
            "d_macro_auroc": pt[0], "d_macro_auroc_ci": list(np.nanpercentile(bs[:, 0], [2.5, 97.5])),
            "mcc_B2": matthews_corrcoef(d.expected, d.call_B2), "mcc_PU": matthews_corrcoef(d.expected, d.call_PU),
            "d_mcc": pt[1], "d_mcc_ci": list(np.nanpercentile(bs[:, 1], [2.5, 97.5])),
            "sens_B2": se_b2, "sens_PU": se_pu, "spec_B2": sp_b2, "spec_PU": sp_pu}


prim = score(pd.read_csv(ROOT / "results/additional/x1_calls.csv"))
hold = score(pd.read_csv(ROOT / "results/additional/hybrid/holdout_calls.csv"))
prim.to_csv(OUT / "newspecies_primary_scores.csv", index=False)
hold.to_csv(OUT / "newspecies_withheld_scores.csv", index=False)
results = {"primary100": compare(prim), "withheld43": compare(hold)}

# ---------------------------------------------------------------- HD4: internal genus CV, macro AUPRC
oof = pd.read_parquet(ROOT / "results/oof/B2.parquet")
oof = oof[oof.level == "genus"].groupby(["trait", "genome_id"]).prob.mean()
base = {t: oof.loc[t].reindex(lab.genome_id).to_numpy() for t in traits}
ug = np.unique(groups)
rng = np.random.default_rng(SEED)


def macro_ap(idx, src):
    v = [safe_ap(Y[idx, traits.index(t)], src[t][idx]) for t in traits if Y[idx, traits.index(t)].sum() > 0]
    return float(np.nanmean(v))


allidx = np.arange(len(lab))
diffs = []
for _ in range(2000):
    idx = np.concatenate([np.where(groups == u)[0] for u in rng.choice(ug, len(ug))])
    diffs.append(macro_ap(idx, outer) - macro_ap(idx, base))
results["cv_genus"] = {"macro_auprc_B2": macro_ap(allidx, base), "macro_auprc_PU": macro_ap(allidx, outer),
                       "d": macro_ap(allidx, outer) - macro_ap(allidx, base),
                       "d_ci": list(np.percentile(diffs, [2.5, 97.5])),
                       "non_inferior_margin_-0.02": bool(np.percentile(diffs, 2.5) > -0.02)}
(OUT / "pu_results.json").write_text(json.dumps(results, indent=2, default=float))
print(json.dumps(results, indent=2, default=lambda x: round(float(x), 4)))
