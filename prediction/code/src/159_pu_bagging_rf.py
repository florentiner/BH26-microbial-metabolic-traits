"""Experiment D (src/112_pu_bagging.py, docs/preregistration_D_pu_bagging.md) for the final model RedoxKO: bagging
positive-unlabeled learning with the random forest as base learner, using the released forest configuration per
category (results/models/RF/model.joblib: min_samples_leaf; 500 trees x 3 seeds, balanced subsample), compared with the
released forest on the newly described species (HD1-HD3; score_RF / call_RF of results/additional/rf_paper/x1_calls_rf.csv
and holdout_calls_rf.csv) and in genus-blocked cross-validation (HD4; the forest's genus-level out-of-fold scores,
results/trees_xai/da/oof/RF.parquet). Same bags (K = 20, all positives + equal-size draw of unlabelled, seed 20261017 + b),
inner folds, threshold and calibrator rules, paper bootstrap (2000, seed 20261017) and genus bootstrap as the original.
usage: python src/159_pu_bagging_rf.py [--smoke] [--jobs N]
Outputs: results/models_pu_rf/model.joblib, results/additional/pu_rf/*
"""
import argparse
import importlib
import json
import os
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import (best_mcc_threshold, fit_calibrator, inner_folds, load_features,  # noqa: E402
                          load_labelled, safe_ap)

tx = importlib.import_module("47_trees_xai")
tx.TREE_JOBS = 1
ap = argparse.ArgumentParser()
ap.add_argument("--smoke", action="store_true")
ap.add_argument("--jobs", type=int, default=10)
args = ap.parse_args()
SEED, K = 20261017, 20
NB = 2000
if args.smoke:   # tiny forests, 2 bags, 50 bootstrap resamples; outputs to $SMOKE_DIR
    _make = tx.make
    tx.make = lambda cfg, seed: _make(cfg, seed).set_params(n_estimators=10)
    K, NB = 2, 50
SCRATCH = Path(os.environ.get("SMOKE_DIR", "/tmp/rf_smoke"))
OUT = SCRATCH / "pu_rf" if args.smoke else ROOT / "results/additional/pu_rf"
MODELS = SCRATCH / "models_pu_rf" if args.smoke else ROOT / "results/models_pu_rf"
OUT.mkdir(parents=True, exist_ok=True)

lab, traits, Y = load_labelled()
X0 = load_features("energy_panel", lab.genome_id)
X0 = X0.loc[:, X0.std(axis=0) > 0]
cols = list(X0.columns)
groups = lab.group_genus.to_numpy()
rf = joblib.load(ROOT / "results/models/RF/model.joblib")
assert all(list(rf[t]["model"].columns) == cols for t in traits)


def fit_bag(X, y, cfg):
    pos, unl = np.where(y == 1)[0], np.where(y == 0)[0]
    n = min(len(pos), len(unl))
    models = []
    for b in range(K):
        rng = np.random.default_rng(SEED + b)
        idx = np.r_[pos, rng.choice(unl, n, replace=False)]
        models += tx.fit_all(cfg, X[idx], y[idx])
    return models


def predict(models, X):
    return np.mean([m.predict_proba(X)[:, 1] for m in models], axis=0)


def one_category(t):
    j = traits.index(t)
    y = Y[:, j]
    X = X0.to_numpy(dtype=np.float32)
    cfg = rf[t]["config"]
    oof = np.zeros(len(y))
    for itr, iva in inner_folds(y, groups):
        oof[iva] = predict(fit_bag(X[itr], y[itr], cfg), X[iva])
    final = fit_bag(X, y, cfg)
    thr = best_mcc_threshold(y, oof)
    cal = fit_calibrator(y, oof)
    outer = np.zeros(len(y))
    fold = lab.fold_genus.to_numpy()
    for k in np.unique(fold):
        te = fold == k
        outer[te] = predict(fit_bag(X[~te], y[~te], cfg), X[te])
    print(t, "done", flush=True)
    return t, {"models": final, "threshold": thr, "calibrator": cal, "config": cfg, "inner_auprc": safe_ap(y, oof)}, outer


res = Parallel(n_jobs=args.jobs)(delayed(one_category)(t) for t in traits)
bundle = {t: b for t, b, _ in res}
outer = {t: o for t, _, o in res}
MODELS.mkdir(parents=True, exist_ok=True)
joblib.dump(bundle, MODELS / "model.joblib", compress=3)
pd.concat([pd.DataFrame({"genome_id": lab.genome_id, "trait": t, "y": Y[:, traits.index(t)], "prob_PU_RF": outer[t]})
           for t in traits]).to_csv(OUT / "cv_genus_oof.csv", index=False)

# ---------------------------------------------------------------- new species
val = (pd.read_parquet(ROOT / "validation_new_organisms/ko_panel_validation.parquet") > 0).astype(np.float32)


def score(calls):
    d = calls.copy()
    for t in d.category.unique():
        sel = d.category == t
        Xv = val.reindex(index=d.loc[sel, "accession"], columns=cols).fillna(0.0).to_numpy(dtype=np.float32)
        d.loc[sel, "score_PU_RF"] = predict(bundle[t]["models"], Xv)
        d.loc[sel, "call_PU_RF"] = (d.loc[sel, "score_PU_RF"] >= bundle[t]["threshold"]).astype(int)
    return d


def macro_auroc(df, col):
    v = [roc_auc_score(g.expected, g[col]) for _, g in df.groupby("category") if g.expected.nunique() == 2]
    return float(np.mean(v)) if v else np.nan


def sens_spec(df, col):
    y, c = df.expected.to_numpy(), df[col].to_numpy()
    return (c[y == 1] == 1).mean(), (c[y == 0] == 0).mean()


def compare(d, B=NB):
    groups_p = {p: g for p, g in d.groupby("paper")}
    papers = np.array(list(groups_p))
    rng = np.random.default_rng(SEED)

    def stat(s):
        return (macro_auroc(s, "score_PU_RF") - macro_auroc(s, "score_RF"),
                matthews_corrcoef(s.expected, s.call_PU_RF) - matthews_corrcoef(s.expected, s.call_RF))
    pt = stat(d)
    bs = np.array([stat(pd.concat([groups_p[p] for p in rng.choice(papers, len(papers))])) for _ in range(B)])
    se_pu, sp_pu = sens_spec(d, "call_PU_RF")
    se_rf, sp_rf = sens_spec(d, "call_RF")
    return {"n_tests": int(len(d)), "macro_auroc_RF": macro_auroc(d, "score_RF"), "macro_auroc_PU_RF": macro_auroc(d, "score_PU_RF"),
            "d_macro_auroc": pt[0], "d_macro_auroc_ci": list(np.nanpercentile(bs[:, 0], [2.5, 97.5])),
            "mcc_RF": matthews_corrcoef(d.expected, d.call_RF), "mcc_PU_RF": matthews_corrcoef(d.expected, d.call_PU_RF),
            "d_mcc": pt[1], "d_mcc_ci": list(np.nanpercentile(bs[:, 1], [2.5, 97.5])),
            "sens_RF": se_rf, "sens_PU_RF": se_pu, "spec_RF": sp_rf, "spec_PU_RF": sp_pu}


prim = score(pd.read_csv(ROOT / "results/additional/rf_paper/x1_calls_rf.csv"))
hold = score(pd.read_csv(ROOT / "results/additional/rf_paper/holdout_calls_rf.csv"))
prim.to_csv(OUT / "newspecies_primary_scores.csv", index=False)
hold.to_csv(OUT / "newspecies_withheld_scores.csv", index=False)
results = {"base_learner": "random forest (released RF config per category)", "primary100": compare(prim),
           "withheld43": compare(hold)}

# ---------------------------------------------------------------- HD4: internal genus CV, macro AUPRC
oof = pd.read_parquet(ROOT / "results/trees_xai/da/oof/RF.parquet")
oof = oof[oof.level == "genus"].groupby(["trait", "genome_id"]).prob.mean()
base = {t: oof.loc[t].reindex(lab.genome_id).to_numpy() for t in traits}
ug = np.unique(groups)
rng = np.random.default_rng(SEED)


def macro_ap(idx, src):
    v = [safe_ap(Y[idx, traits.index(t)], src[t][idx]) for t in traits if Y[idx, traits.index(t)].sum() > 0]
    return float(np.nanmean(v))


allidx = np.arange(len(lab))
diffs = []
for _ in range(NB):
    idx = np.concatenate([np.where(groups == u)[0] for u in rng.choice(ug, len(ug))])
    diffs.append(macro_ap(idx, outer) - macro_ap(idx, base))
results["cv_genus"] = {"macro_auprc_RF": macro_ap(allidx, base), "macro_auprc_PU_RF": macro_ap(allidx, outer),
                       "d": macro_ap(allidx, outer) - macro_ap(allidx, base),
                       "d_ci": list(np.percentile(diffs, [2.5, 97.5])),
                       "non_inferior_margin_-0.02": bool(np.percentile(diffs, 2.5) > -0.02)}
(OUT / "pu_results.json").write_text(json.dumps(results, indent=2, default=float))
print(json.dumps(results, indent=2, default=lambda x: round(float(x), 4)))
