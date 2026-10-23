"""Experiment F (src/117_full_profile_train_eval.py, docs/preregistration_F_full_ko_profile.md) for the final model
RedoxKO (random forest): the forest procedure of src/143_final_rf.py (min_samples_leaf 1 or 5 by grouped inner CV AUPRC,
500 trees x 3 seeds, balanced subsample, threshold = best inner-CV MCC, Platt calibrator) on the full KO profile
(10,517 EnergyKO-M KOs, varying in training) = "RF-full", compared with the released forest on the 852-KO panel
(score_RF / call_RF of results/additional/rf_paper/x1_calls_rf.csv and holdout_calls_rf.csv) on the newly described
species (HF1-HF3) and in genus-blocked cross-validation (HF4; panel forest = results/trees_xai/da/oof/RF.parquet, genus
level). Same features, tests, paper bootstrap (2000, seed 20261017) and genus bootstrap as the original.
usage: python src/160_full_profile_train_eval_rf.py [--smoke] [--jobs N]
Outputs: results/models_full_rf/model.joblib, results/additional/full_profile_rf/results.json and score files.
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
from redox_common import (FittedModel, best_mcc_threshold, fit_calibrator, inner_folds, load_labelled,  # noqa: E402
                          safe_ap)

tx = importlib.import_module("47_trees_xai")
tx.TREE_JOBS = 1
ap = argparse.ArgumentParser()
ap.add_argument("--smoke", action="store_true")
ap.add_argument("--jobs", type=int, default=10)
args = ap.parse_args()
SRC = ROOT / "results/additional/full_profile"
SEED, NB = 20261017, 2000
if args.smoke:
    _make = tx.make
    tx.make = lambda cfg, seed: _make(cfg, seed).set_params(n_estimators=10)
    NB = 50
SCRATCH = Path(os.environ.get("SMOKE_DIR", "/tmp/rf_smoke"))
OUT = SCRATCH / "full_profile_rf" if args.smoke else ROOT / "results/additional/full_profile_rf"
MODELS = SCRATCH / "models_full_rf" if args.smoke else ROOT / "results/models_full_rf"
OUT.mkdir(parents=True, exist_ok=True)

lab, traits, Y = load_labelled()
KOS = (SRC / "ko_list.txt").read_text().split()
F = pd.concat([pd.read_parquet(f) for f in sorted((SRC / "features").glob("chunk_*.parquet"))])
F = pd.concat([F, pd.read_parquet(ROOT / "validation_new_organisms/pathway100_ko_matrix.parquet")])
F = (F[~F.index.duplicated()].reindex(columns=KOS).fillna(0) > 0).astype(np.float32)
Xtr = F.loc[lab.genome_id]
Xtr = Xtr.loc[:, Xtr.std(axis=0) > 0]
cols = list(Xtr.columns)
groups = lab.group_genus.to_numpy()
print(f"training genomes {len(Xtr)}, varying KOs {len(cols)}", flush=True)


def fit(Xdf, y, g):
    """src/143_final_rf.py procedure on one category (columns: those varying in the training data)."""
    Xv = Xdf.loc[:, Xdf.std(axis=0) > 0]
    c = list(Xv.columns)
    X = Xv.to_numpy(dtype=np.float32)
    folds = inner_folds(y, g)
    best = None
    for cfg in tx.CONFIGS["rf"]:
        oof = np.zeros(len(y))
        for itr, iva in folds:
            oof[iva] = tx.fit_predict(cfg, X[itr], y[itr], X[iva])
        a = safe_ap(y, oof)
        if best is None or a > best[0]:
            best = (a, cfg, oof)
    a, cfg, oof = best
    return FittedModel(c, cfg, tx.fit_all(cfg, X, y)), best_mcc_threshold(y, oof), cfg, a, fit_calibrator(y, oof)


def final_one(j, t):
    model, thr, cfg, inner, cal = fit(Xtr, Y[:, j], groups)
    print("final", t, cfg, flush=True)
    return t, {"type": "estimator", "model": model, "threshold": thr, "config": cfg, "inner_cv_auprc": inner,
               "calibrator": cal}


bundle = dict(Parallel(n_jobs=args.jobs)(delayed(final_one)(j, t) for j, t in enumerate(traits)))
MODELS.mkdir(parents=True, exist_ok=True)
joblib.dump(bundle, MODELS / "model.joblib", compress=3)
print("final models:", {t: b["config"] for t, b in bundle.items()}, flush=True)


# ---------------------------------------------------------------- new species (HF1-HF3)
def score(calls):
    d = calls.copy()
    for t in d.category.unique():
        sel = d.category == t
        Xv = F.reindex(index=d.loc[sel, "accession"], columns=cols).fillna(0.0)
        assert not Xv.isna().any().any()
        d.loc[sel, "score_full_RF"] = bundle[t]["model"].predict_proba(Xv)
        d.loc[sel, "call_full_RF"] = (d.loc[sel, "score_full_RF"] >= bundle[t]["threshold"]).astype(int)
    return d


def macro(df, col):
    v = [roc_auc_score(g.expected, g[col]) for _, g in df.groupby("category") if g.expected.nunique() == 2]
    return float(np.mean(v)) if v else np.nan


def sens_spec(df, col):
    y, c = df.expected.to_numpy(), df[col].to_numpy()
    return float((c[y == 1] == 1).mean()), float((c[y == 0] == 0).mean())


def compare(d, B=NB):
    gp = {p: g for p, g in d.groupby("paper")}
    papers = np.array(list(gp))
    rng = np.random.default_rng(SEED)

    def stat(s):
        return (macro(s, "score_full_RF") - macro(s, "score_RF"),
                matthews_corrcoef(s.expected, s.call_full_RF) - matthews_corrcoef(s.expected, s.call_RF))
    pt = stat(d)
    bs = np.array([stat(pd.concat([gp[p] for p in rng.choice(papers, len(papers))])) for _ in range(B)])

    def ci(k):
        return [float(x) for x in np.nanpercentile(bs[:, k], [2.5, 97.5])]
    se_f, sp_f = sens_spec(d, "call_full_RF")
    se_b, sp_b = sens_spec(d, "call_RF")
    return {"n_tests": int(len(d)), "macro_auroc_RF_panel": macro(d, "score_RF"), "macro_auroc_RF_full": macro(d, "score_full_RF"),
            "d_macro_auroc": pt[0], "d_macro_auroc_ci": ci(0), "success_lower_gt_0": bool(ci(0)[0] > 0),
            "mcc_RF_panel": float(matthews_corrcoef(d.expected, d.call_RF)),
            "mcc_RF_full": float(matthews_corrcoef(d.expected, d.call_full_RF)), "d_mcc": pt[1], "d_mcc_ci": ci(1),
            "sens_RF_panel": se_b, "sens_RF_full": se_f, "spec_RF_panel": sp_b, "spec_RF_full": sp_f}


prim = score(pd.read_csv(ROOT / "results/additional/rf_paper/x1_calls_rf.csv"))
hold = score(pd.read_csv(ROOT / "results/additional/rf_paper/holdout_calls_rf.csv"))
prim.to_csv(OUT / "newspecies_primary_scores.csv", index=False)
hold.to_csv(OUT / "newspecies_withheld_scores.csv", index=False)
res = {"n_training_genomes": len(Xtr), "n_varying_kos": len(cols),
       "HF1_HF2_primary100": compare(prim), "HF3_withheld43": compare(hold)}
res["HF3_same_sign"] = {k: bool(np.sign(res["HF3_withheld43"][k]) == np.sign(res["HF1_HF2_primary100"][k]))
                        for k in ("d_macro_auroc", "d_mcc")}
print(json.dumps(res, indent=1, default=float), flush=True)
(OUT / "results.json").write_text(json.dumps(res, indent=2, default=float))   # interim, completed below

# ---------------------------------------------------------------- HF4: genus-blocked CV
fold = lab.fold_genus.to_numpy()


def outer_one(j, k):
    te = fold == k
    m, *_ = fit(Xtr[~te], Y[~te, j], groups[~te])
    print("outer", traits[j], k, flush=True)
    return j, k, m.predict_proba(Xtr[te])


outs = Parallel(n_jobs=args.jobs)(delayed(outer_one)(j, k) for j in range(len(traits)) for k in np.unique(fold))
pred = {t: np.zeros(len(lab)) for t in traits}
for j, k, p in outs:
    pred[traits[j]][fold == k] = p
pd.concat([pd.DataFrame({"genome_id": lab.genome_id, "trait": t, "y": Y[:, j], "prob_RF_full": pred[t]})
           for j, t in enumerate(traits)]).to_csv(OUT / "cv_genus_oof.csv", index=False)
oof = pd.read_parquet(ROOT / "results/trees_xai/da/oof/RF.parquet")
oof = oof[oof.level == "genus"].groupby(["trait", "genome_id"]).prob.mean()
base = {t: oof.loc[t].reindex(lab.genome_id).to_numpy() for t in traits}


def macro_ap(idx, src):
    v = [safe_ap(Y[idx, j], src[t][idx]) for j, t in enumerate(traits) if Y[idx, j].sum() > 0]
    return float(np.nanmean(v))


allidx = np.arange(len(lab))
ug = np.unique(groups)
rng = np.random.default_rng(SEED)
diffs = []
for _ in range(NB):
    idx = np.concatenate([np.where(groups == u)[0] for u in rng.choice(ug, len(ug))])
    diffs.append(macro_ap(idx, pred) - macro_ap(idx, base))
res["HF4_cv_genus"] = {"macro_auprc_RF_panel": macro_ap(allidx, base), "macro_auprc_RF_full": macro_ap(allidx, pred),
                       "d": macro_ap(allidx, pred) - macro_ap(allidx, base),
                       "d_ci": [float(x) for x in np.percentile(diffs, [2.5, 97.5])],
                       "non_inferior_margin_-0.02": bool(np.percentile(diffs, 2.5) > -0.02),
                       "per_trait": {t: {"RF_panel": safe_ap(Y[:, j], base[t]), "RF_full": safe_ap(Y[:, j], pred[t])}
                                     for j, t in enumerate(traits)}}
(OUT / "results.json").write_text(json.dumps(res, indent=2, default=float))
print(json.dumps(res["HF4_cv_genus"], indent=1, default=float))
