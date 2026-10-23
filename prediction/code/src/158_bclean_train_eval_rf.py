"""B-clean (src/111_bclean_train_eval.py, docs/preregistration_Bclean_A.md) for the final model RedoxKO (random forest on
the 852-KO panel): the forest procedure (src/143_final_rf.py: min_samples_leaf 1 or 5 by grouped inner CV AUPRC,
500 trees x 3 seeds, balanced subsample, threshold = best inner-CV MCC, Platt calibrator) retrained for acceptor_nitrate
and acceptor_CO2_methanogen_acetogen on the 358 training genomes plus the 546 added non-FAPROTAX genomes, then tested
on the newly described species (H1-H3) and in genus-blocked cross-validation (H4). Same data, added genomes, labels,
tests, paper bootstrap (2000, seed 20261017) and genus bootstrap as the original; the comparator is the released
forest (score_RF / call_RF of results/additional/rf_paper/x1_calls_rf.csv and holdout_calls_rf.csv; for H4 its
genus-level out-of-fold scores, results/trees_xai/da/oof/RF.parquet) instead of B2.
usage: python src/158_bclean_train_eval_rf.py [--smoke]
Outputs: results/models_bclean_rf/model.joblib, results/additional/bclean_rf/
"""
import importlib
import json
import os
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import (FittedModel, best_mcc_threshold, fit_calibrator, inner_folds, load_features,  # noqa: E402
                          load_labelled, safe_ap)

tx = importlib.import_module("47_trees_xai")
tx.TREE_JOBS = 10
SMOKE = "--smoke" in sys.argv
if SMOKE:   # tiny forests and bootstrap, outputs to $SMOKE_DIR
    _make = tx.make
    tx.make = lambda cfg, seed: _make(cfg, seed).set_params(n_estimators=10)
SRC = ROOT / "results/additional/bclean"
SCRATCH = Path(os.environ.get("SMOKE_DIR", "/tmp/rf_smoke"))   # smoke-test outputs stay out of results/
OUT = SCRATCH / "bclean_rf" if SMOKE else ROOT / "results/additional/bclean_rf"
MODELS = SCRATCH / "models_bclean_rf" if SMOKE else ROOT / "results/models_bclean_rf"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261017
NB = 50 if SMOKE else 2000
CATS = {"acceptor_nitrate": ["nitrate_reduction", "denitrification", "dnra"],
        "acceptor_CO2_methanogen_acetogen": ["methanogenesis", "acetogenesis"]}


def rf_select_and_fit(Xdf, y, groups):
    """src/143_final_rf.py procedure on one category: returns model, threshold, config, inner AUPRC, calibrator."""
    Xv = Xdf.loc[:, Xdf.std(axis=0) > 0]
    cols = list(Xv.columns)
    X = Xv.to_numpy(dtype=np.float32)
    folds = inner_folds(y, groups)
    best = None
    for cfg in tx.CONFIGS["rf"]:
        oof = np.zeros(len(y))
        for itr, iva in folds:
            oof[iva] = tx.fit_predict(cfg, X[itr], y[itr], X[iva])
        ap = safe_ap(y, oof)
        if best is None or ap > best[0]:
            best = (ap, cfg, oof)
    ap, cfg, oof = best
    return FittedModel(cols, cfg, tx.fit_all(cfg, X, y)), best_mcc_threshold(y, oof), cfg, ap, fit_calibrator(y, oof)


# ---------------------------------------------------------------- data (as src/111)
lab, traits, Y = load_labelled()
X0 = load_features("energy_panel", lab.genome_id)
panel = list(X0.columns)
add = pd.read_csv(SRC / "added_genomes.tsv", sep="\t")
pw = pd.read_csv(ROOT / "data/pathways/labels.tsv", sep="\t").set_index("genome_id")
feats = pd.concat([pd.read_parquet(f) for f in sorted((SRC / "features").glob("chunk_*.parquet"))])
feats = (feats > 0).astype(np.float32).reindex(columns=panel).fillna(0.0)
add = add[add.genome_id.isin(feats.index)].reset_index(drop=True)
Xa = feats.loc[add.genome_id]
train_genus_group = dict(zip(lab.genus.astype(str).str.lower(), lab.group_genus))
ga = np.array([train_genus_group.get(g.lower(), f"add:{g}") for g in add.genus])
print(f"added genomes with features: {len(add)} (of {len(pd.read_csv(SRC / 'added_genomes.tsv', sep=chr(9)))})")
Xall = pd.concat([X0, Xa])
gall = np.r_[lab.group_genus.to_numpy(), ga]
yadd = {c: (pw.loc[add.genome_id, cols].fillna(0).max(axis=1).to_numpy() > 0).astype(float) for c, cols in CATS.items()}
yall = {c: np.r_[Y[:, traits.index(c)], yadd[c]] for c in CATS}
summary = {"model": "random forest (RedoxKO procedure, src/143_final_rf.py)", "n_added_with_features": int(len(add)),
           "added_positives": {c: int(yadd[c].sum()) for c in CATS},
           "training_positives": {c: int(Y[:, traits.index(c)].sum()) for c in CATS}}

# ---------------------------------------------------------------- final B-clean forests
bundle = {}
for c in CATS:
    model, thr, cfg, inner, cal = rf_select_and_fit(Xall, yall[c], gall)
    bundle[c] = {"type": "estimator", "model": model, "threshold": thr, "config": cfg,
                 "inner_cv_auprc_all_labelled": inner, "calibrator": cal}
    print(c, cfg, f"thr={thr:.3f}", f"inner AUPRC={inner:.3f}", flush=True)
MODELS.mkdir(parents=True, exist_ok=True)
joblib.dump(bundle, MODELS / "model.joblib", compress=3)

# ---------------------------------------------------------------- new species (H1-H3)
val = (pd.read_parquet(ROOT / "validation_new_organisms/ko_panel_validation.parquet") > 0).astype(np.float32)


def score_new(calls):
    rows = calls.copy()
    for c in CATS:
        sel = rows.category == c
        Xv = val.reindex(index=rows.loc[sel, "accession"], columns=panel).fillna(0.0)
        rows.loc[sel, "score_Bclean_RF"] = bundle[c]["model"].predict_proba(Xv)
        rows.loc[sel, "call_Bclean_RF"] = (rows.loc[sel, "score_Bclean_RF"] >= bundle[c]["threshold"]).astype(int)
    return rows[rows.category.isin(CATS)]


def paired(d, B=NB):
    groups = {p: g for p, g in d.groupby("paper")}
    papers = np.array(list(groups))
    rng = np.random.default_rng(SEED)

    def stats(s):
        if s.expected.nunique() < 2:
            return None
        return (roc_auc_score(s.expected, s.score_Bclean_RF) - roc_auc_score(s.expected, s.score_RF),
                matthews_corrcoef(s.expected, s.call_Bclean_RF) - matthews_corrcoef(s.expected, s.call_RF))
    pt = stats(d)
    bs = np.array([r for r in (stats(pd.concat([groups[p] for p in rng.choice(papers, len(papers))])) for _ in range(B)) if r])
    return {"n_tests": int(len(d)), "n_pos": int(d.expected.sum()), "n_papers": int(len(papers)),
            "auroc_RF": roc_auc_score(d.expected, d.score_RF), "auroc_Bclean_RF": roc_auc_score(d.expected, d.score_Bclean_RF),
            "mcc_RF": matthews_corrcoef(d.expected, d.call_RF), "mcc_Bclean_RF": matthews_corrcoef(d.expected, d.call_Bclean_RF),
            "d_auroc": pt[0], "d_auroc_ci": list(np.percentile(bs[:, 0], [2.5, 97.5])),
            "d_mcc": pt[1], "d_mcc_ci": list(np.percentile(bs[:, 1], [2.5, 97.5]))}


prim = score_new(pd.read_csv(ROOT / "results/additional/rf_paper/x1_calls_rf.csv"))
hold = score_new(pd.read_csv(ROOT / "results/additional/rf_paper/holdout_calls_rf.csv"))
res = {}
for name, d in [("primary100", prim), ("withheld43", hold)]:
    for c in CATS:
        dc = d[d.category == c]
        if len(dc) and dc.expected.nunique() == 2:
            res[f"{name}:{c}"] = paired(dc)
        else:
            res[f"{name}:{c}"] = {"n_tests": int(len(dc)), "note": "too few tests or one class only"}
prim.to_csv(OUT / "newspecies_primary_scores.csv", index=False)
hold.to_csv(OUT / "newspecies_withheld_scores.csv", index=False)
print(json.dumps(res, indent=1, default=float), flush=True)

# ---------------------------------------------------------------- H4: genus-blocked CV on frozen folds
oof = pd.read_parquet(ROOT / "results/trees_xai/da/oof/RF.parquet")
oof = oof[oof.level == "genus"].groupby(["trait", "genome_id"]).prob.mean()
fold = lab.fold_genus.to_numpy()
cv = {}
for c in CATS:
    pred = np.zeros(len(lab))
    for k in np.unique(fold):
        te = fold == k
        test_groups = set(lab.group_genus[te])
        keep_add = np.array([g not in test_groups for g in ga])
        Xtr = pd.concat([X0[~te], Xa[keep_add]])
        ytr = np.r_[Y[~te, traits.index(c)], yadd[c][keep_add]]
        gtr = np.r_[lab.group_genus.to_numpy()[~te], ga[keep_add]]
        m, *_ = rf_select_and_fit(Xtr, ytr, gtr)
        pred[te] = m.predict_proba(X0[te])
        print(c, "outer fold", k, "done", flush=True)
    y = Y[:, traits.index(c)]
    base = oof.loc[c].reindex(lab.genome_id).to_numpy()
    g = lab.group_genus.to_numpy()
    ug = np.unique(g)
    rng = np.random.default_rng(SEED)
    diffs = []
    for _ in range(NB):
        idx = np.concatenate([np.where(g == u)[0] for u in rng.choice(ug, len(ug))])
        if y[idx].sum() > 0:
            diffs.append(safe_ap(y[idx], pred[idx]) - safe_ap(y[idx], base[idx]))
    cv[c] = {"auprc_RF": safe_ap(y, base), "auprc_Bclean_RF": safe_ap(y, pred),
             "d_auprc": safe_ap(y, pred) - safe_ap(y, base), "d_auprc_ci": list(np.percentile(diffs, [2.5, 97.5])),
             "non_inferior_margin_-0.02": bool(np.percentile(diffs, 2.5) > -0.02)}
    pd.DataFrame({"genome_id": lab.genome_id, "trait": c, "y": y, "prob_Bclean_RF": pred, "prob_RF": base}).to_csv(
        OUT / f"cv_genus_oof_{c}.csv", index=False)
    print(c, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in cv[c].items()}, flush=True)

summary.update({"newspecies": res, "cv_genus": cv})
(OUT / "bclean_results.json").write_text(json.dumps(summary, indent=2, default=float))
print(json.dumps(summary, indent=2, default=lambda x: round(float(x), 4)))
