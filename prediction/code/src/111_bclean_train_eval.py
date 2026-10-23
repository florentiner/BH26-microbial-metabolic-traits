"""B-clean (docs/preregistration_Bclean_A.md): retrain the KO model (B2 procedure, unchanged) for acceptor_nitrate and
acceptor_CO2_methanogen_acetogen on the 358 training genomes plus the added non-FAPROTAX genomes, then test on the
newly described species (H1-H3) and in genus-blocked cross-validation (H4).
Released models are not touched; outputs go to results/models_bclean/ and results/additional/bclean/.
"""
import importlib
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import load_features, load_labelled, safe_ap  # noqa: E402

tf = importlib.import_module("09_train_final")
arms = importlib.import_module("06_run_arms")
OUT = ROOT / "results/additional/bclean"
SEED = 20261017
CATS = {"acceptor_nitrate": ["nitrate_reduction", "denitrification", "dnra"],
        "acceptor_CO2_methanogen_acetogen": ["methanogenesis", "acetogenesis"]}

# ---------------------------------------------------------------- data
lab, traits, Y = load_labelled()
X0 = load_features("energy_panel", lab.genome_id)
panel = list(X0.columns)
add = pd.read_csv(OUT / "added_genomes.tsv", sep="\t")
pw = pd.read_csv(ROOT / "data/pathways/labels.tsv", sep="\t").set_index("genome_id")
feats = pd.concat([pd.read_parquet(f) for f in sorted((OUT / "features").glob("chunk_*.parquet"))])
feats = (feats > 0).astype(np.float32).reindex(columns=panel).fillna(0.0)
add = add[add.genome_id.isin(feats.index)].reset_index(drop=True)
Xa = feats.loc[add.genome_id]
train_genus_group = dict(zip(lab.genus.astype(str).str.lower(), lab.group_genus))
ga = np.array([train_genus_group.get(g.lower(), f"add:{g}") for g in add.genus])
print(f"added genomes with features: {len(add)} (of {len(pd.read_csv(OUT / 'added_genomes.tsv', sep=chr(9)))})")

Xall = pd.concat([X0, Xa])
gall = np.r_[lab.group_genus.to_numpy(), ga]
yadd = {c: (pw.loc[add.genome_id, cols].fillna(0).max(axis=1).to_numpy() > 0).astype(float) for c, cols in CATS.items()}
yall = {c: np.r_[Y[:, traits.index(c)], yadd[c]] for c in CATS}
summary = {"n_added_with_features": int(len(add)),
           "added_positives": {c: int(yadd[c].sum()) for c in CATS},
           "training_positives": {c: int(Y[:, traits.index(c)].sum()) for c in CATS}}


def varying(df):
    return df.loc[:, df.std(axis=0) > 0]


def fit(Xdf, y, groups):
    Xv = varying(Xdf)
    cols = list(Xv.columns)
    return tf.select_and_fit(Xv, y, groups, arms.ko_configs(), lambda cfg, _c=cols: _c)


# ---------------------------------------------------------------- final B-clean models
bundle = {}
for c in CATS:
    model, thr, cfg, inner, cal = fit(Xall, yall[c], gall)
    bundle[c] = {"type": "estimator", "model": model, "threshold": thr, "config": cfg,
                 "inner_cv_auprc_all_labelled": inner, "calibrator": cal}
    print(c, cfg, f"thr={thr:.3f}", f"inner AUPRC={inner:.3f}", flush=True)
(ROOT / "results/models_bclean").mkdir(parents=True, exist_ok=True)
joblib.dump(bundle, ROOT / "results/models_bclean/model.joblib")

# ---------------------------------------------------------------- new species (H1-H3)
val = pd.read_parquet(ROOT / "validation_new_organisms/ko_panel_validation.parquet")
val = (val > 0).astype(np.float32)
b2 = joblib.load(ROOT / "results/models/B2/model.joblib")


def score_new(calls):
    rows = calls.copy()
    for c in CATS:
        sel = rows.category == c
        Xv = val.reindex(index=rows.loc[sel, "accession"], columns=panel).fillna(0.0)
        rows.loc[sel, "score_Bclean"] = bundle[c]["model"].predict_proba(Xv)
        rows.loc[sel, "call_Bclean"] = (rows.loc[sel, "score_Bclean"] >= bundle[c]["threshold"]).astype(int)
    return rows[rows.category.isin(CATS)]


def paired(d, B=2000):
    groups = {p: g for p, g in d.groupby("paper")}
    papers = np.array(list(groups))
    rng = np.random.default_rng(SEED)

    def stats(s):
        if s.expected.nunique() < 2:
            return None
        return (roc_auc_score(s.expected, s.score_Bclean) - roc_auc_score(s.expected, s.score_B2),
                matthews_corrcoef(s.expected, s.call_Bclean) - matthews_corrcoef(s.expected, s.call_B2))
    pt = stats(d)
    bs = [r for r in (stats(pd.concat([groups[p] for p in rng.choice(papers, len(papers))])) for _ in range(B)) if r]
    bs = np.array(bs)
    return {"n_tests": int(len(d)), "n_pos": int(d.expected.sum()), "n_papers": int(len(papers)),
            "auroc_B2": roc_auc_score(d.expected, d.score_B2), "auroc_Bclean": roc_auc_score(d.expected, d.score_Bclean),
            "mcc_B2": matthews_corrcoef(d.expected, d.call_B2), "mcc_Bclean": matthews_corrcoef(d.expected, d.call_Bclean),
            "d_auroc": pt[0], "d_auroc_ci": list(np.percentile(bs[:, 0], [2.5, 97.5])),
            "d_mcc": pt[1], "d_mcc_ci": list(np.percentile(bs[:, 1], [2.5, 97.5]))}


prim = score_new(pd.read_csv(ROOT / "results/additional/x1_calls.csv"))
hold = score_new(pd.read_csv(ROOT / "results/additional/hybrid/holdout_calls.csv"))
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

# ---------------------------------------------------------------- H4: genus-blocked CV on frozen folds
oof = pd.read_parquet(ROOT / "results/oof/B2.parquet")
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
        m, *_ = fit(Xtr, ytr, gtr)
        pred[te] = m.predict_proba(X0[te])
    y = Y[:, traits.index(c)]
    base = oof.loc[c].reindex(lab.genome_id).to_numpy()
    g = lab.group_genus.to_numpy()
    ug = np.unique(g)
    rng = np.random.default_rng(SEED)
    diffs = []
    for _ in range(2000):
        idx = np.concatenate([np.where(g == u)[0] for u in rng.choice(ug, len(ug))])
        if y[idx].sum() > 0:
            diffs.append(safe_ap(y[idx], pred[idx]) - safe_ap(y[idx], base[idx]))
    cv[c] = {"auprc_B2": safe_ap(y, base), "auprc_Bclean": safe_ap(y, pred),
             "d_auprc": safe_ap(y, pred) - safe_ap(y, base), "d_auprc_ci": list(np.percentile(diffs, [2.5, 97.5])),
             "non_inferior_margin_-0.02": bool(np.percentile(diffs, 2.5) > -0.02)}
    print(c, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in cv[c].items()}, flush=True)

summary.update({"newspecies": res, "cv_genus": cv})
(OUT / "bclean_results.json").write_text(json.dumps(summary, indent=2, default=float))
print(json.dumps(summary, indent=2, default=lambda x: round(float(x), 4)))
