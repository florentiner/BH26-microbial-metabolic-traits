"""Experiment F (docs/preregistration_F_full_ko_profile.md): RedoxKO-full, the released RedoxKO procedure on the full
KO profile (10,517 EnergyKO-M KOs, varying in training), compared with RedoxKO on the newly described species
(HF1-HF3) and in genus-blocked cross-validation (HF4). Released models are not touched.
Outputs: results/models_full/model.joblib, results/additional/full_profile/results.json and score files.
"""
import importlib
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
from redox_common import load_labelled, safe_ap  # noqa: E402

tf = importlib.import_module("09_train_final")
arms = importlib.import_module("06_run_arms")
OUT = ROOT / "results/additional/full_profile"
SEED = 20261017

lab, traits, Y = load_labelled()
KOS = (OUT / "ko_list.txt").read_text().split()
F = pd.concat([pd.read_parquet(f) for f in sorted((OUT / "features").glob("chunk_*.parquet"))])
F = pd.concat([F, pd.read_parquet(ROOT / "validation_new_organisms/pathway100_ko_matrix.parquet")])
F = (F[~F.index.duplicated()].reindex(columns=KOS).fillna(0) > 0).astype(np.float32)
Xtr = F.loc[lab.genome_id]
Xtr = Xtr.loc[:, Xtr.std(axis=0) > 0]
cols = list(Xtr.columns)
groups = lab.group_genus.to_numpy()
print(f"training genomes {len(Xtr)}, varying KOs {len(cols)}", flush=True)


def fit(X, y, g):
    return tf.select_and_fit(X, y, g, arms.ko_configs(), lambda cfg, _c=list(X.columns): _c)


def final_one(j, t):
    model, thr, cfg, inner, cal = fit(Xtr, Y[:, j], groups)
    return t, {"type": "estimator", "model": model, "threshold": thr, "config": cfg, "inner_cv_auprc": inner,
               "calibrator": cal}


bundle = dict(Parallel(n_jobs=7)(delayed(final_one)(j, t) for j, t in enumerate(traits)))
(ROOT / "results/models_full").mkdir(parents=True, exist_ok=True)
joblib.dump(bundle, ROOT / "results/models_full/model.joblib")
print("final models:", {t: b["config"] for t, b in bundle.items()}, flush=True)


# ---------------------------------------------------------------- new species (HF1-HF3)
def score(calls):
    d = calls.copy()
    for t in d.category.unique():
        sel = d.category == t
        Xv = F.reindex(index=d.loc[sel, "accession"], columns=cols).fillna(0.0)
        assert not Xv.isna().any().any()
        d.loc[sel, "score_full"] = bundle[t]["model"].predict_proba(Xv)
        d.loc[sel, "call_full"] = (d.loc[sel, "score_full"] >= bundle[t]["threshold"]).astype(int)
    return d


def macro(df, col):
    v = [roc_auc_score(g.expected, g[col]) for _, g in df.groupby("category") if g.expected.nunique() == 2]
    return float(np.mean(v)) if v else np.nan


def sens_spec(df, col):
    y, c = df.expected.to_numpy(), df[col].to_numpy()
    return float((c[y == 1] == 1).mean()), float((c[y == 0] == 0).mean())


def compare(d, B=2000):
    gp = {p: g for p, g in d.groupby("paper")}
    papers = np.array(list(gp))
    rng = np.random.default_rng(SEED)
    stat = lambda s: (macro(s, "score_full") - macro(s, "score_B2"),
                      matthews_corrcoef(s.expected, s.call_full) - matthews_corrcoef(s.expected, s.call_B2))
    pt = stat(d)
    bs = np.array([stat(pd.concat([gp[p] for p in rng.choice(papers, len(papers))])) for _ in range(B)])
    ci = lambda k: [float(x) for x in np.nanpercentile(bs[:, k], [2.5, 97.5])]
    se_f, sp_f = sens_spec(d, "call_full")
    se_b, sp_b = sens_spec(d, "call_B2")
    return {"n_tests": int(len(d)), "macro_auroc_RedoxKO": macro(d, "score_B2"), "macro_auroc_full": macro(d, "score_full"),
            "d_macro_auroc": pt[0], "d_macro_auroc_ci": ci(0), "success_lower_gt_0": bool(ci(0)[0] > 0),
            "mcc_RedoxKO": float(matthews_corrcoef(d.expected, d.call_B2)),
            "mcc_full": float(matthews_corrcoef(d.expected, d.call_full)), "d_mcc": pt[1], "d_mcc_ci": ci(1),
            "sens_RedoxKO": se_b, "sens_full": se_f, "spec_RedoxKO": sp_b, "spec_full": sp_f}


prim = score(pd.read_csv(ROOT / "results/additional/x1_calls.csv"))
hold = score(pd.read_csv(ROOT / "results/additional/hybrid/holdout_calls.csv"))
prim.to_csv(OUT / "newspecies_primary_scores.csv", index=False)
hold.to_csv(OUT / "newspecies_withheld_scores.csv", index=False)
res = {"n_training_genomes": len(Xtr), "n_varying_kos": len(cols),
       "HF1_HF2_primary100": compare(prim), "HF3_withheld43": compare(hold)}
res["HF3_same_sign"] = {k: bool(np.sign(res["HF3_withheld43"][k]) == np.sign(res["HF1_HF2_primary100"][k]))
                        for k in ("d_macro_auroc", "d_mcc")}
print(json.dumps(res, indent=1, default=float), flush=True)

# ---------------------------------------------------------------- HF4: genus-blocked CV
fold = lab.fold_genus.to_numpy()


def outer_one(j, k):
    te = fold == k
    m, *_ = fit(Xtr[~te], Y[~te, j], groups[~te])
    return j, k, m.predict_proba(Xtr[te])


outs = Parallel(n_jobs=7)(delayed(outer_one)(j, k) for j in range(len(traits)) for k in np.unique(fold))
pred = {t: np.zeros(len(lab)) for t in traits}
for j, k, p in outs:
    pred[traits[j]][fold == k] = p
oof = pd.read_parquet(ROOT / "results/oof/B2.parquet")
oof = oof[oof.level == "genus"].groupby(["trait", "genome_id"]).prob.mean()
base = {t: oof.loc[t].reindex(lab.genome_id).to_numpy() for t in traits}


def macro_ap(idx, src):
    v = [safe_ap(Y[idx, j], src[t][idx]) for j, t in enumerate(traits) if Y[idx, j].sum() > 0]
    return float(np.nanmean(v))


allidx = np.arange(len(lab))
ug = np.unique(groups)
rng = np.random.default_rng(SEED)
diffs = []
for _ in range(2000):
    idx = np.concatenate([np.where(groups == u)[0] for u in rng.choice(ug, len(ug))])
    diffs.append(macro_ap(idx, pred) - macro_ap(idx, base))
res["HF4_cv_genus"] = {"macro_auprc_RedoxKO": macro_ap(allidx, base), "macro_auprc_full": macro_ap(allidx, pred),
                       "d": macro_ap(allidx, pred) - macro_ap(allidx, base),
                       "d_ci": [float(x) for x in np.percentile(diffs, [2.5, 97.5])],
                       "non_inferior_margin_-0.02": bool(np.percentile(diffs, 2.5) > -0.02),
                       "per_trait": {t: {"RedoxKO": safe_ap(Y[:, j], base[t]), "full": safe_ap(Y[:, j], pred[t])}
                                     for j, t in enumerate(traits)}}
(OUT / "results.json").write_text(json.dumps(res, indent=2, default=float))
print(json.dumps(res["HF4_cv_genus"], indent=1, default=float))
