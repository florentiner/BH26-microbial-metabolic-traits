"""A+B (docs/preregistration_hybrid_AB.md): internal CV comparison B2M vs B2, final B2M fit on all 358 genomes
(results/models_hybrid/B2M; released models untouched), and evaluation of B2, B2M, the decision hybrid B2M+M and
METABOLIC on the 43 held-out new species (primary), the external 562 and the primary 100 (reported only).

usage: python src/104_hybrid_eval.py cv | fit | eval
"""
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import RES, load_features, load_labelled, oof_matrix, paired_macro_delta, safe_ap  # noqa: E402

tf = __import__("09_train_final")
m95 = __import__("95_metabolic_score")
OUT = ROOT / "results/additional/hybrid"
MOD = ROOT / "results/models_hybrid"
VAL = ROOT / "validation_new_organisms"
SEED = 20261015
CATCH = ["donor_other", "acceptor_other"]
RELEASED = ["acceptor_oxygen", "acceptor_sulfate", "acceptor_sulfur", "donor_sulfur_oxidation",
            "acceptor_CO2_methanogen_acetogen"]


def cv():
    OUT.mkdir(parents=True, exist_ok=True)
    lab, traits, Y = load_labelled()
    rows = []
    for level in ["genus", "family"]:
        Pa = oof_matrix(RES / "oof", "B2M", level, lab, traits)
        Pb = oof_matrix(RES / "oof", "B2", level, lab, traits)
        ma = np.nanmean([safe_ap(Y[:, j], Pa[:, j]) for j in range(len(traits))])
        mb = np.nanmean([safe_ap(Y[:, j], Pb[:, j]) for j in range(len(traits))])
        d = paired_macro_delta(lab, traits, Y, Pa, Pb, level)
        rows.append(dict(level=level, B2M=ma, B2=mb, delta=ma - mb, ci=str(d)))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "cv_B2M_vs_B2.csv", index=False)
    print(R.to_string(index=False))


def features_train(lab):
    b2 = load_features("energy_panel", lab.genome_id)
    mh = load_features("metabolic_hmm", lab.genome_id)
    X = b2.join(mh)
    return X.loc[:, X.std(axis=0) > 0]


def fit():
    lab, traits, Y = load_labelled()
    X = features_train(lab)
    groups = lab.group_genus.to_numpy()
    bundle = {}
    for j, t in enumerate(traits):
        model, thr, cfg, inner, cal = tf.select_and_fit(X, Y[:, j].astype(int), groups, tf.arms_mod.ko_configs(),
                                                        lambda cfg, _c=list(X.columns): _c)
        bundle[t] = {"type": "estimator", "model": model, "threshold": thr, "config": cfg,
                     "inner_cv_auprc_all_labelled": inner, "calibrator": cal}
        print("B2M", t, cfg, f"thr={thr:.3f}", flush=True)
    (MOD / "B2M").mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, MOD / "B2M" / "model.joblib")


def stats(d, col):
    y, c = d.expected.astype(int), d[col].astype(int)
    tp, tn = ((y == 1) & (c == 1)).sum(), ((y == 0) & (c == 0)).sum()
    return dict(calls=len(d), accuracy=(tp + tn) / len(d), mcc=matthews_corrcoef(y, c),
                sensitivity=tp / max((y == 1).sum(), 1), specificity=tn / max((y == 0).sum(), 1))


def boot(d, a, b, cluster="paper", n=1000):
    rng = np.random.default_rng(SEED)
    ids = d[cluster].unique()
    ix = {i: np.where(d[cluster].to_numpy() == i)[0] for i in ids}
    v = []
    for _ in range(n):
        dd = d.iloc[np.concatenate([ix[i] for i in rng.choice(ids, len(ids))])]
        v.append(matthews_corrcoef(dd.expected, dd[b]) - matthews_corrcoef(dd.expected, dd[a]))
    return np.percentile(v, [2.5, 97.5]).round(3).tolist()


def evaluate():
    OUT.mkdir(parents=True, exist_ok=True)
    lab, traits, _ = load_labelled()
    cols = features_train(lab).columns
    keep_b2 = [c for c in cols if not c.startswith("metabolic_")]
    b2m, b2 = joblib.load(MOD / "B2M/model.joblib"), joblib.load(ROOT / "results/models/B2/model.joblib")
    # new species features
    ko_v = (pd.read_parquet(VAL / "ko_panel_validation.parquet") > 0).astype(np.float32)
    mh_v = pd.read_parquet(ROOT / "data/additional/metabolic_hmm/newspecies.parquet")
    Xv = ko_v.reindex(columns=keep_b2).fillna(0).join(mh_v, how="inner").reindex(columns=cols).fillna(0)
    blind = pd.read_csv(VAL / "scored178_B2.csv")
    blind = blind[(blind.blind == True) & ~blind.category.isin(CATCH)].copy()  # noqa: E712
    idx = {a: i for i, a in enumerate(Xv.index)}
    pv = {t: np.asarray(b2m[t]["model"].predict_proba(Xv)) for t in traits}
    blind["score_B2M"] = [pv[c][idx[a]] for a, c in zip(blind.accession, blind.category)]
    blind["call_B2M"] = [int(s >= b2m[c]["threshold"]) for s, c in zip(blind.score_B2M, blind.category)]
    blind = blind.rename(columns={"call": "call_B2", "score": "score_B2"})
    primary = set((VAL / "balanced_prospective_primary_100.txt").read_text().split())
    hold = blind[~blind.accession.isin(primary)].copy()
    # real METABOLIC on the held-out 43 (decision hybrid and comparator)
    m95.W2 = ROOT / "data/additional/metabolic/holdout_out/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet2.tsv"
    m95.OUT = OUT
    hold[["accession", "category"]].assign(expected=hold.expected).to_csv(OUT / "x1_calls.csv", index=False)
    m95.main()
    met = pd.read_csv(OUT / "x1_metabolic_calls.csv")
    hold = hold.merge(met, on=["accession", "category"], how="left")
    hold["call_hybrid"] = np.where(hold.category.isin(RELEASED), hold.call_B2M, hold.call_METABOLIC)
    rows = []
    for col, name in [("call_B2", "B2"), ("call_B2M", "B2M"), ("call_hybrid", "B2M+M (decision hybrid)"),
                      ("call_METABOLIC", "METABOLIC")]:
        r = dict(test_set="held-out 43", model=name, **stats(hold, col))
        if name != "B2":
            r["mcc_minus_B2"] = stats(hold, col)["mcc"] - stats(hold, "call_B2")["mcc"]
            r["ci_vs_B2"] = boot(hold, "call_B2", col)
        if name != "METABOLIC":
            r["mcc_minus_METABOLIC"] = stats(hold, col)["mcc"] - stats(hold, "call_METABOLIC")["mcc"]
            r["ci_vs_METABOLIC"] = boot(hold, "call_METABOLIC", col)
        rows.append(r)
    prim = blind[blind.accession.isin(primary)]
    for col, name in [("call_B2", "B2"), ("call_B2M", "B2M")]:
        rows.append(dict(test_set="primary 100 (not independent)", model=name, **stats(prim, col)))
    hold.to_csv(OUT / "holdout_calls.csv", index=False)
    # external 562
    ext = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t")
    ko_e = pd.read_csv(ROOT / "data/external/ko_matrix_panel.tsv", sep="\t", index_col=0)
    mh_e = pd.read_parquet(ROOT / "data/additional/metabolic_hmm/external.parquet")
    Xe = ko_e.reindex(columns=keep_b2).fillna(0).join(mh_e, how="inner").reindex(columns=cols).fillna(0)
    Xe2 = ko_e.reindex(columns=keep_b2).loc[Xe.index]
    e = []
    for t in [t for t in traits if t in ext.columns]:
        pm = np.asarray(b2m[t]["model"].predict_proba(Xe))
        pb = np.asarray(b2[t]["model"].predict_proba(Xe2))
        y = ext.set_index("genome_id").reindex(Xe.index)[t]
        for g, yy, a, b in zip(Xe.index, y, pm, pb):
            if pd.notna(yy):
                e.append((g, t, int(yy), b, a, int(b >= b2[t]["threshold"]), int(a >= b2m[t]["threshold"])))
    E = pd.DataFrame(e, columns=["genome", "category", "expected", "score_B2", "score_B2M", "call_B2", "call_B2M"])
    E["paper"] = E.genome.map(ext.set_index("genome_id").genus)
    for col, sc, name in [("call_B2", "score_B2", "B2"), ("call_B2M", "score_B2M", "B2M")]:
        au = np.mean([roc_auc_score(g.expected, g[sc]) for _, g in E.groupby("category") if g.expected.nunique() == 2])
        ap = np.mean([average_precision_score(g.expected, g[sc]) for _, g in E.groupby("category")
                      if g.expected.nunique() == 2])
        rows.append(dict(test_set="external 562", model=name, macro_auroc=au, macro_auprc=ap, **stats(E, col)))
    rows[-1]["ci_vs_B2"] = boot(E, "call_B2", "call_B2M")
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "hybrid_results.csv", index=False)
    pd.set_option("display.width", 250)
    print(R.round(3).to_string(index=False))


if __name__ == "__main__":
    {"cv": cv, "fit": fit, "eval": evaluate}[sys.argv[1]]()
