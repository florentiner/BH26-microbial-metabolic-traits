"""D (docs/preregistration_label_recuration_D.md): refit B2 and B2XLE_large on the re-curated labels
(data/processed/labels_curated_D.tsv; NaN = unknown, dropped per trait) with exactly the selection/fit functions of
src/09_train_final.py, save to results/models_curated_D/ (released models untouched), and compare with the original
released models on data whose labels did not change: the 43 held-out new species (primary), the external 562 and
the primary 100 (reported only).

usage: python src/102_train_eval_curated.py fit | eval
"""
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import load_features, load_labelled  # noqa: E402

tf = __import__("09_train_final")
MOD = ROOT / "results/models_curated_D"
OUT = ROOT / "results/additional/label_recuration"
VAL = ROOT / "validation_new_organisms"
SEED = 20261015
CATCH = ["donor_other", "acceptor_other"]


def _fit():
    lab, traits, Y = load_labelled(label_mode="curated_D")
    groups = lab.group_genus.to_numpy()
    b2 = load_features("energy_panel", lab.genome_id)
    b2 = b2.loc[:, b2.std(axis=0) > 0]
    emb = load_features("bacformer_large_mean", lab.genome_id)
    arms = tf.arms_mod
    for arm in ["B2", "B2XLE_large"]:
        bundle = {}
        for j, t in enumerate(traits):
            m = ~np.isnan(Y[:, j])
            y = Y[m, j].astype(int)
            if arm == "B2":
                model, thr, cfg, inner, cal = tf.select_and_fit(b2.loc[m], y, groups[m], arms.ko_configs(),
                                                                lambda cfg, _c=list(b2.columns): _c)
            else:
                model, thr, cfg, inner, cal = tf.select_and_fit_xle_large(b2.loc[m], emb.loc[m], y, groups[m])
            bundle[t] = {"type": "estimator", "model": model, "threshold": thr, "config": cfg,
                         "inner_cv_auprc_all_labelled": inner, "calibrator": cal, "n_labelled": int(m.sum())}
            print(arm, t, "labelled", int(m.sum()), "pos", int(y.sum()), cfg, f"thr={thr:.3f}", flush=True)
        (MOD / arm).mkdir(parents=True, exist_ok=True)
        joblib.dump(bundle, MOD / arm / "model.joblib")


def predict(bundle, X):
    return {t: np.asarray(b["model"].predict_proba(X)) for t, b in bundle.items()}


def evaluate():
    OUT.mkdir(parents=True, exist_ok=True)
    lab, traits, _ = load_labelled()
    keep = load_features("energy_panel", lab.genome_id)
    keep = keep.columns[keep.std(axis=0) > 0]
    # new species (blind), features
    ko_v = pd.read_parquet(VAL / "ko_panel_validation.parquet")
    ko_v = (ko_v > 0).astype(np.float32).reindex(columns=keep).fillna(0)
    emb_v = pd.read_parquet(VAL / "embeddings_all.parquet")
    # external
    ext = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t")
    ko_e = pd.read_csv(ROOT / "data/external/ko_matrix_panel.tsv", sep="\t", index_col=0).reindex(ext.genome_id)[keep]
    emb_e = pd.read_parquet(ROOT / "data/kaggle/large_external/genome_features_bacformer_large_mean.parquet")
    blind = pd.read_csv(VAL / "scored178_B2.csv")
    blind = blind[(blind.blind == True) & ~blind.category.isin(CATCH)]  # noqa: E712
    primary = set((VAL / "balanced_prospective_primary_100.txt").read_text().split())
    rows, frames = [], {}
    for arm in ["B2", "B2XLE_large"]:
        orig = joblib.load(ROOT / f"results/models/{arm}/model.joblib")
        cur = joblib.load(MOD / arm / "model.joblib")
        Xv = ko_v if arm == "B2" else ko_v.join(emb_v, how="inner")
        Xe = ko_e if arm == "B2" else ko_e.join(emb_e, how="inner")
        for name, model in [("original", orig), ("curated_D", cur)]:
            pv, pe = predict(model, Xv), predict(model, Xe)
            # new species: per call
            d = blind.copy()
            idx = {a: i for i, a in enumerate(Xv.index)}
            d["score"] = [pv[c][idx[a]] for a, c in zip(d.accession, d.category)]
            d["callm"] = [int(s >= model[c]["threshold"]) for s, c in zip(d.score, d.category)]
            for sname, dd in [("held-out 43", d[~d.accession.isin(primary)]), ("primary 100", d[d.accession.isin(primary)])]:
                rows.append(dict(model=arm, labels=name, test_set=sname, **scores_calls(dd)))
                frames[(arm, name, sname)] = dd.reset_index(drop=True)
            # external: per genome x trait
            e = []
            eidx = {g: i for i, g in enumerate(Xe.index)}
            for t in [t for t in traits if t in ext.columns]:
                for g, y in zip(ext.genome_id, ext[t]):
                    if pd.notna(y) and g in eidx:
                        e.append((g, t, int(y), pe[t][eidx[g]], int(pe[t][eidx[g]] >= model[t]["threshold"])))
            E = pd.DataFrame(e, columns=["genome", "category", "expected", "score", "callm"])
            E["paper"] = E.genome.map(ext.set_index("genome_id").genus)
            rows.append(dict(model=arm, labels=name, test_set="external 562", **scores_calls(E)))
            frames[(arm, name, "external 562")] = E.reset_index(drop=True)
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "curated_vs_original.csv", index=False)
    deltas = []
    for arm in ["B2", "B2XLE_large"]:
        for sname in ["held-out 43", "external 562", "primary 100"]:
            a, b = frames[(arm, "original", sname)], frames[(arm, "curated_D", sname)]
            rng = np.random.default_rng(SEED)
            ids = a.paper.unique()
            ix = {i: np.where(a.paper.to_numpy() == i)[0] for i in ids}
            v = {"mcc": [], "macro_auroc": [], "sensitivity": [], "specificity": []}
            for _ in range(1000):
                take = np.concatenate([ix[i] for i in rng.choice(ids, len(ids))])
                try:
                    sa, sb = scores_calls(a.iloc[take]), scores_calls(b.iloc[take])
                except ValueError:
                    continue
                for k in v:
                    v[k].append(sb[k] - sa[k])
            sa, sb = scores_calls(a), scores_calls(b)
            deltas.append(dict(model=arm, test_set=sname, **{f"{k}_change": sb[k] - sa[k] for k in v},
                               **{f"{k}_ci": np.nanpercentile(v[k], [2.5, 97.5]).round(3).tolist() for k in v}))
    Dl = pd.DataFrame(deltas)
    Dl.to_csv(OUT / "curated_minus_original.csv", index=False)
    print(Dl.round(3).to_string(index=False))
    pd.set_option("display.width", 250)
    print(R.round(3).to_string(index=False))


def scores_calls(d):
    y, c, s = d.expected.astype(int), d.callm.astype(int), d.score
    tp, tn = ((y == 1) & (c == 1)).sum(), ((y == 0) & (c == 0)).sum()
    aucs = [roc_auc_score(g.expected, g.score) for _, g in d.groupby("category") if g.expected.nunique() == 2]
    aps = [average_precision_score(g.expected, g.score) for _, g in d.groupby("category") if g.expected.nunique() == 2]
    return dict(calls=len(d), accuracy=(tp + tn) / len(d), mcc=matthews_corrcoef(y, c),
                sensitivity=tp / max((y == 1).sum(), 1), specificity=tn / max((y == 0).sum(), 1),
                macro_auroc=float(np.mean(aucs)), macro_auprc=float(np.mean(aps)))


if __name__ == "__main__":
    {"fit": _fit, "eval": evaluate}[sys.argv[1]]()
