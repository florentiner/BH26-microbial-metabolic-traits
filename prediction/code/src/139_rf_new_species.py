"""Random forest on the 852-KO panel tested on the 100 primary new species (255 growth experiments), next to RedoxKO.

The forest is fitted on all 358 labeled genomes with min_samples_leaf and the MCC threshold chosen by inner grouped CV
(as for RedoxKO). Macro AUROC over categories and pooled MCC, with paired 95% bootstrap intervals over papers.
usage: python src/139_rf_new_species.py
Output: results/additional/rf_new_species.csv
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import best_mcc_threshold, inner_folds, load_features, load_labelled, safe_ap  # noqa: E402

tx = importlib.import_module("47_trees_xai")
tx.TREE_JOBS = 8
V = ROOT / "validation_new_organisms"


def main():
    lab, traits, Y = load_labelled()
    X = load_features("energy_panel", lab.genome_id)
    keep = X.columns[X.std(axis=0) > 0]
    Xi = X[keep].to_numpy()
    d = pd.read_csv(V / "scored178_B2_balanced100.csv")
    d = d[~d.category.isin(["donor_other", "acceptor_other"])].copy()
    Xv = (pd.read_parquet(V / "ko_panel_validation.parquet") > 0).astype(np.float32)
    Xv = pd.concat([Xv, (pd.read_parquet(V / "ko_panel_auto.parquet") > 0).astype(np.float32)])
    Xv = Xv[~Xv.index.duplicated()].reindex(columns=keep).fillna(0)
    groups = lab.group_genus.to_numpy()
    sc, th = {}, {}
    for t in sorted(d.category.unique()):
        j = traits.index(t)
        y = Y[:, j]
        m = ~np.isnan(y)
        folds = inner_folds(y[m], groups[m])
        best = None
        for cfg in tx.CONFIGS["rf"]:
            oof = np.zeros(m.sum())
            for itr, iva in folds:
                oof[iva] = tx.fit_predict(cfg, Xi[m][itr], y[m][itr], Xi[m][iva])
            ap = safe_ap(y[m], oof)
            if best is None or ap > best[0]:
                best = (ap, cfg, oof)
        _, cfg, oof = best
        th[t] = best_mcc_threshold(y[m], oof)
        sc[t] = pd.Series(tx.fit_predict(cfg, Xi[m], y[m], Xv.to_numpy()), index=Xv.index)
        print(t, cfg, flush=True)
    d["score_rf"] = [sc[t].get(a, np.nan) for a, t in zip(d.accession, d.category)]
    d["call_rf"] = [int(s >= th[t]) for s, t in zip(d.score_rf, d.category)]

    def stats(x, col, call):
        au = [roc_auc_score(g.expected, g[col]) for _, g in x.groupby("category") if 0 < g.expected.sum() < len(g)]
        return np.mean(au), matthews_corrcoef(x.expected, x[call])
    rf, b2 = stats(d, "score_rf", "call_rf"), stats(d, "score", "call")
    rng = np.random.default_rng(20261021)
    papers = d.paper.unique()
    idx = {p: np.where(d.paper.to_numpy() == p)[0] for p in papers}
    diffs = []
    for _ in range(2000):
        x = d.iloc[np.concatenate([idx[p] for p in rng.choice(papers, len(papers))])]
        try:
            a, b = stats(x, "score_rf", "call_rf"), stats(x, "score", "call")
            diffs.append((a[0] - b[0], a[1] - b[1]))
        except ValueError:
            pass
    diffs = np.array(diffs)
    out = pd.DataFrame([dict(model="RedoxKO", macro_auroc=b2[0], mcc=b2[1]),
                        dict(model="Random forest", macro_auroc=rf[0], mcc=rf[1],
                             auroc_diff=rf[0] - b2[0], auroc_lo=np.percentile(diffs[:, 0], 2.5),
                             auroc_hi=np.percentile(diffs[:, 0], 97.5), mcc_diff=rf[1] - b2[1],
                             mcc_lo=np.percentile(diffs[:, 1], 2.5), mcc_hi=np.percentile(diffs[:, 1], 97.5))])
    out.to_csv(ROOT / "results/additional/rf_new_species.csv", index=False, float_format="%.4f")
    d.to_csv(ROOT / "results/additional/rf_new_species_scores.csv", index=False)
    print(out.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
