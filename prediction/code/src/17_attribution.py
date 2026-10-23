"""Feature attribution of the released KO model vs known marker genes (protocol stopping rule: if top features are
housekeeping/lineage genes the model reads clade, not metabolism).

For every trait, a model is refit on each outer training fold (genus hold-out) with the config selected there
(from results/oof/<arm>.parquet), and permutation importance is measured on that fold's *test* genomes (AUPRC drop,
10 repeats). Importances are averaged over folds. Marker enrichment: rank of the trait's rule-marker KOs
(RULE_GROUPS) among all features; AUPRC-based importance is label-dependent, so it is computed only on held-out
genomes.
Output: results/attribution_top_features.csv, results/attribution_marker_enrichment.csv
"""
import argparse
import json

import numpy as np
import pandas as pd
from sklearn.inspection import permutation_importance
from sklearn.metrics import average_precision_score

from redox_common import RES, RULE_GROUPS, fit_models, load_features, load_labelled

TOP = 15


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="B2")
    args = ap.parse_args()
    lab, traits, Y = load_labelled()
    feats = {"B2": "energy_panel", "B1r": "markers_recomputed"}[args.arm]
    Xdf = load_features(feats, lab.genome_id)
    Xdf = Xdf.loc[:, Xdf.std(axis=0) > 0]
    X, cols = Xdf.to_numpy(), np.array(Xdf.columns)
    panel = pd.read_csv(RES.parent / "data/processed/ko_panel.tsv", sep="\t").set_index("KO").definition
    oof = pd.read_parquet(RES / "oof" / f"{args.arm}.parquet")
    oof = oof[oof.level == "genus"]
    fold = lab.fold_genus.to_numpy()
    top_rows, enr_rows = [], []
    for j, t in enumerate(traits):
        imp = np.zeros(len(cols))
        n_used = 0
        for k in sorted(np.unique(fold)):
            tr, te = fold != k, fold == k
            if Y[te, j].sum() == 0:
                continue
            cfg = json.loads(oof[(oof.trait == t) & (oof.fold == k)].config.iloc[0])
            models = fit_models(cfg, X[tr], Y[tr, j])

            class Avg:  # sklearn-compatible averaged ensemble for permutation_importance
                def fit(self, *a):
                    return self

                def predict_proba(self, Z):
                    p = np.mean([m.predict_proba(Z)[:, 1] for m in models], axis=0)
                    return np.c_[1 - p, p]

            r = permutation_importance(Avg(), X[te], Y[te, j], n_repeats=10, random_state=0, n_jobs=6,
                                       scoring=lambda est, Z, y: average_precision_score(y, est.predict_proba(Z)[:, 1]))
            imp += r.importances_mean
            n_used += 1
        imp /= max(n_used, 1)
        order = np.argsort(-imp)
        markers = {ko for kos in RULE_GROUPS[t].values() for ko in kos}
        for rank, i in enumerate(order[:TOP], 1):
            top_rows.append(dict(trait=t, rank=rank, KO=cols[i], definition=panel.get(cols[i], ""),
                                 importance_auprc_drop=imp[i], is_rule_marker=cols[i] in markers))
        present = [m for m in markers if m in set(cols)]
        ranks = [int(np.where(cols[order] == m)[0][0]) + 1 for m in present]
        enr_rows.append(dict(trait=t, n_rule_markers_in_features=len(present),
                             n_rule_markers_in_top15=sum(r <= TOP for r in ranks),
                             best_marker_rank=min(ranks) if ranks else np.nan,
                             n_features=len(cols), n_folds_used=n_used))
        print(t, "top:", ", ".join(cols[order[:5]]), "| markers in top15:", enr_rows[-1]["n_rule_markers_in_top15"],
              flush=True)
    pd.DataFrame(top_rows).to_csv(RES / "attribution_top_features.csv", index=False, float_format="%.4f")
    pd.DataFrame(enr_rows).to_csv(RES / "attribution_marker_enrichment.csv", index=False)


if __name__ == "__main__":
    main()
