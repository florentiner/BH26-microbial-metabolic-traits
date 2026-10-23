"""Refit the released pathway arm on all labelled pathway genomes and save the deployable model.

Released arm = highest genus-level macro AUPRC in results/pathways/metrics_macro.csv among the learned arms
(PB1/PB2/PB3/PE); the rule arm PA is a baseline, never a release candidate. Per trait, config and decision
threshold are chosen by grouped (genus) inner CV over all labelled genomes, then the model is refitted on all of
them. Performance claims come from the outer CV, never from this refit.

A trait is flagged reliable for reporting on new genomes when, at family level (the harder split): the MCC
bootstrap CI lower bound > 0.2 and it beats the taxonomic null at q < 0.05 — the same gate as the
donor/acceptor release (docs/deviations.md #16).

Output: results/models/pathways/<arm>/model.joblib, metadata.json, results/models/pathways/RELEASED.json
usage: python src/35_pathway_train_final.py [--arm PB2]
"""
import argparse
import json
from datetime import date

import joblib
import numpy as np
import pandas as pd

from pathway_common import PATH, RES, energy_panel_subset, ko_matrix, load_pathways, module_completeness
from redox_common import FittedModel, best_mcc_threshold, fit_calibrator, fit_models, inner_folds, safe_ap, sha256

LR_C = [0.01, 0.1, 1.0, 10.0]
CONFIGS = [{"model": "logreg", "C": c} for c in LR_C] + [{"model": "xgb", "max_depth": d} for d in (2, 4)]
FEATURES = {"PB1": "bh26_energy_panel_852ko", "PB2": "bh26_full_ko_profile",
            "PB3": "bh26_full_ko_profile + kegg_module_completeness"}


def features(arm, genome_ids):
    if arm == "PB1":
        return energy_panel_subset(genome_ids)
    if arm == "PB2":
        return ko_matrix(genome_ids)
    if arm == "PB3":
        return pd.concat([ko_matrix(genome_ids), module_completeness(genome_ids).add_prefix("mod_")], axis=1)
    raise SystemExit(f"arm {arm} is not a release candidate")


def select_and_fit(Xdf, y, groups):
    folds = inner_folds(y, groups)
    X = Xdf.to_numpy(dtype=np.float32)
    scores, oofs = {}, {}
    for ci, cfg in enumerate(CONFIGS):
        oof = np.zeros(len(y))
        for itr, iva in folds:
            oof[iva] = np.mean([m.predict_proba(X[iva])[:, 1] for m in fit_models(cfg, X[itr], y[itr])], axis=0)
        scores[ci], oofs[ci] = safe_ap(y, oof), oof
    bi = max(scores, key=lambda c: (np.nan_to_num(scores[c], nan=-1), -c))
    cfg = CONFIGS[bi]
    model = FittedModel(list(Xdf.columns), cfg, fit_models(cfg, X, y))
    return model, best_mcc_threshold(y, oofs[bi]), cfg, scores[bi], fit_calibrator(y, oofs[bi])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default=None, help="override the arm selected from the metrics table")
    args = ap.parse_args()
    lab, traits, Y = load_pathways()
    groups = lab.group_genus.to_numpy()
    macro = pd.read_csv(RES / "pathways/metrics_macro.csv")
    g = macro[(macro.split_level == "genus") & macro.arm.isin(FEATURES)].sort_values("macro_auprc", ascending=False)
    arm = args.arm or g.arm.iloc[0]
    print(g[["arm", "macro_auprc", "macro_auprc_ci_low", "macro_auprc_ci_high"]].to_string(index=False))
    print("released arm:", arm, flush=True)

    Xdf = features(arm, lab.genome_id)
    Xdf = Xdf.loc[:, Xdf.std(axis=0) > 0]
    out = RES / "models/pathways" / arm
    out.mkdir(parents=True, exist_ok=True)
    bundle = {}
    for j, t in enumerate(traits):
        model, thr, cfg, inner, cal = select_and_fit(Xdf, Y[:, j], groups)
        bundle[t] = {"type": "estimator", "model": model, "threshold": thr, "config": cfg,
                     "inner_cv_auprc_all_labelled": inner, "calibrator": cal}
        print(arm, t, cfg, f"thr={thr:.3f}", f"inner AUPRC={inner:.3f}", flush=True)
    joblib.dump(bundle, out / "model.joblib")

    pt = pd.read_csv(RES / "pathways/metrics_per_trait.csv")
    fam = pt[(pt.arm == arm) & (pt.split_level == "family")].set_index("trait")
    gen = pt[(pt.arm == arm) & (pt.split_level == "genus")].set_index("trait")
    rel = {t: bool(fam.at[t, "mcc_ci_low"] > 0.2 and fam.at[t, "delta_vs_null_q_bh"] < 0.05
                   and fam.at[t, "delta_vs_null"] > 0) for t in traits if t in fam.index}
    meta = {"arm": arm, "task": "energy metabolism pathways", "date": str(date.today()),
            "trained_on": f"{len(lab)} labelled genomes (all)", "traits": traits, "feature_set": FEATURES[arm],
            "feature_source": "RDF Portal BH26 KofamScan KO profiles (the store exposed through TogoMCP)",
            "labels": "Madin et al. 2020 `pathways` field, merged per docs/preregistration_pathways.md",
            "labels_sha256": sha256(PATH / "labels.tsv"),
            "reliable_traits": sorted(t for t, v in rel.items() if v),
            "auprc_genus": {t: float(gen.at[t, "auprc"]) for t in traits if t in gen.index},
            "auprc_family": {t: float(fam.at[t, "auprc"]) for t in traits if t in fam.index},
            "cv_performance_file": "results/pathways/metrics_per_trait.csv (outer CV; the refit has no held-out "
                                   "estimate)"}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2))
    (RES / "models/pathways/RELEASED.json").write_text(json.dumps(
        {"arm": arm, "path": f"results/models/pathways/{arm}/model.joblib",
         "selected_by": "highest genus-level macro AUPRC among PB1/PB2/PB3/PE",
         "n_reliable_traits": int(sum(rel.values())), "n_traits": len(traits)}, indent=2))
    print(f"reliable traits: {sorted(t for t, v in rel.items() if v)}")


if __name__ == "__main__":
    main()
