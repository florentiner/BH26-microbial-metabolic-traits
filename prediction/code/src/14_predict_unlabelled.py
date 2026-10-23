"""Apply a saved model to the 228 unlabelled ("not in database") genomes using precomputed features.

Equivalent to src/predict.py (same features, same saved models) without re-downloading/annotating.
usage: python src/14_predict_unlabelled.py [--arm released|B2|B1r|...]
Output: results/predictions_unlabelled_<arm>.csv (long) and results/predictions_unlabelled_<arm>_calls.csv (wide 0/1)
"""
import argparse
import json

import joblib
import numpy as np
import pandas as pd

from redox_common import PROC, RES, RULE_GROUPS, apply_calibrator, load_features

FEATURES = {"B1r": "markers_recomputed", "B2": "energy_panel", "B2XL": "energy_panel", "B2cnt": "energy_panel_log1p_counts",
            "E": "bacformer_mean", "E_MAG": "bacformer_mag_mean"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="released")
    args = ap.parse_args()
    arm = json.loads((RES / "models/RELEASED.json").read_text())["released_arm"] if args.arm == "released" else args.arm
    labels = pd.read_csv(PROC / "labels.tsv", sep="\t")
    unl = labels[~labels.labelled].reset_index(drop=True)
    train_genera = set(labels.loc[labels.labelled, "genus"])
    train_families = set(labels.loc[labels.labelled, "family"].dropna())
    bundle = joblib.load(RES / "models" / arm / "model.joblib")
    if arm == "B2E":
        X = pd.concat([load_features("energy_panel", unl.genome_id), load_features("bacformer_mean", unl.genome_id)],
                      axis=1)
    elif arm == "C":
        X = pd.concat([load_features("esm2_mean", unl.genome_id), load_features("esm2_max", unl.genome_id)], axis=1)
    elif arm == "A":
        X = load_features("markers_recomputed", unl.genome_id)
    else:
        X = load_features(FEATURES[arm], unl.genome_id)
    rows = []
    for t, b in bundle.items():
        if b["type"] == "rule":
            g = RULE_GROUPS[t]
            score = np.max([X[k].mean(axis=1).to_numpy() for k in g.values()], axis=0) if g else np.zeros(len(X))
        else:
            score = b["model"].predict_proba(X)
        cal_score = apply_calibrator(b.get("calibrator"), score) if b["type"] == "estimator" else score
        for i, s in enumerate(score):
            rows.append(dict(genome_id=unl.genome_id[i], organism=unl.organism_name[i], phylum=unl.phylum[i],
                             genus=unl.genus[i], trait=t, score=float(s), calibrated_score=float(cal_score[i]),
                             threshold=b["threshold"],
                             predicted=int(s >= b["threshold"]), arm=arm,
                             genus_in_training_set=unl.genus[i] in train_genera,
                             family_in_training_set=unl.family[i] in train_families))
    out = pd.DataFrame(rows)
    # per-trait reliability from the outer CV. Genomes from a family absent from training use family-level hold-out
    # metrics, the others genus-level. A call is reliable when its thresholded calls agree at least moderately with
    # labels on held-out clades (MCC 95% cluster-bootstrap CI lower bound > 0.2; MCC is prevalence-aware). Whether
    # the model beats the taxonomic null and the external-validation AUPRC are separate columns. Scores always kept.
    mt = pd.read_csv(RES / "metrics_per_trait.csv")
    mt = mt[mt.arm == arm].set_index(["split_level", "trait"])
    level = np.where(out.family_in_training_set, "genus", "family")
    key = list(zip(level, out.trait))
    out["cv_level_used"] = level
    out["cv_auprc"] = [mt.auprc.get(k, np.nan) for k in key]
    out["cv_mcc"] = [mt.mcc.get(k, np.nan) for k in key]
    out["cv_beats_taxonomic_null_q05"] = [bool(mt.delta_vs_null_q_bh.get(k, 1) < 0.05 and mt.delta_vs_null.get(k, 0) > 0)
                                          for k in key]
    # audit round 3 (M2): a call is only emitted when the model shows out-of-clade skill *beyond the taxonomic
    # null*, not just non-random agreement; external evidence (where the trait is externally evaluable) must not
    # contradict it (AUPRC lift > 1).
    out["call_reliable"] = [bool(mt.mcc_ci_low.get(k, -1) > 0.2 and mt.delta_vs_null_q_bh.get(k, 1) < 0.05
                                 and mt.delta_vs_null.get(k, 0) > 0) for k in key]
    ext_file = RES / "external_validation_per_trait.csv"
    if ext_file.exists():
        ev = pd.read_csv(ext_file)
        ev = ev[(ev.subset == "all") & (ev.arm == arm) & (ev.evaluable == True)].set_index("trait")  # noqa: E712
        out["external_auprc"] = out.trait.map(ev.auprc)
        out["external_auprc_lift"] = out.trait.map(ev.auprc_lift)
        contradicted = out.external_auprc_lift.notna() & (out.external_auprc_lift <= 1)
        out.loc[contradicted, "call_reliable"] = False
    out.loc[~out.call_reliable.astype(bool), "predicted"] = np.nan
    out.to_csv(RES / f"predictions_unlabelled_{arm}.csv", index=False, float_format="%.4f")
    wide = out.pivot(index=["genome_id", "organism", "genus_in_training_set"], columns="trait", values="predicted")
    print("traits with reliable calls:", sorted(out.loc[out.call_reliable.astype(bool), "trait"].unique()))
    wide.to_csv(RES / f"predictions_unlabelled_{arm}_calls.csv")
    print(len(unl), "genomes;", "genus seen in training:", int(unl.genus.isin(train_genera).sum()))
    print(wide.sum().sort_values(ascending=False).to_string())


if __name__ == "__main__":
    main()
