"""Predictions of the final model RedoxKO (random forest on the 852-KO panel, results/models/RF/model.joblib) for the
228 unlabelled ("not in database") genomes: src/14_predict_unlabelled.py with the forest arm.

14_predict_unlabelled.py cannot take --arm RF as it stands (no RF entry in its FEATURES map, and it reads the
reliability gate from results/metrics_per_trait.csv and results/external_validation_per_trait.csv, which have no
forest rows). This wrapper applies the identical rules with the forest's own numbers:
  * features: energy_panel (852-KO panel), as for B2;
  * cv columns and the call gate (MCC 95% CI low > 0.2, beats the taxonomic baseline at BH q < 0.05, delta > 0; genus
    level when the genome's family is in training, family level otherwise) from the forest's outer-CV metrics,
    results/trees_xai/da/metrics_per_trait.csv (arm RF; q column q_vs_null);
  * external AUPRC and lift from the forest's external scores (results/external_validation_missing_rows_predictions.csv,
    arm "Random forest"; category evaluable with >= 10 positives and >= 10 negatives, as src/143_final_rf.py), lift <= 1
    withdraws the call.
The resulting reliable categories are checked against results/models/RF/metadata.json (reliable_traits_by_level).
usage: python src/162_predict_unlabelled_rf.py
Output: results/predictions_unlabelled_RF.csv (long), results/predictions_unlabelled_RF_calls.csv (wide 0/1)
"""
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import PROC, RES, apply_calibrator, load_features, safe_ap  # noqa: E402

ARM = "RF"


def external_lift():
    ev = pd.read_csv(RES / "external_validation_missing_rows_predictions.csv")
    ev = ev[ev.arm == "Random forest"]
    ext = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t").set_index("genome_id")
    out = {}
    for t, g in ev.groupby("trait"):
        if t not in ext.columns:
            continue
        y = ext[t].reindex(g.genome_id).to_numpy(dtype=float)
        m = ~np.isnan(y)
        if y[m].sum() >= 10 and (y[m] == 0).sum() >= 10:
            ap = safe_ap(y[m], g.score.to_numpy()[m])
            out[t] = (ap, ap / y[m].mean())
    return out


def main():
    labels = pd.read_csv(PROC / "labels.tsv", sep="\t")
    unl = labels[~labels.labelled].reset_index(drop=True)
    train_genera = set(labels.loc[labels.labelled, "genus"])
    train_families = set(labels.loc[labels.labelled, "family"].dropna())
    bundle = joblib.load(RES / "models" / ARM / "model.joblib")
    X = load_features("energy_panel", unl.genome_id)
    rows = []
    for t, b in bundle.items():
        score = b["model"].predict_proba(X)
        cal_score = apply_calibrator(b.get("calibrator"), score)
        for i, s in enumerate(score):
            rows.append(dict(genome_id=unl.genome_id[i], organism=unl.organism_name[i], phylum=unl.phylum[i],
                             genus=unl.genus[i], trait=t, score=float(s), calibrated_score=float(cal_score[i]),
                             threshold=b["threshold"], predicted=int(s >= b["threshold"]), arm=ARM,
                             genus_in_training_set=unl.genus[i] in train_genera,
                             family_in_training_set=unl.family[i] in train_families))
    out = pd.DataFrame(rows)
    mt = pd.read_csv(RES / "trees_xai/da/metrics_per_trait.csv")
    mt = mt[mt.arm == ARM].set_index(["split_level", "trait"])
    level = np.where(out.family_in_training_set, "genus", "family")
    key = list(zip(level, out.trait))
    out["cv_level_used"] = level
    out["cv_auprc"] = [mt.auprc.get(k, np.nan) for k in key]
    out["cv_mcc"] = [mt.mcc.get(k, np.nan) for k in key]
    out["cv_beats_taxonomic_null_q05"] = [bool(mt.q_vs_null.get(k, 1) < 0.05 and mt.delta_vs_null.get(k, 0) > 0)
                                          for k in key]
    out["call_reliable"] = [bool(mt.mcc_ci_low.get(k, -1) > 0.2 and mt.q_vs_null.get(k, 1) < 0.05
                                 and mt.delta_vs_null.get(k, 0) > 0) for k in key]
    lift = external_lift()
    out["external_auprc"] = out.trait.map({t: v[0] for t, v in lift.items()})
    out["external_auprc_lift"] = out.trait.map({t: v[1] for t, v in lift.items()})
    contradicted = out.external_auprc_lift.notna() & (out.external_auprc_lift <= 1)
    out.loc[contradicted, "call_reliable"] = False
    # fourth condition (22 October 2026): external MCC of the calls above 0.2 where external labels exist
    sc = pd.read_csv(RES / "additional/metabolic_bench/rf/scorecard.csv")
    ext_mcc = sc[(sc.set == "external") & (sc.metric == "mcc")].set_index("category").RedoxKO
    out["external_mcc"] = out.trait.map(ext_mcc)
    out.loc[out.external_mcc.notna() & (out.external_mcc <= 0.2), "call_reliable"] = False
    meta = json.loads((RES / "models" / ARM / "metadata.json").read_text())["reliable_traits_by_level"]
    for lv in ("genus", "family"):
        got = sorted(out.loc[(out.cv_level_used == lv) & out.call_reliable, "trait"].unique())
        if (out.cv_level_used == lv).any():
            assert got == sorted(meta[lv]), (lv, got, meta[lv])
    out.loc[~out.call_reliable.astype(bool), "predicted"] = np.nan
    out.to_csv(RES / f"predictions_unlabelled_{ARM}.csv", index=False, float_format="%.4f")
    wide = out.pivot(index=["genome_id", "organism", "genus_in_training_set"], columns="trait", values="predicted")
    print("traits with reliable calls:", sorted(out.loc[out.call_reliable.astype(bool), "trait"].unique()))
    wide.to_csv(RES / f"predictions_unlabelled_{ARM}_calls.csv")
    print(len(unl), "genomes;", "genus seen in training:", int(unl.genus.isin(train_genera).sum()),
          "; family seen:", int(unl.family.isin(train_families).sum()))
    print(wide.sum().sort_values(ascending=False).to_string())


if __name__ == "__main__":
    main()
