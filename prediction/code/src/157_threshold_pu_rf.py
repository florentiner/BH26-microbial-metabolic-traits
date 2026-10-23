"""Threshold correction by the estimated reporting rate (src/99_threshold_pu.py, docs/preregistration_threshold_correction.md)
for the final model RedoxKO (random forest on the 852-KO panel). Same method, sets, seed and bootstrap as the original;
only the model changes:
  * label frequency c and the Elkan & Noto score map g: fitted on the forest's family-held-out out-of-fold scores
    (results/trees_xai/da/oof/RF.parquet, level family) instead of results/oof/B2.parquet;
  * original threshold: the forest's released threshold (results/models/RF/model.joblib);
  * new-species scores and calls: score_RF / call_RF of results/additional/rf_paper/scored178_rf.csv (the 43 withheld
    species = blind species not in the primary 100, identical to results/additional/rf_paper/holdout_calls_rf.csv);
  * external 562: forest scores of results/external_validation_missing_rows_predictions.csv (arm "Random forest"),
    labels from results/external_validation_predictions.csv (same genome/category rows as for B2).
No model is refitted.
usage: python src/157_threshold_pu_rf.py
Output: results/additional/threshold_correction_rf/threshold_correction_results.csv, label_frequency_c.csv
"""
import importlib
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
t99 = importlib.import_module("99_threshold_pu")
OUT = ROOT / "results/additional/threshold_correction_rf"
VAL = ROOT / "validation_new_organisms"
ARM = "RF"


def calibrators():
    oof = pd.read_parquet(ROOT / "results/trees_xai/da/oof/RF.parquet")
    oof = oof[(oof.arm == ARM) & (oof.level == "family")]
    bundle = joblib.load(ROOT / "results/models/RF/model.joblib")
    out = {}
    for t, d in oof.groupby("trait"):
        s, y = d.prob.to_numpy().reshape(-1, 1), d.y.to_numpy().astype(int)
        g = LogisticRegression(C=1e6).fit(s, y)
        c = float(np.clip(g.predict_proba(s[y == 1])[:, 1].mean(), 0.3, 1.0))
        thr = float(bundle[t]["threshold"])
        tau = float(g.predict_proba([[thr]])[0, 1])
        out[t] = {"g": g, "c": c, "threshold_raw": thr, "tau": tau}
    return out


def row(name, d, cluster):
    a, b = t99.metrics(d.y, d.c_orig), t99.metrics(d.y, d.c_corr)
    ci = t99.boot(d, cluster)
    return dict(model=ARM, test_set=name, calls=len(d), organisms=(d.accession if "accession" in d else d.genome_id).nunique(),
                **{f"{k}_original": v for k, v in a.items() if k != "n"},
                **{f"{k}_corrected": v for k, v in b.items() if k != "n"},
                mcc_change_ci=ci["mcc"], sensitivity_change_ci=ci["sensitivity"], specificity_change_ci=ci["specificity"])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    blind = pd.read_csv(VAL / "scored178_B2.csv")
    blind = sorted(blind[blind.blind == True].accession.unique())  # noqa: E712
    primary = set((VAL / "balanced_prospective_primary_100.txt").read_text().split())
    holdout = [a for a in blind if a not in primary]
    cal = calibrators()
    cvals = [dict(model=ARM, category=t, c=v["c"], threshold_raw=v["threshold_raw"], tau=v["tau"],
                  corrected_calibrated_cutoff=v["c"] * v["tau"]) for t, v in cal.items()]
    sc = pd.read_csv(ROOT / "results/additional/rf_paper/scored178_rf.csv")
    sc = sc[(sc.blind == True) & ~sc.category.isin(t99.CATCH_ALL)]  # noqa: E712
    assert (sc.call_RF == (sc.score_RF >= sc.category.map({t: v["threshold_raw"] for t, v in cal.items()})).astype(int)).all()
    hc = pd.read_csv(ROOT / "results/additional/rf_paper/holdout_calls_rf.csv")
    h = sc[sc.accession.isin(holdout)]
    assert len(h) == len(hc) and set(h.accession) == set(hc.accession), "43 withheld species differ from holdout_calls_rf.csv"
    rows = []
    for name, d in {"held-out 43 new species (primary)": h,
                    "primary 100 new species (not independent)": sc[sc.accession.isin(primary)]}.items():
        d = d.assign(y=d.expected, c_orig=d.call_RF,
                     c_corr=[t99.corrected(cal, t, [s])[0] for t, s in zip(d.category, d.score_RF)])
        rows.append(row(name, d, "paper"))
    ext = pd.read_csv(ROOT / "results/external_validation_predictions.csv")
    e = ext[ext.arm == "B2"].dropna(subset=["y"])[["genome_id", "trait", "y"]]
    rf = pd.read_csv(ROOT / "results/external_validation_missing_rows_predictions.csv")
    rf = rf[rf.arm == "Random forest"][["genome_id", "trait", "score"]]
    e = e.merge(rf, on=["genome_id", "trait"], how="left", validate="1:1")
    assert e.score.notna().all()
    e = e[e.trait.isin(cal)].copy()
    ext_lab = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t").set_index("genome_id")
    e["c_orig"] = [int(s >= cal[t]["threshold_raw"]) for t, s in zip(e.trait, e.score)]
    e["c_corr"] = [t99.corrected(cal, t, [s])[0] for t, s in zip(e.trait, e.score)]
    e["genus"] = e.genome_id.map(ext_lab.genus)
    rows.append(row("external 562 (Madin labels)", e, "genus"))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "threshold_correction_results.csv", index=False)
    pd.DataFrame(cvals).to_csv(OUT / "label_frequency_c.csv", index=False)
    pd.set_option("display.width", 250)
    show = ["model", "test_set", "calls", "accuracy_original", "accuracy_corrected", "mcc_original", "mcc_corrected",
            "mcc_change_ci", "sensitivity_original", "sensitivity_corrected", "specificity_original",
            "specificity_corrected"]
    print(R[show].round(3).to_string(index=False))
    print(pd.DataFrame(cvals)[["category", "c", "tau", "corrected_calibrated_cutoff"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
