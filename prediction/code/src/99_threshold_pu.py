"""C (docs/preregistration_threshold_correction.md): positive-unlabeled (Elkan & Noto e1) threshold correction for
the released donor/acceptor models, estimated on the family-held-out out-of-fold scores of the 358 training genomes
only, then applied unchanged to the 43 held-out new species (primary), the external 562 genomes (secondary) and the
primary 100 new species (reported only). No model is refitted.

usage: python src/99_threshold_pu.py
"""
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import matthews_corrcoef

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/additional/threshold_correction"
VAL = ROOT / "validation_new_organisms"
SEED = 20261014
MODELS = ["B2", "B2XLE_large"]
CATCH_ALL = ["donor_other", "acceptor_other"]


def calibrators(arm):
    oof = pd.read_parquet(ROOT / f"results/oof/{arm}.parquet")
    oof = oof[(oof.arm == arm) & (oof.level == "family")]
    bundle = joblib.load(ROOT / f"results/models/{arm}/model.joblib")
    out = {}
    for t, d in oof.groupby("trait"):
        s, y = d.prob.to_numpy().reshape(-1, 1), d.y.to_numpy().astype(int)
        g = LogisticRegression(C=1e6).fit(s, y)
        c = float(np.clip(g.predict_proba(s[y == 1])[:, 1].mean(), 0.3, 1.0))
        thr = float(bundle[t]["threshold"])
        tau = float(g.predict_proba([[thr]])[0, 1])
        out[t] = {"g": g, "c": c, "threshold_raw": thr, "tau": tau}
    return out


def corrected(cal, trait, score):
    k = cal[trait]
    return (k["g"].predict_proba(np.asarray(score, float).reshape(-1, 1))[:, 1] / k["c"] >= k["tau"]).astype(int)


def metrics(y, c):
    y, c = np.asarray(y, int), np.asarray(c, int)
    tp, tn = ((y == 1) & (c == 1)).sum(), ((y == 0) & (c == 0)).sum()
    return dict(n=len(y), accuracy=(tp + tn) / len(y),
                mcc=matthews_corrcoef(y, c) if len(set(c)) > 1 and len(set(y)) > 1 else 0.0,
                sensitivity=tp / max((y == 1).sum(), 1), specificity=tn / max((y == 0).sum(), 1))


def boot(df, cluster, n=1000):
    rng = np.random.default_rng(SEED)
    ids = df[cluster].unique()
    idx = {i: np.where(df[cluster].to_numpy() == i)[0] for i in ids}
    out = {"mcc": [], "sensitivity": [], "specificity": []}
    for _ in range(n):
        d = df.iloc[np.concatenate([idx[i] for i in rng.choice(ids, len(ids))])]
        a, b = metrics(d.y, d.c_orig), metrics(d.y, d.c_corr)
        for k in out:
            out[k].append(b[k] - a[k])
    return {k: np.percentile(v, [2.5, 97.5]).round(3).tolist() for k, v in out.items()}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    blind = pd.read_csv(VAL / "scored178_B2.csv")
    blind = sorted(blind[blind.blind == True].accession.unique())  # noqa: E712
    primary = set((VAL / "balanced_prospective_primary_100.txt").read_text().split())
    holdout = [a for a in blind if a not in primary]
    ext = pd.read_csv(ROOT / "results/external_validation_predictions.csv")
    ext_lab = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t").set_index("genome_id")
    rows, cvals = [], []
    for arm in MODELS:
        cal = calibrators(arm)
        cvals += [dict(model=arm, category=t, c=v["c"], threshold_raw=v["threshold_raw"], tau=v["tau"],
                       corrected_calibrated_cutoff=v["c"] * v["tau"]) for t, v in cal.items()]
        sc = pd.read_csv(VAL / f"scored178_{arm}.csv")
        sc = sc[(sc.blind == True) & ~sc.category.isin(CATCH_ALL)]  # noqa: E712
        sets = {"held-out 43 new species (primary)": sc[sc.accession.isin(holdout)],
                "primary 100 new species (not independent)": sc[sc.accession.isin(primary)]}
        for name, d in sets.items():
            d = d.assign(y=d.expected, c_orig=d.call, c_corr=[corrected(cal, t, [s])[0]
                                                               for t, s in zip(d.category, d.score)])
            a, b = metrics(d.y, d.c_orig), metrics(d.y, d.c_corr)
            ci = boot(d, "paper")
            rows.append(dict(model=arm, test_set=name, calls=len(d), organisms=d.accession.nunique(),
                             **{f"{k}_original": v for k, v in a.items() if k != "n"},
                             **{f"{k}_corrected": v for k, v in b.items() if k != "n"},
                             mcc_change_ci=ci["mcc"], sensitivity_change_ci=ci["sensitivity"],
                             specificity_change_ci=ci["specificity"]))
        e = ext[ext.arm == arm].dropna(subset=["y"]).copy()
        e = e[e.trait.isin(cal)]
        e["c_orig"] = [int(s >= cal[t]["threshold_raw"]) for t, s in zip(e.trait, e.score)]
        e["c_corr"] = [corrected(cal, t, [s])[0] for t, s in zip(e.trait, e.score)]
        e["genus"] = e.genome_id.map(ext_lab.genus)
        a, b = metrics(e.y, e.c_orig), metrics(e.y, e.c_corr)
        ci = boot(e.rename(columns={"y": "y"}), "genus")
        rows.append(dict(model=arm, test_set="external 562 (Madin labels)", calls=len(e),
                         organisms=e.genome_id.nunique(),
                         **{f"{k}_original": v for k, v in a.items() if k != "n"},
                         **{f"{k}_corrected": v for k, v in b.items() if k != "n"},
                         mcc_change_ci=ci["mcc"], sensitivity_change_ci=ci["sensitivity"],
                         specificity_change_ci=ci["specificity"]))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "threshold_correction_results.csv", index=False)
    pd.DataFrame(cvals).to_csv(OUT / "label_frequency_c.csv", index=False)
    pd.set_option("display.width", 250)
    show = ["model", "test_set", "calls", "accuracy_original", "accuracy_corrected", "mcc_original", "mcc_corrected",
            "mcc_change_ci", "sensitivity_original", "sensitivity_corrected", "specificity_original",
            "specificity_corrected"]
    print(R[show].round(3).to_string(index=False))
    print(pd.DataFrame(cvals)[lambda d: d.model == "B2"][["category", "c", "tau", "corrected_calibrated_cutoff"]]
          .round(3).to_string(index=False))


if __name__ == "__main__":
    main()
