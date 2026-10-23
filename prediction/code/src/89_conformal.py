"""X3 (docs/preregistration_additional_experiments.md): per-prediction uncertainty for B2 by Mondrian split-conformal
prediction. Calibration: family-held-out out-of-fold B2 scores of the 358 training genomes. Non-conformity
1 - p(class), class-conditional quantiles per category, target coverage 90 %. Evaluated on the external 562 set
(Madin labels) and on the primary 100 newly described species (blind growth tests). Nothing is refitted.

usage: python src/89_conformal.py
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/additional"
EPS = 0.10


def quantiles(oof):
    q = {}
    for t, d in oof.groupby("trait"):
        q[t] = {}
        for c in (0, 1):
            a = np.sort((1 - d.prob[d.y == 1]) if c == 1 else d.prob[d.y == 0])  # 1 - p(class c)
            n = len(a)
            k = math.ceil((n + 1) * (1 - EPS))
            q[t][c] = float(a[min(k, n) - 1]) if n else 1.0
    return q


def sets(df, q):
    inc1 = np.array([(1 - p) <= q[t][1] for t, p in zip(df.trait, df.p)])
    inc0 = np.array([p <= q[t][0] for t, p in zip(df.trait, df.p)])
    df = df.copy()
    df["in1"], df["in0"] = inc1, inc0
    df["set"] = np.select([inc1 & ~inc0, inc0 & ~inc1, inc0 & inc1], ["{1}", "{0}", "{0,1}"], "{}")
    df["covered"] = np.where(df.y == 1, df.in1, df.in0)
    return df


def summarize(df, name):
    single = df[df.set.isin(["{0}", "{1}"])]
    pred = (single.set == "{1}").astype(int)
    row = dict(evaluation=name, calls=len(df), coverage=df.covered.mean(),
               coverage_pos=df[df.y == 1].covered.mean(), coverage_neg=df[df.y == 0].covered.mean(),
               frac_singleton=len(single) / len(df), frac_uncertain=(df.set == "{0,1}").mean(),
               frac_empty=(df.set == "{}").mean(),
               singleton_accuracy=float((pred.to_numpy() == single.y.to_numpy()).mean()) if len(single) else np.nan)
    return row


def main():
    oof = pd.read_parquet(ROOT / "results/oof/B2.parquet")
    oof = oof[(oof.arm == "B2") & (oof.level == "family")][["trait", "genome_id", "y", "prob"]]
    q = quantiles(oof)

    ext = pd.read_csv(ROOT / "results/external_validation_predictions.csv")
    ext = ext[ext.arm == "B2"].rename(columns={"score": "p"})[["genome_id", "trait", "y", "p"]].dropna()
    new = pd.read_csv(ROOT / "validation_new_organisms/scored178_B2_balanced100.csv")
    new = new[new.blind & ~new.category.isin(["donor_other", "acceptor_other"])]
    new = new.rename(columns={"category": "trait", "expected": "y", "score": "p", "accession": "genome_id"})
    new = new[["genome_id", "trait", "y", "p"]]

    res, per = [], []
    for name, df in [("training OOF (family, in-sample check)", oof.rename(columns={"prob": "p"})),
                     ("external 562", ext), ("new species 100", new)]:
        s = sets(df, q)
        res.append(summarize(s, name))
        for t, d in s.groupby("trait"):
            per.append(dict(summarize(d, name), trait=t))
        if name == "new species 100":
            s.to_csv(OUT / "x3_new_species_sets.csv", index=False)
    res, per = pd.DataFrame(res), pd.DataFrame(per)
    res.to_csv(OUT / "x3_conformal_summary.csv", index=False)
    per.to_csv(OUT / "x3_conformal_per_trait.csv", index=False)
    (OUT / "x3_quantiles.json").write_text(json.dumps(q, indent=1))
    pd.set_option("display.width", 200)
    print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
