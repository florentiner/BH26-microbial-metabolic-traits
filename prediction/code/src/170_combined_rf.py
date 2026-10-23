"""Combined models built on the final model RedoxKO (random forest on the 852-KO panel). Not registered.

Every combination of RedoxKO with one or more of the other models (RedoxKO-LX, RedoxKO-ens, logistic regression on
the Bacformer-large, ESM-C or ESM-2 genome vector) is scored by averaging the probabilities, as the combined models
of Table 1 do. Cross-validation uses the existing out-of-fold predictions (same folds for every model); the external
genomes use the existing external scores. Differences from RedoxKO are paired genus- or family-cluster bootstraps.
Because the best combination is chosen on the same out-of-fold data it is scored on, its cross-validation value is
optimistic; the external genomes are the check.
Output: results/combined_rf.csv
"""
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import RES, load_labelled, oof_matrix, paired_macro_delta, safe_ap  # noqa: E402

POOL = {"RedoxKO-LX": "B2", "RedoxKO-ens": "B2XL", "Bacformer-large": "E_large", "ESM-C": "C_large", "ESM-2": "C"}
OOF_DIR = {"RF": RES / "trees_xai/da/oof"}


def ext_scores():
    e = pd.read_csv(RES / "external_validation_predictions.csv")
    m = pd.read_csv(RES / "external_validation_missing_rows_predictions.csv")
    m["arm"] = m.arm.map({"Random forest": "RF", "ESM-2": "C"}).fillna(m.arm)
    y = e[e.arm == "B2"][["genome_id", "trait", "y"]]
    s = pd.concat([e[["genome_id", "trait", "arm", "score"]], m[["genome_id", "trait", "arm", "score"]]])
    w = s.pivot_table(index=["genome_id", "trait"], columns="arm", values="score").reset_index()
    return w.merge(y, on=["genome_id", "trait"])


def ext_macro(w, arms):
    p = w[list(arms)].mean(axis=1)
    return float(np.mean([safe_ap(g.y.to_numpy(), p[g.index].to_numpy()) for _, g in w.groupby("trait")]))


def main():
    lab, traits, Y = load_labelled()
    obs = ~np.isnan(Y)
    W = ext_scores()
    P = {lv: {a: oof_matrix(OOF_DIR.get(a, RES / "oof"), a, lv, lab, traits) for a in ["RF", *POOL.values()]}
         for lv in ("genus", "family")}
    rows = []
    names = list(POOL)
    for k in range(0, len(names) + 1):
        for combo in itertools.combinations(names, k):
            arms = ["RF", *[POOL[n] for n in combo]]
            row = {"model": " + ".join(["RedoxKO", *combo])}
            for lv in ("genus", "family"):
                Pa = np.mean([P[lv][a] for a in arms], axis=0)
                row[f"{lv}_cv"] = float(np.nanmean([safe_ap(Y[obs[:, j], j], Pa[obs[:, j], j]) for j in range(len(traits))]))
                if k:
                    d, lo, hi = paired_macro_delta(lab, traits, Y, Pa, P[lv]["RF"], lv)
                    row.update({f"{lv}_minus_RedoxKO": d, f"{lv}_ci_low": lo, f"{lv}_ci_high": hi})
            row["external"] = ext_macro(W, arms)
            rows.append(row)
            print(row, flush=True)
    d = pd.DataFrame(rows).sort_values("genus_cv", ascending=False)
    d.to_csv(RES / "combined_rf.csv", index=False, float_format="%.4f")
    print(d.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
