"""Balanced 100-genome external subset (src/78_balanced_subsets.py external; Table N5) for the final model RedoxKO
(random forest on the 852-KO panel). The subset itself is the one already on disk
(validation_new_organisms/balanced_external_100.txt, chosen from labels only); nothing is re-selected.

Macro AUPRC over the 9 evaluable external categories, all 562 genomes and the balanced 100, for the forest, the
taxonomic baseline (null_taxNN) and RedoxKO-LX (B2) as in validation_new_organisms/external_balanced_comparison.csv.
Forest scores: results/external_validation_missing_rows_predictions.csv (arm "Random forest": fitted on all 358
labelled genomes with min_samples_leaf by inner CV, src/138_external_missing_rows.py). The script checks that these
scores agree with the final bundle results/models/RF/model.joblib applied to the external KO panel.
Added (B2 table had no intervals): 95% genus-bootstrap intervals (1000 resamples, seed 20260917 as in
src/22_external_eval.py) for each arm's macro AUPRC and for the margin over the taxonomic baseline.
usage: python src/156_balanced_external_rf.py
Output: validation_new_organisms/external_balanced_comparison_rf.csv
"""
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import safe_ap  # noqa: E402

V, RES, EXT = ROOT / "validation_new_organisms", ROOT / "results", ROOT / "data/external"
SEED, N_BOOT = 20260917, 1000


def main():
    p = pd.read_csv(RES / "external_validation_predictions.csv")
    ev = pd.read_csv(RES / "external_validation_per_trait.csv")
    traits = sorted(set(ev[(ev.arm == "B2") & (ev.subset == "all") & ev.evaluable].trait))
    p = p[p.trait.isin(traits) & p.y.notna()]
    y = p[p.arm == "B2"].pivot(index="genome_id", columns="trait", values="y")
    S = {a: p[p.arm == a].pivot(index="genome_id", columns="trait", values="score") for a in ("B2", "null_taxNN")}
    rf = pd.read_csv(RES / "external_validation_missing_rows_predictions.csv")
    rf = rf[rf.arm == "Random forest"].pivot(index="genome_id", columns="trait", values="score")
    S["RF"] = rf.reindex(index=y.index, columns=traits)
    assert not S["RF"].isna().any().any()

    # check: the stored forest scores equal the final bundle applied to the external KO panel
    bundle = joblib.load(RES / "models/RF/model.joblib")
    K = pd.read_csv(EXT / "ko_matrix_panel.tsv", sep="\t", index_col=0).reindex(y.index)
    K = (K.fillna(0) > 0).astype(np.float32)
    chk = {t: np.corrcoef(bundle[t]["model"].predict_proba(K.reindex(columns=bundle[t]["model"].columns).fillna(0)),
                          S["RF"][t].to_numpy())[0, 1] for t in traits if t in bundle}
    mad = {t: float(np.abs(bundle[t]["model"].predict_proba(K.reindex(columns=bundle[t]["model"].columns).fillna(0))
                           - S["RF"][t].to_numpy()).max()) for t in traits if t in bundle}
    print("stored forest scores vs final bundle: Pearson r", {t: round(v, 4) for t, v in chk.items()})
    print("max |difference|", {t: round(v, 4) for t, v in mad.items()})

    lin = pd.read_csv(EXT / "lineage.tsv", sep="\t").set_index("genome_id")
    bal = [g for g in (V / "balanced_external_100.txt").read_text().split() if g in y.index]
    print(f"balanced subset: {len(bal)} genomes")
    rows = []
    for subset, ids in (("all 562", list(y.index)), ("balanced 100", bal)):
        Ys = y.loc[ids]
        genus = lin.ncbi_genus.reindex(ids).fillna("NA").to_numpy()
        ug = np.unique(genus)
        mem = {g: np.where(genus == g)[0] for g in ug}
        rng = np.random.default_rng(SEED)
        boots = [np.concatenate([mem[g] for g in rng.choice(ug, len(ug))]) for _ in range(N_BOOT)]

        def macro(arm, idx):
            v = []
            for t in traits:
                yy, ss = Ys[t].to_numpy()[idx], S[arm].loc[ids, t].to_numpy()[idx]
                m = ~np.isnan(yy)
                v.append(safe_ap(yy[m], ss[m]))
            return float(np.nanmean(v))
        allidx = np.arange(len(ids))
        B = {a: np.array([macro(a, b) for b in boots]) for a in S}
        for a in S:
            r = dict(subset=subset, arm=a, n_genomes=len(ids), macro_auprc=macro(a, allidx),
                     ci_low=np.percentile(B[a], 2.5), ci_high=np.percentile(B[a], 97.5))
            if a != "null_taxNN":
                d = B[a] - B["null_taxNN"]
                r.update(margin_vs_taxonomic=r["macro_auprc"] - macro("null_taxNN", allidx),
                         margin_ci_low=np.percentile(d, 2.5), margin_ci_high=np.percentile(d, 97.5))
            if a == "RF":
                d = B["RF"] - B["B2"]
                r.update(rf_minus_b2=r["macro_auprc"] - macro("B2", allidx), rf_minus_b2_ci_low=np.percentile(d, 2.5),
                         rf_minus_b2_ci_high=np.percentile(d, 97.5))
            rows.append(r)
    out = pd.DataFrame(rows)
    old = pd.read_csv(V / "external_balanced_comparison.csv").set_index("arm")
    for a, col in (("B2", "all 562"), ("B2", "balanced 100"), ("null_taxNN", "all 562"), ("null_taxNN", "balanced 100")):
        new = out[(out.arm == a) & (out.subset == col)].macro_auprc.iloc[0]
        assert abs(new - old.loc[a, col]) < 1e-6, (a, col, new, old.loc[a, col])
    print("B2 and taxonomic-baseline values reproduce validation_new_organisms/external_balanced_comparison.csv")
    out.to_csv(V / "external_balanced_comparison_rf.csv", index=False, float_format="%.4f")
    print(out.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
