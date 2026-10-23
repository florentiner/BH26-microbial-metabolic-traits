"""Metrics of RedoxKO and RedoxKO-ens-BF recomputed on autotroph-like test genomes only, next to the full test sets.

Autotroph-like: at least one carbon-fixation KEGG module >= 80% complete (src/133_autotroph_subset.py).
External genomes: macro AUPRC, macro AUROC and macro MCC over the 9 categories, calls at the final models' thresholds.
New species: macro AUROC over categories, pooled MCC, sensitivity and specificity over the 255 growth experiments, calls
as scored for the paper. Intervals: 95% bootstrap over genera (external) or papers (new species), 2,000 resamples.
usage: python src/134_autotroph_only_metrics.py
Output: results/additional/autotroph_only_metrics.csv
"""
import importlib
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import PROC, load_labelled, safe_ap  # noqa: E402

a133 = importlib.import_module("133_autotroph_subset")
SEED, B = 20261021, 2000
ARMS = {"RedoxKO": "B2", "RedoxKO-ens-BF": "B2XLE_large"}


def ext_metrics(d):
    ap, au, mc = [], [], []
    for _, g in d.groupby("trait"):
        if 0 < g.y.sum() < len(g):
            ap.append(safe_ap(g.y.to_numpy(), g.score.to_numpy()))
            au.append(roc_auc_score(g.y, g.score))
            mc.append(matthews_corrcoef(g.y, g.call))
    return dict(macro_auprc=np.mean(ap), macro_auroc=np.mean(au), macro_mcc=np.mean(mc))


def new_metrics(d):
    au = [roc_auc_score(g.y, g.score) for _, g in d.groupby("trait") if 0 < g.y.sum() < len(g)]
    return dict(macro_auroc=np.mean(au), mcc=matthews_corrcoef(d.y, d.call),
                sensitivity=(d.call[d.y == 1] == 1).mean(), specificity=(d.call[d.y == 0] == 0).mean())


def boot(d, cl, fn):
    rng = np.random.default_rng(SEED)
    u = d[cl].unique()
    idx = {c: np.where(d[cl].to_numpy() == c)[0] for c in u}
    res = [fn(d.iloc[np.concatenate([idx[c] for c in rng.choice(u, len(u))])]) for _ in range(B)]
    return {k: np.nanpercentile([r[k] for r in res], [2.5, 97.5]) for k in res[0]}


def main():
    lab, _, _ = load_labelled()
    train = pd.read_csv(PROC / "ko_matrix_panel.tsv", sep="\t", index_col=0).loc[lab.genome_id]
    cols = train.columns
    a_ex, _ = a133.autotroph_like(pd.read_csv(ROOT / "data/external/ko_matrix_panel.tsv", sep="\t", index_col=0)
                                  .reindex(columns=cols).fillna(0))
    a_nw, _ = a133.autotroph_like(pd.read_parquet(ROOT / "validation_new_organisms/ko_panel_validation.parquet")
                                  .reindex(columns=cols).fillna(0))
    lin = pd.read_csv(ROOT / "data/external/lineage.tsv", sep="\t")[["genome_id", "ncbi_genus"]]
    allx = pd.read_csv(ROOT / "results/external_validation_predictions.csv")
    rows = []
    for name, arm in ARMS.items():
        thr = {t: float(b["threshold"]) for t, b in joblib.load(ROOT / f"results/models/{arm}/model.joblib").items()}
        e = allx[allx.arm == arm].dropna(subset=["y"]).merge(lin, on="genome_id", how="left")
        e["call"] = [int(s >= thr[t]) for s, t in zip(e.score, e.trait)]
        e["auto"] = e.genome_id.map(a_ex).fillna(False).astype(bool)
        n = pd.read_csv(ROOT / f"validation_new_organisms/scored178_{arm}_balanced100.csv")
        n = n[~n.category.isin(["donor_other", "acceptor_other"])].rename(columns={"category": "trait", "expected": "y"})
        n["auto"] = n.accession.map(a_nw).fillna(False).astype(bool)
        for setn, d, gid, cl, fn in [("external genomes", e, "genome_id", "ncbi_genus", ext_metrics),
                                     ("new species", n, "accession", "paper", new_metrics)]:
            for sub, s in [("all", d), ("autotroph-like only", d[d.auto])]:
                m, ci = fn(s), boot(s, cl, fn)
                for k, v in m.items():
                    rows.append(dict(model=name, set=setn, subset=sub, genomes=s[gid].nunique(), tests=len(s),
                                     metric=k, value=v, ci_low=ci[k][0], ci_high=ci[k][1]))
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results/additional/autotroph_only_metrics.csv", index=False, float_format="%.4f")
    w = out.pivot_table(index=["model", "set", "metric"], columns="subset", values="value").round(3)
    w["change"] = (w["autotroph-like only"] - w["all"]).round(3)
    pd.set_option("display.width", 200)
    print(w.to_string())


if __name__ == "__main__":
    main()
