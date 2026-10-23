"""Pre-registration AT (docs/preregistration_AT_autotroph_test.md): RedoxKO on 100 autotroph-like new species.

Species: the autotroph-like species among the 167 scored new species (blind evidence, catch-all categories excluded),
plus newly harvested autotroph-like species (validation_new_organisms/pool_auto, audited evidence) in harvest order
until 100 species are reached. Primary endpoint: macro AUROC of RedoxKO with a 95% bootstrap interval over papers.
Secondary: pooled MCC, sensitivity, specificity.
usage: python src/142_autotroph_test.py
Output: results/additional/rf_paper/autotroph_test_100_rf.csv (+ _tests.csv)
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import PROC  # noqa: E402

a133 = importlib.import_module("133_autotroph_subset")
V = ROOT / "validation_new_organisms"


def stats(x):
    au = [roc_auc_score(g.expected, g.score) for _, g in x.groupby("category") if 0 < g.expected.sum() < len(g)]
    return dict(macro_auroc=np.mean(au), categories=len(au), mcc=matthews_corrcoef(x.expected, x.call),
                sensitivity=(x.call[x.expected == 1] == 1).mean(), specificity=(x.call[x.expected == 0] == 0).mean())


def main():
    cols = pd.read_csv(PROC / "ko_matrix_panel.tsv", sep="\t", index_col=0).columns
    K = pd.concat([pd.read_parquet(V / "ko_panel_validation.parquet"), pd.read_parquet(V / "ko_panel_auto.parquet")])
    K = K[~K.index.duplicated()]
    auto, _ = a133.autotroph_like(K.reindex(columns=cols).fillna(0))
    old = pd.read_csv(ROOT / "results/additional/rf_paper/scored178_rf.csv").drop(columns=["score", "call"]).rename(
        columns={"score_RF": "score", "call_RF": "call"})
    old = old[old.blind & ~old.category.isin(["donor_other", "acceptor_other"])]
    old = old[old.accession.map(auto).fillna(False).astype(bool)][["accession", "paper", "category", "expected", "score", "call"]]
    old["source"] = "existing"
    new_t = pd.read_csv(V / "pool_auto/tests_auto.csv")
    import joblib
    Ka = (pd.read_parquet(V / "ko_panel_auto.parquet") > 0).astype(np.float32)
    sc = []
    for t_, b in joblib.load(ROOT / "results/models/RF/model.joblib").items():
        p_ = np.asarray(b["model"].predict_proba(Ka.reindex(columns=b["model"].columns).fillna(0)))
        sc += [dict(accession=a, category=t_, score=v, call=int(v >= b["threshold"])) for a, v in zip(Ka.index, p_)]
    sc = pd.DataFrame(sc)
    new = new_t.merge(sc[["accession", "category", "score", "call"]], on=["accession", "category"])
    new["source"] = "harvested"
    order = [r["accession"] for f in ("pool.json", "pool2.json") for r in json.load(open(V / "pool_auto" / f))]
    need = 100 - old.accession.nunique()
    take = [a for a in order if a in set(new.accession)][:max(need, 0)]
    d = pd.concat([old, new[new.accession.isin(take)][old.columns]], ignore_index=True)
    print(f"species {d.accession.nunique()} (existing {old.accession.nunique()}, harvested {len(take)}), tests {len(d)}")
    rng = np.random.default_rng(20261021)
    papers = d.paper.unique()
    idx = {p: np.where(d.paper.to_numpy() == p)[0] for p in papers}
    boots = []
    for _ in range(2000):
        try:
            boots.append(stats(d.iloc[np.concatenate([idx[p] for p in rng.choice(papers, len(papers))])]))
        except ValueError:
            pass
    rows = []
    for sub, x in [("100 autotroph-like species", d), ("existing only", d[d.source == "existing"]),
                   ("harvested only", d[d.source == "harvested"])]:
        s = stats(x)
        r = dict(subset=sub, species=x.accession.nunique(), tests=len(x), **s)
        if sub.startswith("100"):
            for k in ("macro_auroc", "mcc", "sensitivity", "specificity"):
                r[k + "_lo"], r[k + "_hi"] = np.nanpercentile([b[k] for b in boots], [2.5, 97.5])
        rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results/additional/rf_paper/autotroph_test_100_rf.csv", index=False, float_format="%.4f")
    d.to_csv(ROOT / "results/additional/rf_paper/autotroph_test_100_tests_rf.csv", index=False)
    pd.set_option("display.width", 250)
    print(out.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
