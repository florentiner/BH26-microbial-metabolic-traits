"""Head-to-head of every saved task-1 model on the same balanced, diverse 100 novel organisms.

The question this answers: which model is best on genomes outside the gold dataset? Every candidate is the saved
final model (trained on all 358 labelled genomes) applied unchanged to the same 100 organisms (one per genus, 18
phyla; validation_new_organisms/balanced_prospective_primary_100.txt), scored against the same blind, quoted growth
tests. Two views:

* threshold-free ranking skill — per-category AUROC and AUPRC over the scored calls, averaged across categories
  that have both a tested positive and a tested negative. This compares models fairly regardless of how each one
  picked its cut-off;
* the yes/no calls each model would actually release — accuracy, MCC, sensitivity, specificity.

Paired differences use a cluster bootstrap over description papers.

usage: python src/82_novel_model_comparison.py
Output: validation_new_organisms/novel100_model_comparison.csv (+ per-call scores)
"""
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
V = ROOT / "validation_new_organisms"
K = ROOT / "data/kaggle/validation178"
ARMS = {"B2": "KO only (released)", "B2XL": "KO only, two learners",
        "B2XLE_large": "KO + Bacformer-large (DL fusion)", "E_large": "Bacformer-large alone (DL)",
        "C_large": "ESM-C 300M alone (DL)"}
CATCH_ALL = {"donor_other", "acceptor_other"}


def features(accs):
    ko = pd.read_parquet(V / "ko_panel_validation.parquet")
    bfl = pd.read_parquet(K / "genome_features_bacformer_large_mean.parquet")
    esm = pd.read_parquet(K / "genome_features_esmc_mean.parquet")
    esx = pd.read_parquet(K / "genome_features_esmc_max.parquet")
    return ko.join(bfl, how="left").join(esm, how="left").join(esx, how="left").loc[accs]


def calls(recs, accs):
    rows = []
    for r in recs:
        if r["refseq_accession"] not in accs or r.get("_source_file") == "candidates.json":
            continue
        for c, ev in r["categories"].items():
            if c in CATCH_ALL or ev.get("status") not in ("tested_positive", "tested_negative"):
                continue
            rows.append(dict(accession=r["refseq_accession"], paper=str(r.get("pmid") or r["species"]),
                             category=c, expected=int(ev["status"] == "tested_positive")))
    return pd.DataFrame(rows)


def mcc(y, p):
    y, p = np.asarray(y), np.asarray(p)
    tp, tn = int(((y == 1) & (p == 1)).sum()), int(((y == 0) & (p == 0)).sum())
    fp, fn = int(((y == 0) & (p == 1)).sum()), int(((y == 1) & (p == 0)).sum())
    d = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    return (tp * tn - fp * fn) / d if d else 0.0


def ranking(d, col):
    """Mean per-category AUROC / AUPRC over categories with both classes."""
    au, ap = [], []
    for _, g in d.groupby("category"):
        if 0 < g.expected.sum() < len(g):
            au.append(roc_auc_score(g.expected, g[col]))
            ap.append(average_precision_score(g.expected, g[col]))
    return float(np.mean(au)), float(np.mean(ap)), len(au)


def boot(d, stat, n=2000, seed=20260921):
    rng = np.random.default_rng(seed)
    papers = d.paper.unique()
    idx = {p: np.where(d.paper.to_numpy() == p)[0] for p in papers}
    out = []
    for _ in range(n):
        take = np.concatenate([idx[p] for p in rng.choice(papers, len(papers), replace=True)])
        try:
            out.append(stat(d.iloc[take]))
        except ValueError:
            continue
    return np.percentile(out, 2.5), np.percentile(out, 97.5)


def main():
    accs = open(V / "balanced_prospective_primary_100.txt").read().split()
    recs = json.loads((V / "validation_set.json").read_text())
    X = features(accs)
    d = calls(recs, set(accs))
    print(f"{d.accession.nunique()} organisms, {len(d)} calls, {d.paper.nunique()} papers, "
          f"{d.category.nunique()} categories")
    for arm in ARMS:
        bundle = joblib.load(ROOT / f"results/models/{arm}/model.joblib")
        sc, th = {}, {}
        for t, b in bundle.items():
            sc[t] = pd.Series(np.asarray(b["model"].predict_proba(X)), index=X.index)
            th[t] = float(b["threshold"])
        d[f"score_{arm}"] = [sc[c][a] if c in sc else np.nan for a, c in zip(d.accession, d.category)]
        d[f"call_{arm}"] = [int(sc[c][a] >= th[c]) if c in sc else 0 for a, c in zip(d.accession, d.category)]
    d.to_csv(V / "novel100_model_scores.csv", index=False)

    rows = []
    for arm, desc in ARMS.items():
        au, ap, nc = ranking(d, f"score_{arm}")
        agree = (d[f"call_{arm}"] == d.expected).astype(int)
        pos, neg = d.expected == 1, d.expected == 0
        rows.append(dict(model=arm, what=desc, auroc=au, auprc=ap, categories=nc, accuracy=agree.mean(),
                         mcc=mcc(d.expected, d[f"call_{arm}"]),
                         sensitivity=agree[pos].mean(), specificity=agree[neg].mean()))
    t = pd.DataFrame(rows).sort_values("auroc", ascending=False)
    ref = "B2XLE_large"
    for arm in ARMS:
        if arm == ref:
            continue
        lo, hi = boot(d, lambda x: ranking(x, f"score_{arm}")[0] - ranking(x, f"score_{ref}")[0])
        t.loc[t.model == arm, f"auroc_minus_{ref}"] = f"{ranking(d, f'score_{arm}')[0] - ranking(d, f'score_{ref}')[0]:+.3f} [{lo:+.3f}, {hi:+.3f}]"
    t.to_csv(V / "novel100_model_comparison.csv", index=False, float_format="%.4f")
    pd.set_option("display.width", 220)
    print(t.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
