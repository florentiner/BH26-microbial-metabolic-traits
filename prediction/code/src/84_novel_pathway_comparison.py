"""Task 2 head-to-head: every saved pathway model on the same balanced 100 novel organisms.

Mirrors src/82_novel_model_comparison.py for the pathway task. Inputs are the full KO matrix saved by
src/79_pathway_validation_predict.py for the balanced 100 (validation_new_organisms/pathway100_ko_matrix.parquet),
the module-completeness features derived from it, and the Bacformer-large vectors already computed for these
genomes. Every saved final pathway model (results/models/pathways/<arm>/model.joblib, all trained on the 3 378
labelled pathway genomes) is applied unchanged and scored against the same blind growth tests, mapped to pathway
traits by the pre-registered exact mapping.

usage: python src/84_novel_pathway_comparison.py
Output: validation_new_organisms/novel100_pathway_comparison.csv (+ per-call scores)
"""
import importlib
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
EXACT = {"acceptor_sulfate": "sulfate_reduction", "acceptor_sulfur": "sulfur_reduction",
         "donor_sulfur_oxidation": "sulfur_oxidation", "acceptor_nitrate": "nitrate_reduction",
         "acceptor_iron": "iron_reduction", "donor_iron_oxidation": "iron_oxidation",
         "donor_hydrogen_oxidation": "hydrogen_oxidation"}
DESC = {"PB3": "full KO profile + module completeness (released)", "PB2": "full KO profile",
        "PB1": "852-KO energy panel", "PE": "PB2 + Bacformer-large probe (DL fusion)"}


def features(accs):
    X = pd.read_parquet(V / "pathway100_ko_matrix.parquet").reindex(accs).fillna(0.0)
    pm = importlib.import_module("29_pathway_modules")
    md = pd.read_csv(ROOT / "data/processed/kegg_modules_energy.tsv", sep="\t")
    M = pd.DataFrame(index=X.index, columns=[f"mod_{m}" for m in md.module], dtype=float)
    for acc in X.index:
        kset = set(X.columns[X.loc[acc] > 0])
        for m, defn in zip(md.module, md.definition):
            M.loc[acc, f"mod_{m}"] = float(np.mean([pm.step_present(st, kset) for st in pm.parse_steps(defn)]))
    bfl = pd.read_parquet(ROOT / "data/kaggle/validation178/genome_features_bacformer_large_mean.parquet")
    extra = pd.read_parquet(V / "embeddings_all.parquet") if (V / "embeddings_all.parquet").exists() else bfl
    bfl = pd.concat([bfl, extra[~extra.index.isin(bfl.index)]])
    return X.join(M.fillna(0.0)).join(bfl, how="left")


def calls(recs, accs):
    rows = []
    for r in recs:
        if r["refseq_accession"] not in accs or r.get("_source_file") == "candidates.json":
            continue
        for c, ev in r["categories"].items():
            if c in EXACT and ev.get("status") in ("tested_positive", "tested_negative"):
                rows.append(dict(accession=r["refseq_accession"], paper=str(r.get("pmid") or r["species"]),
                                 trait=EXACT[c], expected=int(ev["status"] == "tested_positive")))
    return pd.DataFrame(rows)


def mcc(y, p):
    y, p = np.asarray(y), np.asarray(p)
    tp, tn = int(((y == 1) & (p == 1)).sum()), int(((y == 0) & (p == 0)).sum())
    fp, fn = int(((y == 0) & (p == 1)).sum()), int(((y == 1) & (p == 0)).sum())
    d = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    return (tp * tn - fp * fn) / d if d else 0.0


def ranking(d, col):
    au, ap = [], []
    for _, g in d.groupby("trait"):
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
    accs = open(V / "balanced_pathway_100.txt").read().split()
    recs = json.loads((V / "validation_set.json").read_text())
    X = features(accs)
    d = calls(recs, set(accs))
    print(f"{d.accession.nunique()} organisms, {len(d)} calls, {d.paper.nunique()} papers, {d.trait.nunique()} traits")
    arms = [a for a in DESC if (ROOT / f"results/models/pathways/{a}/model.joblib").exists()]
    for arm in arms:
        bundle = joblib.load(ROOT / f"results/models/pathways/{arm}/model.joblib")
        sc, th = {}, {}
        for t, b in bundle.items():
            cols = b["model"].columns if hasattr(b["model"], "columns") else None
            Xa = X.reindex(columns=cols).fillna(0.0) if cols is not None else X
            try:
                sc[t] = pd.Series(np.asarray(b["model"].predict_proba(Xa)), index=X.index)
            except AssertionError:
                sc[t] = pd.Series(np.asarray(b["model"].predict_proba(X)), index=X.index)
            th[t] = float(b["threshold"])
        d[f"score_{arm}"] = [sc[t][a] if t in sc else np.nan for a, t in zip(d.accession, d.trait)]
        d[f"call_{arm}"] = [int(sc[t][a] >= th[t]) if t in sc else 0 for a, t in zip(d.accession, d.trait)]
    d.to_csv(V / "novel100_pathway_scores.csv", index=False)
    rows = []
    for arm in arms:
        au, ap, nt = ranking(d, f"score_{arm}")
        agree = (d[f"call_{arm}"] == d.expected).astype(int)
        rows.append(dict(model=arm, what=DESC[arm], auroc=au, auprc=ap, traits=nt, accuracy=agree.mean(),
                         mcc=mcc(d.expected, d[f"call_{arm}"]), sensitivity=agree[d.expected == 1].mean(),
                         specificity=agree[d.expected == 0].mean()))
    t = pd.DataFrame(rows).sort_values("auroc", ascending=False)
    ref = "PB3"
    for arm in arms:
        if arm != ref:
            lo, hi = boot(d, lambda x: ranking(x, f"score_{arm}")[0] - ranking(x, f"score_{ref}")[0])
            t.loc[t.model == arm, f"auroc_minus_{ref}"] = \
                f"{ranking(d, f'score_{arm}')[0] - ranking(d, f'score_{ref}')[0]:+.3f} [{lo:+.3f}, {hi:+.3f}]"
    t.to_csv(V / "novel100_pathway_comparison.csv", index=False, float_format="%.4f")
    pd.set_option("display.width", 220)
    print(t.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
