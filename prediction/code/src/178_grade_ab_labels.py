"""Review round 3: do the cross-validation results depend on positives that may rest on genome-based inference?

Each training positive is traced to the donor or acceptor cell of the audited phenotype collection
(phenotype_data_GRADED_VALIDATED.csv, matched by binomial, categories assigned with the rules of
src/92_label_audit_impact.py). The cell grade is the lowest grade across its entries, so a positive counts as
growth-backed (grade A/B) only when every entry of its cell is graded A or B. Positives without such support are set
to missing, negatives are kept, and the saved out-of-fold predictions of RedoxKO (RF) and the taxonomic baseline are
re-scored at genus and family level. Paired difference with a cluster bootstrap over genera or families.
Output: results/additional/rf_paper/grade_ab.json
"""
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
m92 = importlib.import_module("92_label_audit_impact")
OUT = ROOT / "results/additional/rf_paper"
SEED, B = 20261022, 2000


def support():
    p = pd.read_csv(m92.AFTER, low_memory=False)
    p["binomial"] = p.organism.astype(str).str.replace(r"\s+sp\. nov\.?", "", regex=True) \
        .str.replace(r"\(.*?\)", "", regex=True).str.split().str[:2].str.join(" ")
    ab = {}
    for b, d in p.groupby("binomial"):
        for _, r in d.iterrows():
            for col, rules, other, gcol in (("electron_donors", m92.DONOR, "donor_other", "electron_donors_confidence_grade"),
                                            ("electron_acceptors", m92.ACCEPTOR, "acceptor_other",
                                             "electron_acceptors_confidence_grade")):
                g = r.get(gcol)
                for c in m92.categorise(r.get(col), rules, other):
                    ab[(b, c)] = ab.get((b, c), False) or (g in ("A", "B"))
    return ab


def macro(d):
    v = [average_precision_score(g.y, g.prob) for _, g in d.groupby("trait") if 0 < g.y.sum() < len(g)]
    return float(np.mean(v))


def main():
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    lab = lab[lab.labelled == True].copy()  # noqa: E712
    lab["binomial"] = lab.species.fillna(lab.organism_name).astype(str).str.split().str[:2].str.join(" ")
    bino = lab.set_index("genome_id").binomial
    ab = support()
    rf = pd.read_parquet(ROOT / "results/trees_xai/da/oof/RF.parquet")
    tx = pd.read_parquet(ROOT / "results/oof/nulls.parquet").query("arm == 'null_taxNN'")
    tax = lab.set_index("genome_id")
    res = {"cells": {}}
    keep_pos = {}
    for (g, t), y in rf.drop_duplicates(["genome_id", "trait"]).set_index(["genome_id", "trait"]).y.items():
        if y == 1:
            keep_pos[(g, t)] = bool(ab.get((bino.get(g), t), False))
    kp = pd.Series(keep_pos)
    per = kp.groupby(level=1).agg(["size", "sum"]).rename(columns={"size": "positives", "sum": "grade_AB"})
    res["cells"] = per.astype(int).to_dict("index")
    res["positives"], res["positives_AB"] = int(kp.size), int(kp.sum())
    for level, rank in (("genus", "genus"), ("family", "family")):
        a = rf[rf.level == level].copy()
        b = tx[tx.level == level].copy()
        drop = {k for k, v in keep_pos.items() if not v}
        for d in (a, b):
            d["keep"] = [not ((g, t) in drop) for g, t in zip(d.genome_id, d.trait)]
        a2, b2 = a[a.keep], b[b.keep]
        r = {"RedoxKO_all": macro(a), "taxonomic_all": macro(b), "RedoxKO_AB": macro(a2), "taxonomic_AB": macro(b2)}
        # paired bootstrap of (RedoxKO - taxonomic) on the A/B-restricted labels, clusters = taxon of the level
        cl = tax[rank].reindex(a2.genome_id).to_numpy()
        a2 = a2.assign(cl=cl)
        b2 = b2.merge(a2[["genome_id", "trait", "cl"]], on=["genome_id", "trait"])
        rng = np.random.default_rng(SEED)
        cls = np.unique(cl.astype(str))
        ia = {c: np.where(a2.cl.astype(str).to_numpy() == c)[0] for c in cls}
        ib = {c: np.where(b2.cl.astype(str).to_numpy() == c)[0] for c in cls}
        diffs = []
        for _ in range(B):
            s = rng.choice(cls, len(cls))
            xa = a2.iloc[np.concatenate([ia[c] for c in s])]
            xb = b2.iloc[np.concatenate([ib[c] for c in s])]
            try:
                diffs.append(macro(xa) - macro(xb))
            except Exception:
                pass
        r["diff_AB"] = r["RedoxKO_AB"] - r["taxonomic_AB"]
        r["diff_AB_ci"] = np.nanpercentile(diffs, [2.5, 97.5]).tolist()
        res[level] = r
    (OUT / "grade_ab.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: v for k, v in res.items() if k != "cells"}, indent=1, default=lambda v: round(float(v), 3)))
    print(per.to_string())


if __name__ == "__main__":
    main()
