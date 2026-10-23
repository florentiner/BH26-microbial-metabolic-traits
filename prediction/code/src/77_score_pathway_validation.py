"""Score the released pathway model on the prospective organisms (preregistration_pathways.md, external addendum).

The mapping from published donor/acceptor growth tests to pathway traits was fixed and committed before any
pathway prediction existed for these genomes. Seven mappings are exact; two are approximate and reported apart.

usage: python src/77_score_pathway_validation.py <pathway_predictions.csv>
Output: validation_new_organisms/pathway_scored.csv, pathway_summary.csv and a printed report
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
V = ROOT / "validation_new_organisms"
EXACT = {"acceptor_sulfate": "sulfate_reduction", "acceptor_sulfur": "sulfur_reduction",
         "donor_sulfur_oxidation": "sulfur_oxidation", "acceptor_nitrate": "nitrate_reduction",
         "acceptor_iron": "iron_reduction", "donor_iron_oxidation": "iron_oxidation",
         "donor_hydrogen_oxidation": "hydrogen_oxidation"}
APPROX = {"acceptor_oxygen": "aerobic_chemoheterotrophy",
          "acceptor_CO2_methanogen_acetogen": "methanogenesis|acetogenesis"}
EXPECT = {"tested_positive": 1, "tested_negative": 0}


def mcc(y, p):
    y, p = np.asarray(y), np.asarray(p)
    tp, tn = int(((y == 1) & (p == 1)).sum()), int(((y == 0) & (p == 0)).sum())
    fp, fn = int(((y == 0) & (p == 1)).sum()), int(((y == 1) & (p == 0)).sum())
    den = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    return (tp * tn - fp * fn) / den if den else 0.0


def boot_ci(d, stat, n_boot=2000, seed=20260921):
    rng = np.random.default_rng(seed)
    papers = d.paper.unique()
    idx = {p: np.where(d.paper.to_numpy() == p)[0] for p in papers}
    vals = [stat(d.iloc[np.concatenate([idx[p] for p in rng.choice(papers, len(papers), replace=True)])])
            for _ in range(n_boot)]
    return np.percentile(vals, 2.5), np.percentile(vals, 97.5)


def report(d, label, prev):
    if d.empty:
        return None
    lo, hi = boot_ci(d, lambda x: x.agree.mean())
    prior = d.trait.map(lambda t: int(prev.get(t, 0) >= 0.5))
    maj = d.groupby("trait").expected.transform(lambda s: int(s.mean() >= 0.5))
    pos, neg = d[d.expected == 1], d[d.expected == 0]
    print(f"\n== {label}: {len(d)} calls, {d.accession.nunique()} organisms, {d.paper.nunique()} papers")
    print(f"   accuracy {d.agree.mean():.1%} [{lo:.1%}, {hi:.1%}] | MCC {mcc(d.expected, d.call):.3f} | "
          f"sensitivity {pos.agree.mean():.1%} ({int(pos.agree.sum())}/{len(pos)}) | "
          f"specificity {neg.agree.mean():.1%} ({int(neg.agree.sum())}/{len(neg)})")
    print(f"   baselines: training-prevalence rule {(prior == d.expected).mean():.1%} "
          f"(MCC {mcc(d.expected, prior):.3f}) | oracle majority {(maj == d.expected).mean():.1%}")
    return dict(label=label, n=len(d), organisms=d.accession.nunique(), accuracy=d.agree.mean(), ci_low=lo,
                ci_high=hi, mcc=mcc(d.expected, d.call), sensitivity=pos.agree.mean(),
                specificity=neg.agree.mean(), baseline_train_prior=(prior == d.expected).mean())


def main():
    pred = pd.read_csv(sys.argv[1])
    recs = json.loads((V / "validation_set.json").read_text())
    lab = pd.read_csv(ROOT / "data/pathways/labels.tsv", sep="\t")
    traits = [c for c in lab.columns if lab[c].dropna().isin([0, 1]).all() and lab[c].sum() > 0
              and c not in ("species_tax_id",)]
    prev = {t: lab[t].mean() for t in traits if t in lab}
    rel = json.loads((ROOT / "results/models/pathways/RELEASED.json").read_text())
    reliable = set(pd.read_csv(V / "pathway_reliable_traits.csv").trait) \
        if (V / "pathway_reliable_traits.csv").exists() else set()
    rows = []
    for r in recs:
        acc = r["refseq_accession"]
        for cat, ev in r["categories"].items():
            exp = EXPECT.get(ev.get("status"))
            if exp is None or (cat not in EXACT and cat not in APPROX):
                continue
            t = EXACT.get(cat, APPROX.get(cat))
            sub = pred[(pred.genome_id == acc) & (pred.trait.isin(t.split("|")))]
            if sub.empty:
                continue
            best = sub.loc[sub.score.idxmax()]                 # max over the alternatives for the combined trait
            rows.append(dict(accession=acc, organism=r["species"], paper=str(r.get("pmid", r["species"])),
                             evidence_category=cat, trait=t, exact=cat in EXACT,
                             blind=r.get("_source_file", "") != "candidates.json", published=ev["status"],
                             expected=exp, call=int(best.predicted), score=round(float(best.score), 4),
                             agree=int(exp == int(best.predicted)), deployed=t in reliable))
    d = pd.DataFrame(rows)
    d.to_csv(V / "pathway_scored.csv", index=False)
    print(f"released pathway model: {rel['arm']}")
    out = []
    for sub, name in [(d[d.blind & d.exact], "PRIMARY: blind evidence, exact mappings only"),
                      (d[d.exact], "exact mappings, all organisms"),
                      (d[~d.exact], "approximate mappings (reported apart)"),
                      (d, "everything")]:
        r = report(sub, name, prev)
        if r:
            out.append(r)
    pd.DataFrame(out).to_csv(V / "pathway_summary.csv", index=False, float_format="%.4f")
    per = []
    for t, g in d.groupby("trait"):
        gp, gn = g[g.expected == 1], g[g.expected == 0]
        per.append(dict(trait=t, n=len(g), accuracy=round(g.agree.mean(), 3), mcc=round(mcc(g.expected, g.call), 3),
                        sensitivity=round(gp.agree.mean(), 3) if len(gp) else np.nan,
                        specificity=round(gn.agree.mean(), 3) if len(gn) else np.nan, n_pos=len(gp), n_neg=len(gn)))
    pc = pd.DataFrame(per).sort_values("n", ascending=False)
    pc.to_csv(V / "pathway_per_trait.csv", index=False)
    print("\nper trait:")
    print(pc.to_string(index=False))


if __name__ == "__main__":
    main()
