"""Score the released models on the prospective validation set (validation_new_organisms/selection_v2.md).

Revised after the audit of 2026-09-21. What changed and why:
* PRIMARY analysis uses only the organisms whose evidence came from the six blind extraction batches. The 24
  organisms carried over from the hand-curated first round are reported separately: their records have no verbatim
  quote, and model output for ten of them existed before this protocol was written.
* Accuracy alone is not reported as the headline. A validation set built by searching for papers that advertise
  redox chemistry is ~60% positive, so three baselines are scored on exactly the same calls: the constant rule a
  deployer could actually apply (training prevalence), an oracle per-category majority, and always-positive.
  MCC and balanced accuracy are reported beside accuracy.
* "Deployed mode" uses each arm's OWN reliable-trait set, recomputed with the released rule from
  src/14_predict_unlabelled.py (mcc_ci_low > 0.2, BH q < 0.05, delta_vs_null > 0, not contradicted externally),
  at the split level that applies to each organism. Hardcoding one list for both arms reversed their order.
* The two catch-all categories (donor_other, acceptor_other) are excluded from the headline: the training labels
  come from the source database's literal "Other" column, the validation schema defines them by substrate lists,
  and the two do not denote the same thing.
* Confidence intervals come from a cluster bootstrap over description papers (organisms sharing a paper are not
  independent), not from Wilson intervals that assume independent calls.

usage: python src/71_score_validation.py predictions178_B2.csv predictions178_B2XLE_large.csv
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
V = ROOT / "validation_new_organisms"
RES = ROOT / "results"
EXPECT = {"tested_positive": 1, "tested_negative": 0}
CATCH_ALL = {"donor_other", "acceptor_other"}


def reliable_traits(arm):
    """The trait sets the released predictor actually emits calls for, per split level (src/14_predict_unlabelled)."""
    mt = pd.read_csv(RES / "metrics_per_trait.csv")
    mt = mt[mt.arm == arm]
    out = {}
    for level in ("genus", "family"):
        d = mt[mt.split_level == level]
        ok = set(d.loc[(d.mcc_ci_low > 0.2) & (d.delta_vs_null_q_bh < 0.05) & (d.delta_vs_null > 0), "trait"])
        ext = RES / "external_validation_per_trait.csv"
        if ext.exists():
            ev = pd.read_csv(ext)
            ev = ev[(ev.subset == "all") & (ev.arm == arm) & (ev.evaluable == True)]          # noqa: E712
            ok -= set(ev.loc[ev.auprc_lift <= 1, "trait"])
        out[level] = ok
    return out


def paper_id(rec):
    for k in ("pmid", "description_reference", "description_url", "species"):
        if rec.get(k):
            return str(rec[k])
    return rec["refseq_accession"]


def build(recs, pred):
    rows = []
    for r in recs:
        acc = r["refseq_accession"]
        for cat, ev in r["categories"].items():
            exp = EXPECT.get(ev.get("status"))
            if exp is None:
                continue
            q = pred[(pred.genome_id == acc) & (pred.trait == cat)]
            if q.empty:
                continue
            q = q.iloc[0]
            rows.append(dict(accession=acc, organism=r["species"], paper=paper_id(r), category=cat,
                             blind=r.get("_source_file", "") != "candidates.json",
                             published=ev["status"], expected=exp, call=int(q.predicted),
                             score=round(float(q.score), 4), threshold=round(float(q.threshold), 4),
                             family_in_training=bool(q.get("family_in_training_set", False)),
                             agree=int(exp == int(q.predicted))))
    return pd.DataFrame(rows)


def mcc(y, p):
    y, p = np.asarray(y), np.asarray(p)
    tp, tn = int(((y == 1) & (p == 1)).sum()), int(((y == 0) & (p == 0)).sum())
    fp, fn = int(((y == 0) & (p == 1)).sum()), int(((y == 1) & (p == 0)).sum())
    den = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    return (tp * tn - fp * fn) / den if den else 0.0


def bacc(y, p):
    y, p = np.asarray(y), np.asarray(p)
    se = (p[y == 1] == 1).mean() if (y == 1).any() else np.nan
    sp = (p[y == 0] == 0).mean() if (y == 0).any() else np.nan
    return np.nanmean([se, sp])


def boot_ci(d, stat, n_boot=2000, seed=20260921):
    """Cluster bootstrap over description papers."""
    rng = np.random.default_rng(seed)
    papers = d.paper.unique()
    idx = {p: np.where(d.paper.to_numpy() == p)[0] for p in papers}
    vals = []
    for _ in range(n_boot):
        take = np.concatenate([idx[p] for p in rng.choice(papers, len(papers), replace=True)])
        vals.append(stat(d.iloc[take]))
    return np.percentile(vals, 2.5), np.percentile(vals, 97.5)


def report(d, label, train_prev):
    if d.empty:
        print(f"\n== {label}: no calls")
        return None
    acc = lambda x: x.agree.mean()                                                   # noqa: E731
    lo, hi = boot_ci(d, acc)
    m = mcc(d.expected, d.call)
    b = bacc(d.expected, d.call)
    # baselines on the same calls
    prior = d.category.map(lambda c: int(train_prev.get(c, 0) >= 0.5))
    maj = d.groupby("category").expected.transform(lambda s: int(s.mean() >= 0.5))
    print(f"\n== {label}: {len(d)} calls, {d.accession.nunique()} organisms, {d.paper.nunique()} papers")
    print(f"   accuracy {acc(d):.1%} [{lo:.1%}, {hi:.1%}] | MCC {m:.3f} | balanced accuracy {b:.1%}")
    print(f"   sensitivity {d[d.expected == 1].agree.mean():.1%} ({int(d[d.expected == 1].agree.sum())}/"
          f"{int((d.expected == 1).sum())}), specificity {d[d.expected == 0].agree.mean():.1%} "
          f"({int(d[d.expected == 0].agree.sum())}/{int((d.expected == 0).sum())})")
    print(f"   baselines on the same calls: training-prevalence rule {(prior == d.expected).mean():.1%} "
          f"(MCC {mcc(d.expected, prior):.3f}) | oracle per-category majority {(maj == d.expected).mean():.1%} "
          f"(MCC {mcc(d.expected, maj):.3f}) | always-positive {(d.expected == 1).mean():.1%}")
    return dict(label=label, n=len(d), organisms=d.accession.nunique(), papers=d.paper.nunique(),
                accuracy=acc(d), ci_low=lo, ci_high=hi, mcc=m, balanced_accuracy=b,
                sensitivity=d[d.expected == 1].agree.mean(), specificity=d[d.expected == 0].agree.mean(),
                baseline_train_prior=(prior == d.expected).mean(), baseline_majority=(maj == d.expected).mean())


def main():
    recs = json.loads((V / "validation_set.json").read_text())
    subset_file = [a for a in sys.argv[1:] if a.startswith("--subset=")]
    if subset_file:
        keep = set(Path(subset_file[0].split("=", 1)[1]).read_text().split())
        recs = [r for r in recs if r["refseq_accession"] in keep]
        print(f"restricted to the balanced subset: {len(recs)} organisms")
    tag = "_balanced100" if subset_file else ""
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    traits = [c for c in lab.columns if c.startswith(("donor_", "acceptor_"))]
    train_prev = {t: lab[t].mean() for t in traits}
    summary, scored = [], {}
    for f in [a for a in sys.argv[1:] if not a.startswith("--")]:
        pred = pd.read_csv(V / f)
        arm = pred.arm.iloc[0]
        rel = reliable_traits(arm)
        d = build(recs, pred)
        missing = {r["refseq_accession"] for r in recs} - set(d.accession)
        d["deployed"] = [c in rel["genus" if fam else "family"] for c, fam in zip(d.category, d.family_in_training)]
        d.to_csv(V / f"scored178_{arm}{tag}.csv", index=False)
        print(f"\n######## {arm}  (reliable traits: genus-level {sorted(rel['genus'])}, "
              f"family-level {sorted(rel['family'])})")
        print(f"   {len(missing)} selected organisms have no prediction and are excluded: {sorted(missing)}")
        blind = d[d.blind]
        core = blind[~blind.category.isin(CATCH_ALL)]
        for sub, name in [(core, "PRIMARY: blind evidence, excluding the two catch-all categories"),
                          (blind, "blind evidence, all categories"),
                          (d, "all organisms including the non-blind first round"),
                          (d[d.deployed], "deployed mode: each arm's own reliable traits")]:
            r = report(sub, name, train_prev)
            if r:
                summary.append(dict(arm=arm, **r))
        per = []
        for c, g in d.groupby("category"):
            gp, gn = g[g.expected == 1], g[g.expected == 0]
            per.append(dict(category=c, n=len(g), accuracy=round(g.agree.mean(), 3), mcc=round(mcc(g.expected, g.call), 3),
                            sensitivity=round(gp.agree.mean(), 3) if len(gp) else np.nan,
                            specificity=round(gn.agree.mean(), 3) if len(gn) else np.nan,
                            n_pos=len(gp), n_neg=len(gn), deployed=bool(g.deployed.any())))
        pc = pd.DataFrame(per).sort_values("n", ascending=False)
        pc.to_csv(V / f"scored178_{arm}{tag}_per_category.csv", index=False)
        print(pc.to_string(index=False))
        scored[arm] = d
    pd.DataFrame(summary).to_csv(V / f"validation_summary{tag}.csv", index=False, float_format="%.4f")
    if len(scored) == 2:
        a, b = list(scored)
        m = scored[a].merge(scored[b], on=["accession", "category"], suffixes=("_a", "_b"))
        m = m[m.blind_a & ~m.category.isin(CATCH_ALL)]
        rng = np.random.default_rng(20260921)
        papers = m.paper_a.unique()
        idx = {p: np.where(m.paper_a.to_numpy() == p)[0] for p in papers}
        delta = []
        for _ in range(2000):
            take = np.concatenate([idx[p] for p in rng.choice(papers, len(papers), replace=True)])
            delta.append(m.agree_a.to_numpy()[take].mean() - m.agree_b.to_numpy()[take].mean())
        nd_a = int(((m.agree_a == 1) & (m.agree_b == 0)).sum())
        nd_b = int(((m.agree_a == 0) & (m.agree_b == 1)).sum())
        print(f"\n{a} − {b} on the primary set: {m.agree_a.mean() - m.agree_b.mean():+.3f} "
              f"[{np.percentile(delta, 2.5):+.3f}, {np.percentile(delta, 97.5):+.3f}] "
              f"(paper-clustered bootstrap, {len(m)} paired calls; discordant {nd_a} vs {nd_b})")


if __name__ == "__main__":
    main()
