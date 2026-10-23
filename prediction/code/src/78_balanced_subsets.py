"""Balanced 100-genome subsets of the evaluation sets, selected from labels only.

Why: the full evaluation sets are uneven. The 562-genome external set ranges from 19 positives (iron oxidation) to
289 (oxygen); the prospective set ranges from 5 scored calls (CO oxidation) to 102 (oxygen). A balanced subset
makes per-category comparison readable and shows the headline numbers are not carried by the common categories.

How the subset is chosen, so it cannot be cherry-picking:
* labels only - model scores, predictions and per-genome performance are never read;
* greedy maximisation of the minimum per-category count, alternating positives and negatives, so rare categories
  fill first;
* at most one genome per genus, which also forces taxonomic spread;
* deterministic (fixed seed) and written to disk, so the same subset is reproduced by anyone.

The subset is reported ALONGSIDE the full set, never instead of it.

usage: python src/78_balanced_subsets.py {external,prospective,pathway} [--n 100]
Output: validation_new_organisms/balanced_<which>_100.txt  (+ a coverage table printed)
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
V = ROOT / "validation_new_organisms"
RES = ROOT / "results"
SEED = 20260921


def greedy_balanced(items, n, seed=SEED):
    """items: {key: ({trait: 1/0}, genus)}. Maximise the minimum per-(trait, class) count, one genome per genus."""
    rng = np.random.default_rng(seed)
    keys = sorted(items)
    rng.shuffle(keys)
    chosen, used_genera, counts = [], set(), {}
    cells = sorted({(t, c) for k in keys for t, c in items[k][0].items()})
    for t, c in cells:
        counts[(t, c)] = 0
    while len(chosen) < n:
        best, best_gain = None, None
        for k in keys:
            labels, genus = items[k]
            if k in chosen or (genus and genus in used_genera):
                continue
            # gain = how much this genome lifts the currently scarcest cells it covers
            gain = sorted(counts[(t, c)] for t, c in labels.items())
            gain = (len(labels) and gain[0], -len(labels))     # scarcest cell first, then prefer more coverage
            if best is None or gain < best_gain:
                best, best_gain = k, gain
        if best is None:
            break
        labels, genus = items[best]
        for t, c in labels.items():
            counts[(t, c)] += 1
        chosen.append(best)
        if genus:
            used_genera.add(genus)
    return chosen, counts


def external(n):
    p = pd.read_csv(RES / "external_validation_predictions.csv")
    p = p[p.arm == "B2"]
    ev = pd.read_csv(RES / "external_validation_per_trait.csv")
    traits = set(ev[(ev.arm == "B2") & (ev.subset == "all") & ev.evaluable].trait)
    p = p[p.trait.isin(traits) & p.y.notna()]
    lab = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t").set_index("genome_id")
    items = {}
    for g, d in p.groupby("genome_id"):
        items[g] = ({t: int(y) for t, y in zip(d.trait, d.y)}, str(lab.genus.get(g, "")) or None)
    return items, "external"


def prospective(n):
    recs = json.loads((V / "validation_set.json").read_text())
    items = {}
    for r in recs:
        lab = {c: int(ev["status"] == "tested_positive") for c, ev in r["categories"].items()
               if ev.get("status") in ("tested_positive", "tested_negative")}
        if lab:
            items[r["refseq_accession"]] = (lab, r.get("genus") or r["species"].split()[0])
    return items, "prospective"


def prospective_primary(n):
    """The headline set: blind evidence only, scored on the categories the primary analysis uses."""
    recs = json.loads((V / "validation_set.json").read_text())
    items = {}
    for r in recs:
        if r.get("_source_file", "") == "candidates.json":          # non-blind first round
            continue
        lab = {c: int(ev["status"] == "tested_positive") for c, ev in r["categories"].items()
               if ev.get("status") in ("tested_positive", "tested_negative")
               and c not in ("donor_other", "acceptor_other")}
        if lab:
            items[r["refseq_accession"]] = (lab, r.get("genus") or r["species"].split()[0])
    return items, "prospective_primary"


EXACT_PW = {"acceptor_sulfate": "sulfate_reduction", "acceptor_sulfur": "sulfur_reduction",
            "donor_sulfur_oxidation": "sulfur_oxidation", "acceptor_nitrate": "nitrate_reduction",
            "acceptor_iron": "iron_reduction", "donor_iron_oxidation": "iron_oxidation",
            "donor_hydrogen_oxidation": "hydrogen_oxidation"}


def pathway_labels(n):
    """Task-2 balanced set from labels alone (no predictions needed): blind evidence, exact mappings only."""
    recs = json.loads((V / "validation_set.json").read_text())
    items = {}
    for r in recs:
        if r.get("_source_file", "") == "candidates.json":
            continue
        lab = {EXACT_PW[c]: int(ev["status"] == "tested_positive") for c, ev in r["categories"].items()
               if c in EXACT_PW and ev.get("status") in ("tested_positive", "tested_negative")}
        if lab:
            items[r["refseq_accession"]] = (lab, r.get("genus") or r["species"].split()[0])
    return items, "pathway"


def pathway(n):
    f = V / "pathway_scored.csv"
    if not f.exists():
        raise SystemExit("run src/77_score_pathway_validation.py first")
    d = pd.read_csv(f)
    d = d[d.exact]
    items = {}
    for g, s in d.groupby("accession"):
        items[g] = ({t: int(e) for t, e in zip(s.trait, s.expected)}, s.organism.iloc[0].split()[0])
    return items, "pathway"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("which", choices=["external", "prospective", "prospective_primary", "pathway", "pathway_labels"])
    ap.add_argument("--n", type=int, default=100)
    args = ap.parse_args()
    items, name = {"external": external, "prospective": prospective, "prospective_primary": prospective_primary, "pathway": pathway, "pathway_labels": pathway_labels}[args.which](args.n)
    print(f"{name}: {len(items)} candidate genomes, {len({g for _, g in items.values() if g})} genera")
    chosen, counts = greedy_balanced(items, args.n)
    out = V / f"balanced_{name}_{args.n}.txt"
    out.write_text("\n".join(sorted(chosen)) + "\n")
    rows = {}
    for (t, c), v in counts.items():
        rows.setdefault(t, {})[("positives" if c == 1 else "negatives")] = v
    tab = pd.DataFrame(rows).T.fillna(0).astype(int).sort_index()
    tab["total"] = tab.sum(axis=1)
    print(f"{len(chosen)} genomes selected -> {out}")
    print(tab.to_string())


if __name__ == "__main__":
    main()
