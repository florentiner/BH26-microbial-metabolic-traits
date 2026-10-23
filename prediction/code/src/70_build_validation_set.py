"""Merge the harvested evidence into the prospective validation set, and report its balance.

Input : validation_new_organisms/pool/evidence_*.json (collected blind by separate agents from the literature)
        validation_new_organisms/candidates.json      (the first, hand-curated round; merged and de-duplicated)
Output: validation_new_organisms/validation_set.json  (every organism with at least one scoreable call)
        validation_new_organisms/validation_accessions.txt
        validation_new_organisms/validation_balance.md (coverage per category, taxonomy, evidence tier)

Selection rule: keep every organism with at least one tested_positive or tested_negative. Nothing is dropped for
being inconvenient, so the set is not curated towards easy cases; balance is reported, not engineered.

usage: python src/70_build_validation_set.py
"""
import json
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
V = ROOT / "validation_new_organisms"
CATS = ["donor_sulfur_oxidation", "donor_ammonia_oxidation", "donor_nitrite_oxidation", "donor_hydrogen_oxidation",
        "donor_CO_oxidation", "donor_iron_oxidation", "donor_other", "acceptor_oxygen", "acceptor_nitrate",
        "acceptor_sulfate", "acceptor_sulfur", "acceptor_CO2_methanogen_acetogen", "acceptor_iron", "acceptor_other"]


def scoreable(rec):
    return sum(1 for c in CATS if rec.get("categories", {}).get(c, {}).get("status") in
               ("tested_positive", "tested_negative"))


def main():
    recs, seen = [], {}
    for f in sorted(V.glob("pool/evidence_*.json")) + [V / "candidates.json"]:
        if not f.exists():
            continue
        for r in json.loads(f.read_text()):
            acc = r.get("refseq_accession")
            if not acc:
                continue
            r["_source_file"] = f.name
            prev = seen.get(acc)
            # keep the record with more scoreable calls (the hand-curated round is usually richer)
            if prev is None or scoreable(r) > scoreable(prev):
                seen[acc] = r
    recs = [r for r in seen.values() if scoreable(r) > 0]
    recs.sort(key=lambda r: r["species"])
    (V / "validation_set.json").write_text(json.dumps(recs, indent=1))
    (V / "validation_accessions.txt").write_text("\n".join(r["refseq_accession"] for r in recs) + "\n")

    pos = Counter()
    neg = Counter()
    for r in recs:
        for c in CATS:
            s = r["categories"].get(c, {}).get("status")
            if s == "tested_positive":
                pos[c] += 1
            elif s == "tested_negative":
                neg[c] += 1
    tier = Counter(r.get("evidence_tier", "?")[:1] for r in recs)
    first10 = set((V / "selected_accessions.txt").read_text().split())
    rows = [f"| {c} | {pos[c]} | {neg[c]} | {pos[c] + neg[c]} |" for c in CATS]
    md = [f"# Prospective validation set — {len(recs)} organisms, "
          f"{sum(pos.values()) + sum(neg.values())} scoreable calls",
          "",
          f"Harvested from {len(list(V.glob('pool/evidence_*.json')))} blind extraction batches over 178 candidate "
          f"species, plus the {len(first10)} hand-curated organisms of the first round. Every organism with at least "
          "one published growth test is kept; none was dropped for being inconvenient.",
          "",
          f"* organisms carrying the first round: {len([r for r in recs if r['refseq_accession'] in first10])}",
          f"* evidence tier: {dict(tier)}",
          f"* organisms with both a tested positive and a tested negative: "
          f"{sum(1 for r in recs if any(r['categories'].get(c, {}).get('status') == 'tested_positive' for c in CATS) and any(r['categories'].get(c, {}).get('status') == 'tested_negative' for c in CATS))}",
          "", "## Calls per category", "",
          "| Category | Tested positive | Tested negative | Total |", "|---|---|---|---|"] + rows
    (V / "validation_balance.md").write_text("\n".join(md) + "\n")
    print(f"{len(recs)} organisms, {sum(pos.values())} tested positives, {sum(neg.values())} tested negatives")
    print(pd.DataFrame({"positive": pos, "negative": neg}).reindex(CATS).fillna(0).astype(int).to_string())


if __name__ == "__main__":
    main()
