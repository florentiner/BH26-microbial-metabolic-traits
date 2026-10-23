"""Demote evidence calls that violate the extraction rules (audit of 2026-09-21, finding 8).

Each entry below was flagged because its quote does not meet the rule it was recorded under:
stimulated growth, a biochemical test without growth coupling, fermentative reduction rather than respiration,
a medium recipe rather than a test, a non-verbatim quote, or a partly genome-derived statement. The status is set
to not_tested and the reason is appended to the quote, so the demotion is auditable rather than silent.

usage: python src/72_fix_evidence.py         (then re-run src/70_build_validation_set.py)
"""
import json
from pathlib import Path

V = Path(__file__).resolve().parents[1] / "validation_new_organisms"

DEMOTE = {
    ("GCF_009910325.1", "donor_other"): "stimulated growth is not a growth test (schema rule)",
    ("GCF_009910325.1", "acceptor_sulfur"): "stimulated growth is not a growth test (schema rule)",
    ("GCF_009910325.1", "acceptor_iron"): "stimulated growth is not a growth test (schema rule)",
    ("GCF_009910325.1", "acceptor_other"): "stimulated growth is not a growth test (schema rule)",
    ("GCF_900890425.1", "acceptor_sulfur"): "reduction during fermentation, not respiration",
    ("GCF_900890425.1", "acceptor_other"): "reduction during fermentation, not respiration",
    ("GCF_900890705.1", "acceptor_sulfur"): "reduction during fermentation, not respiration",
    ("GCF_900890705.1", "acceptor_other"): "reduction during fermentation, not respiration",
    ("GCF_001748225.1", "donor_other"): "quote is a medium recipe (yeast extract present), not a donor test",
    ("GCF_022818495.1", "acceptor_iron"): "quote is not verbatim; the outcome is the extractor's assertion",
    ("GCF_042219855.1", "acceptor_iron"): "statement is partly genomic analysis",
    ("GCF_042219855.1", "donor_other"): "statement is partly genomic analysis",
}
# Legionella maioricensis: nitrate reduction from an API-strip biochemical test, no growth coupling
LEGIONELLA = ("acceptor_nitrate", "biochemical nitrate-reduction test without growth coupling")


def main():
    n = 0
    for f in sorted(V.glob("pool/evidence_*.json")):
        recs = json.loads(f.read_text())
        for r in recs:
            acc = r.get("refseq_accession")
            for cat, ev in r.get("categories", {}).items():
                reason = DEMOTE.get((acc, cat))
                if reason is None and r.get("species", "").startswith("Legionella maioricensis") \
                        and cat == LEGIONELLA[0]:
                    reason = LEGIONELLA[1]
                if reason and ev.get("status") in ("tested_positive", "tested_negative"):
                    ev["quote"] = f"[DEMOTED by audit: {reason}] was {ev['status']}: " + str(ev.get("quote", ""))
                    ev["status"] = "not_tested"
                    n += 1
        f.write_text(json.dumps(recs, indent=1))
    print(f"demoted {n} calls to not_tested")


if __name__ == "__main__":
    main()
