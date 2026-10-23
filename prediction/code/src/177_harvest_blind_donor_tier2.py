"""Pre-registration BL2, amendment A1 (docs/preregistration_BL2_donor_blind_test.md): tier-2 candidates.

Same PubMed query as tier 1 (src/175_harvest_blind.py donor). Species and accession absent from every earlier set;
genus absent from the 358 training genomes of RedoxKO (it may occur in other project data); one species per genus;
annotated RefSeq genome required; the 80 newest candidates by PubMed ID. Blind to any model output.
Output: validation_new_organisms/pool_blind_donor/pool_tier2.json
"""
import importlib
import json
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
h = importlib.import_module("68_harvest_candidates")
hb = importlib.import_module("175_harvest_blind")
V = ROOT / "validation_new_organisms"
OUT = V / "pool_blind_donor"


def main():
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    train_genera = set(lab[lab.labelled.astype(str) == "True"].genus.dropna())
    excl = json.load(open(V / "exclusions.json"))
    sp_ex, acc_ex = set(excl["species"]), set(excl["accessions"])
    tier1 = json.load(open(OUT / "pool.json"))
    for f in [V / "validation_set.json", *V.glob("pool/pool.json"), *V.glob("pool_auto/pool*.json"),
              *V.glob("pool100/pool.json"), *V.glob("candidates.json"), *V.glob("pool_blind/pool.json")]:
        for r in json.load(open(f)):
            sp_ex.add(r.get("species", ""))
            acc_ex.add(r.get("refseq_accession") or r.get("accession", ""))
    seen = {r["genus"] for r in tier1}
    for r in tier1:
        sp_ex.add(r["species"]); acc_ex.add(r["accession"])
    q = f"nov[ti] AND {hb.DONOR} AND 1800/01/01:2026/10/22[edat]"
    d = h.fetch(f"{h.EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax=9999&term={h.urllib.parse.quote(q)}")
    pmids = sorted(set(d["esearchresult"]["idlist"]))
    recs = h.abstracts(pmids)
    cands = []
    for r in sorted(recs.values(), key=lambda r: int(r["pmid"]), reverse=True):
        for sp in h.species_from_title(r["title"]):
            genus = sp.split()[0]
            if genus in train_genera or sp in sp_ex or genus in seen:
                continue
            seen.add(genus)
            cands.append(dict(species=sp, genus=genus, tier=2, **r))
    print(f"{len(pmids)} papers, {len(cands)} tier-2 candidate species", flush=True)
    keep = []
    for c in cands:
        if len(keep) >= 80:
            break
        g = h.refseq_genome(c["species"])
        if g and g["accession"] not in acc_ex and g["annotated"]:
            keep.append({**c, **g})
        time.sleep(0.34)
    (OUT / "pool_tier2.json").write_text(json.dumps(keep, indent=1))
    print(f"{len(keep)} tier-2 species with an annotated RefSeq genome")
    for k in keep:
        print(" ", k["species"], k["accession"], k.get("pmid", ""))


if __name__ == "__main__":
    main()
