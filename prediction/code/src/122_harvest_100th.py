"""Pathway species 100 (docs/preregistration_pathway_100th_species.md): re-run the original harvest
(src/68_harvest_candidates.py functions, unchanged) with the extended exclusions, then order the new candidates at
random (seed 20261018). Output: validation_new_organisms/pool100/{papers.json, pool.json, order.json}
"""
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
h = importlib.import_module("68_harvest_candidates")
V = ROOT / "validation_new_organisms"
OUT = V / "pool100"
OUT.mkdir(parents=True, exist_ok=True)

excl = json.load(open(V / "exclusions.json"))
gen_ex, sp_ex, acc_ex = set(excl["genera"]), set(excl["species"]), set(excl["accessions"])
for r in json.load(open(V / "pool/pool.json")) + json.load(open(V / "validation_set.json")):
    gen_ex.add(r.get("genus") or r["species"].split()[0])
    sp_ex.add(r["species"])
    acc_ex.add(r.get("accession") or r.get("refseq_accession") or "")

pmids = h.pubmed_ids("2015:2026", 1200)
print(f"{len(pmids)} description papers", flush=True)
recs = h.abstracts(pmids)
cands, seen = [], set()
for r in recs.values():
    for sp in h.species_from_title(r["title"]):
        g = sp.split()[0]
        if g in gen_ex or sp in sp_ex or sp in seen:
            continue
        seen.add(sp)
        cands.append(dict(species=sp, genus=g, **r))
print(f"{len(cands)} candidate species outside every genus already used", flush=True)
(OUT / "papers.json").write_text(json.dumps(cands, indent=1))
keep = []
for i, c in enumerate(cands):
    gnm = h.refseq_genome(c["species"])
    if gnm and gnm["accession"] not in acc_ex and gnm["annotated"]:
        keep.append({**c, **gnm})
    time.sleep(0.34)
(OUT / "pool.json").write_text(json.dumps(keep, indent=1))
print(f"{len(keep)} with an annotated RefSeq genome", flush=True)
order = [keep[i]["species"] for i in np.random.default_rng(20261018).permutation(len(keep))]
(OUT / "order.json").write_text(json.dumps(order, indent=1))
print("order:", order)
