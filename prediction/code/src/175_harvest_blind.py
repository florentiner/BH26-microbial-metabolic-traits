"""Pre-registration BL (docs/preregistration_BL_blind_test.md): harvest species described in papers that entered
PubMed from 15 September to 22 October 2026, with the query and genome lookup of src/68_harvest_candidates.py.

Every genus, species and accession anywhere in the project's data or earlier new-species pools is excluded; one
species per genus; an annotated RefSeq genome is required. Blind to any model output.
usage: python src/175_harvest_blind.py [edat range, default 2026/09/15:2026/10/22] [broad|donor]
       (BL amendment A1: 2026/01/01:2026/10/22 broad; BL2: 1800/01/01:2026/10/22 donor -> pool_blind_donor/)
Output: validation_new_organisms/pool_blind/pool.json
"""
import importlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
h = importlib.import_module("68_harvest_candidates")
V = ROOT / "validation_new_organisms"
OUT = V / "pool_blind"
BROAD = ('("nitrate"[tiab] OR anaerob*[tiab] OR "aerobic"[tiab] OR "electron acceptor"[tiab] OR "electron donor"[tiab]'
         ' OR thiosulfate[tiab] OR sulfate[tiab] OR "Fe(III)"[tiab] OR methanogen*[tiab] OR "hydrogen"[tiab])')

DONOR = ('("sulfur-oxidizing"[tiab] OR "sulphur-oxidizing"[tiab] OR thiosulfate[tiab] OR "sulfide-oxidizing"[tiab]'
         ' OR "hydrogen-oxidizing"[tiab] OR "hydrogen oxidation"[tiab] OR "Fe(II)"[tiab] OR "iron-oxidizing"[tiab]'
         ' OR "ammonia-oxidizing"[tiab] OR "nitrite-oxidizing"[tiab] OR "carbon monoxide"[tiab]'
         ' OR chemolithoautotroph*[tiab] OR chemolithotroph*[tiab] OR lithoautotroph*[tiab] OR "electron donor"[tiab])')


def main():
    edat = sys.argv[1] if len(sys.argv) > 1 else "2026/09/15:2026/10/22"
    mode = sys.argv[2] if len(sys.argv) > 2 else ""
    broad = mode == "broad"   # BL amendment A1
    donor = mode == "donor"   # pre-registration BL2
    cap = 60 if broad else 80 if donor else None
    global OUT
    if donor:
        OUT = V / "pool_blind_donor"
    OUT.mkdir(parents=True, exist_ok=True)
    excl = json.load(open(V / "exclusions.json"))
    gen_ex, sp_ex, acc_ex = set(excl["genera"]), set(excl["species"]), set(excl["accessions"])
    for f in [V / "validation_set.json", *V.glob("pool/pool.json"), *V.glob("pool_auto/pool*.json"),
              *V.glob("pool100/pool.json"), *V.glob("candidates.json"),
              *V.glob("pool_blind/pool.json")]:
        for r in json.load(open(f)):
            gen_ex.add(r.get("genus") or r.get("species", " ").split()[0])
            sp_ex.add(r.get("species", ""))
            acc_ex.add(r.get("refseq_accession") or r.get("accession", ""))
    terms = DONOR if donor else BROAD if broad else h.PHENO
    q = f"nov[ti] AND {terms} AND {edat}[edat]"
    d = h.fetch(f"{h.EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax=9999&term={h.urllib.parse.quote(q)}")
    pmids = sorted(set(d["esearchresult"]["idlist"]))
    print(f"query: {q}\n{len(pmids)} description papers", flush=True)
    recs = h.abstracts(pmids) if pmids else {}
    cands, seen = [], set()
    for r in sorted(recs.values(), key=lambda r: r["pmid"], reverse=True):
        for sp in h.species_from_title(r["title"]):
            genus = sp.split()[0]
            if genus in gen_ex or sp in sp_ex or genus in seen:
                continue
            seen.add(genus)
            cands.append(dict(species=sp, genus=genus, **r))
    print(f"{len(cands)} candidate species, one per new genus", flush=True)
    keep = []
    for c in cands:
        if cap and len(keep) >= cap:
            break
        g = h.refseq_genome(c["species"])
        if g and g["accession"] not in acc_ex and g["annotated"]:
            keep.append({**c, **g})
        time.sleep(0.34)
    (OUT / "pool.json").write_text(json.dumps(keep, indent=1))
    print(f"{len(keep)} candidates with an annotated RefSeq genome -> {OUT / 'pool.json'}")
    for k in keep:
        print(" ", k["species"], k["accession"], k.get("journal", ""), k.get("pmid", ""))


if __name__ == "__main__":
    main()
