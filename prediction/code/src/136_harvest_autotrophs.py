"""Pre-registration AT (docs/preregistration_AT_autotroph_test.md): harvest further newly described species for the
autotroph test, with the procedure of src/68_harvest_candidates.py and a broader query.

Differences from src/68: the query adds autotrophy and carbon-fixation terms, the years 2010-2026 are searched, and
every genus already among the scored new species (validation_set.json, pool/) is excluded, so that the set keeps one
species per genus. Blind to any model output.
usage: python src/136_harvest_autotrophs.py [years,years] [output name]   (second round: 2000:2009 pool2.json)
Output: validation_new_organisms/pool_auto/pool.json
"""
import importlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
h = importlib.import_module("68_harvest_candidates")
OUT = ROOT / "validation_new_organisms/pool_auto"
AUTO = ('(autotroph*[tiab] OR chemolithoautotroph*[tiab] OR chemoautotroph*[tiab] OR "carbon fixation"[tiab]'
        ' OR "CO2 fixation"[tiab] OR "Calvin"[tiab] OR RuBisCO[tiab] OR "Wood-Ljungdahl"[tiab]'
        ' OR "reductive TCA"[tiab] OR methanogen*[tiab] OR acetogen*[tiab] OR "sulfur-oxidizing"[tiab]'
        ' OR "iron-oxidizing"[tiab] OR "hydrogen-oxidizing"[tiab] OR "ammonia-oxidizing"[tiab]'
        ' OR "nitrite-oxidizing"[tiab] OR thiosulfate[tiab])')


def main():
    years_list = sys.argv[1].split(",") if len(sys.argv) > 1 else ["2010:2014", "2015:2026"]
    out_name = sys.argv[2] if len(sys.argv) > 2 else "pool.json"
    OUT.mkdir(parents=True, exist_ok=True)
    excl = json.load(open(ROOT / "validation_new_organisms/exclusions.json"))
    gen_ex, sp_ex, acc_ex = set(excl["genera"]), set(excl["species"]), set(excl["accessions"])
    for f in ["validation_new_organisms/validation_set.json", "validation_new_organisms/pool/pool.json",
              "validation_new_organisms/pool_auto/pool.json"]:
        for r in json.load(open(ROOT / f)):
            gen_ex.add(r.get("genus") or r.get("species", " ").split()[0])
            acc_ex.add(r.get("refseq_accession") or r.get("accession", ""))
    pmids = []
    for years in years_list:
        q = f"nov[ti] AND ({h.PHENO} OR {AUTO}) AND {years}[dp]"
        d = h.fetch(f"{h.EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax=3000&term={h.urllib.parse.quote(q)}")
        pmids += d["esearchresult"]["idlist"]
        time.sleep(0.5)
    pmids = sorted(set(pmids))
    print(f"{len(pmids)} description papers", flush=True)
    recs = h.abstracts(pmids)
    cands, seen = [], set()
    for r in sorted(recs.values(), key=lambda r: r["pmid"], reverse=True):     # newest first
        for sp in h.species_from_title(r["title"]):
            genus = sp.split()[0]
            if genus in gen_ex or sp in sp_ex or genus in seen:
                continue
            seen.add(genus)
            cands.append(dict(species=sp, genus=genus, **r))
    print(f"{len(cands)} candidate species, one per new genus", flush=True)
    keep = []
    for i, c in enumerate(cands):
        g = h.refseq_genome(c["species"])
        if g and g["accession"] not in acc_ex and g["annotated"]:
            keep.append({**c, **g})
        if (i + 1) % 50 == 0:
            print(f"  genome lookup {i + 1}/{len(cands)}: {len(keep)} with an annotated RefSeq genome", flush=True)
        time.sleep(0.34)
    (OUT / out_name).write_text(json.dumps(keep, indent=1))
    print(f"{len(keep)} candidates with an annotated RefSeq genome -> {OUT / out_name}")


if __name__ == "__main__":
    main()
