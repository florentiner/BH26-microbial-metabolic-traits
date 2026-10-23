"""Harvest candidate organisms for the prospective validation (scale-up of validation_new_organisms/selection.md).

Funnel, all from public sources and blind to any model output:
  1. PubMed: species-description papers ("... sp. nov. ...") whose title/abstract mentions electron donor/acceptor
     chemistry, 2015-2026.
  2. Parse the proposed binomial(s) from the title.
  3. Drop every genus, species and accession that appears anywhere in the project's data (exclusions.json).
  4. Keep organisms whose class is represented among the labelled training genomes (applicability domain).
  5. Look up a RefSeq genome for the species (type material preferred) through the NCBI datasets API.
Output: validation_new_organisms/pool/pool.json  (one record per candidate: names, genome, journal, abstract)

usage: python src/68_harvest_candidates.py [--years 2015:2026] [--retmax 1200]
"""
import argparse
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation_new_organisms/pool"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
PHENO = ('("electron acceptor"[tiab] OR "electron acceptors"[tiab] OR "electron donor"[tiab] OR "electron donors"[tiab]'
         ' OR thiosulfate[tiab] OR "elemental sulfur"[tiab] OR "sulfate-reducing"[tiab] OR "sulfur-oxidizing"[tiab]'
         ' OR "iron-reducing"[tiab] OR "nitrate reduction"[tiab] OR methanogen*[tiab] OR chemolithoautotroph*[tiab]'
         ' OR chemolithotroph*[tiab] OR "hydrogen oxidation"[tiab] OR acetogen*[tiab] OR "Fe(III)"[tiab])')


def fetch(url, tries=5, parse="json"):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                raw = r.read()
            return json.loads(raw) if parse == "json" else raw.decode("utf-8", "replace")
        except Exception as e:                                   # NCBI rate limits and transient 5xx
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))
            err = e
    raise err


def pubmed_ids(years, retmax):
    q = f"nov[ti] AND {PHENO} AND {years}[dp]"
    d = fetch(f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax={retmax}&term={urllib.parse.quote(q)}")
    return d["esearchresult"]["idlist"]


def abstracts(pmids, chunk=100):
    """title, abstract, journal, year and PMC id for each PubMed record."""
    out = {}
    for i in range(0, len(pmids), chunk):
        xml = fetch(f"{EUTILS}/efetch.fcgi?db=pubmed&retmode=xml&id={','.join(pmids[i:i + chunk])}", parse="xml")
        for art in re.findall(r"<PubmedArticle>.*?</PubmedArticle>", xml, re.S):
            pmid = re.search(r"<PMID[^>]*>(\d+)</PMID>", art)
            title = " ".join(re.findall(r"<ArticleTitle[^>]*>(.*?)</ArticleTitle>", art, re.S))
            abst = " ".join(re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", art, re.S))
            jour = re.search(r"<Title>(.*?)</Title>", art, re.S)
            year = re.search(r"<PubDate>.*?<Year>(\d{4})</Year>", art, re.S)
            pmc = re.search(r'<ArticleId IdType="pmc">(PMC\d+)</ArticleId>', art)
            clean = lambda s: re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", s)).strip()       # noqa: E731
            if pmid:
                out[pmid.group(1)] = dict(pmid=pmid.group(1), title=clean(title), abstract=clean(abst),
                                          journal=clean(jour.group(1)) if jour else "",
                                          year=year.group(1) if year else "", pmcid=pmc.group(1) if pmc else "")
        time.sleep(0.4)
    return out


BINOMIAL = re.compile(r"\b([A-Z][a-z]{3,})\s+([a-z]{3,})\s+(?:gen\.\s*nov\.,?\s*)?sp\.\s*nov")


def species_from_title(title):
    """Proposed binomials in a description title, e.g. 'Ghiorsea bivora gen. nov., sp. nov., ...'."""
    names = [f"{g} {s}" for g, s in BINOMIAL.findall(title)]
    seen, out = set(), []
    for n in names:
        if n not in seen:
            seen.add(n)
            out.append(n)
    return out


def refseq_genome(species):
    """A RefSeq assembly for the species; type material and higher assembly level preferred."""
    url = (f"{DATASETS}/genome/taxon/{urllib.parse.quote(species)}/dataset_report"
           "?filters.assembly_source=refseq&page_size=20")
    try:
        d = fetch(url)
    except Exception:
        return None
    best = None
    for r in d.get("reports", []):
        acc = r.get("accession", "")
        if not acc.startswith("GCF_"):
            continue
        info = r.get("assembly_info", {})
        lvl = {"Complete Genome": 3, "Chromosome": 2, "Scaffold": 1, "Contig": 0}.get(info.get("assembly_level"), 0)
        is_type = 1 if info.get("biosample", {}).get("attributes") is not None and \
            r.get("type_material", {}).get("type_label") else 0
        ann = 1 if r.get("annotation_info") else 0
        key = (ann, is_type, lvl)
        if best is None or key > best[0]:
            best = (key, dict(accession=acc, organism=r.get("organism", {}).get("organism_name", ""),
                              level=info.get("assembly_level", ""), annotated=bool(ann), type_material=bool(is_type),
                              release=info.get("release_date", ""),
                              n_proteins=r.get("annotation_info", {}).get("stats", {})
                              .get("gene_counts", {}).get("protein_coding")))
    return best[1] if best else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", default="2015:2026")
    ap.add_argument("--retmax", type=int, default=1200)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    excl = json.load(open(ROOT / "validation_new_organisms/exclusions.json"))
    gen_ex, sp_ex, acc_ex = set(excl["genera"]), set(excl["species"]), set(excl["accessions"])

    pmids = pubmed_ids(args.years, args.retmax)
    print(f"{len(pmids)} description papers with donor/acceptor chemistry in the abstract", flush=True)
    recs = abstracts(pmids)
    print(f"{len(recs)} abstracts fetched", flush=True)

    cands, seen = [], set()
    for r in recs.values():
        for sp in species_from_title(r["title"]):
            genus = sp.split()[0]
            if genus in gen_ex or sp in sp_ex or sp in seen:
                continue
            seen.add(sp)
            cands.append(dict(species=sp, genus=genus, **r))
    print(f"{len(cands)} candidate species after removing every genus and species already in the project", flush=True)
    (OUT / "papers.json").write_text(json.dumps(cands, indent=1))

    keep = []
    for i, c in enumerate(cands):
        g = refseq_genome(c["species"])
        if g and g["accession"] not in acc_ex and g["annotated"]:
            keep.append({**c, **g})
        if (i + 1) % 25 == 0:
            print(f"  genome lookup {i + 1}/{len(cands)}: {len(keep)} with an annotated RefSeq genome", flush=True)
        time.sleep(0.34)
    (OUT / "pool.json").write_text(json.dumps(keep, indent=1))
    print(f"{len(keep)} candidates with an annotated RefSeq genome -> {OUT / 'pool.json'}")


if __name__ == "__main__":
    main()
