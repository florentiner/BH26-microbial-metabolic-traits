"""Build the pathway-prediction dataset (docs/preregistration_pathways.md).

Species with literature-curated pathway annotations (Madin et al. 2020) -> one RefSeq genome each (NCBI Datasets)
-> full significant-KO profile from the RDF Portal BH26 store (the resource TogoMCP exposes). No genome downloads
and no local HMM search: the KO profiles come from the same RDF resource the project is built on, and our local
KofamScan reproduces them to 1 cell in 43 344 (results/qc/togo_concordance.json).

Steps (each cached, safe to re-run):
  1. species -> RefSeq accession            data/pathways/species_accessions.tsv
  2. accession -> BH26 KO profile           data/pathways/togo_kofam.parquet, ko_matrix.tsv
  3. labels, taxonomy, frozen splits        data/pathways/labels.tsv, splits.tsv (+ .sha256)
usage: python src/30_pathway_dataset.py [--step 1|2|3|all] [--max-species N]
"""
import hashlib
import importlib
import json
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/pathways"
LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
SEED = 20260918
N_FOLDS = 5
MIN_POS = 25

TRAITS = {
    "denitrification": ["denitrification", "nitrous_oxide_reduction", "nitric_oxide_reduction"],
    "dnra": ["nitrate_reduction_to_ammonia", "nitrite_reduction_to_ammonia"],
    "nitrate_reduction": ["nitrate_reduction", "nitrite_reduction", "nitrate_reduction_to_ammonia",
                          "nitrite_reduction_to_ammonia", "denitrification"],
    "nitrogen_fixation": ["nitrogen_fixation"],
    "nitrification": ["nitrification", "ammonia_oxidation"],
    "methanogenesis": ["methanogenesis"],
    "acetogenesis": ["acetogenesis"],
    "sulfate_reduction": ["sulfate_reduction"],
    "sulfur_reduction": ["sulfur_reduction", "thiosulfate_reduction", "sulfite_reduction"],
    "sulfur_oxidation": ["sulfur_oxidation", "thiosulfate_oxidation", "sulfide_oxidation", "sulfite_oxidation",
                         "anoxygenic_photoautotrophy_sulfur_oxidation"],
    "iron_reduction": ["iron_reduction"],
    "iron_oxidation": ["iron_oxidation"],
    "manganese_reduction": ["manganese_reduction"],
    "arsenate_reduction": ["arsenate_reduction"],
    "fumarate_reduction": ["fumarate_reduction"],
    "hydrogen_oxidation": ["hydrogen_oxidation", "anoxygenic_photoautotrophy_hydrogen_oxidation"],
    "fermentation": ["fermentation"],
    "aerobic_chemoheterotrophy": ["aerobic_chemo_heterotrophy", "aerobic_heterotrophy"],
    "photoautotrophy": ["photoautotrophy", "anoxygenic_photoautotrophy", "oxygenic_photoautotrophy"],
    "photoheterotrophy": ["photoheterotrophy", "aerobic_anoxygenic_phototrophy"],
    "methylotrophy": ["methylotrophy", "methanol_oxidation", "methylamine_oxidation"],
    "methanotrophy": ["methane_oxidation"],
}


def madin():
    d = pd.read_csv(ROOT / "data/raw/madin/condensed_species_NCBI.csv", low_memory=False)
    return d[d.pathways.notna()].copy()


def best_accession(taxid):
    url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{int(taxid)}/dataset_report"
           "?filters.assembly_source=refseq&filters.has_annotation=true&page_size=100")
    for i in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                reps = json.load(r).get("reports", [])
            break
        except Exception:
            time.sleep(4 * (i + 1))
    else:
        return None
    best = None
    for rep in reps:
        a = rep.get("assembly_info", {})
        key = (0 if a.get("refseq_category") == "reference genome" else 1,
               LEVEL_RANK.get(a.get("assembly_level"), 9),
               -int((a.get("release_date") or "0000-00-00").replace("-", "")[:8] or 0))
        if best is None or key < best[0]:
            best = (key, rep["accession"])
    time.sleep(0.2)
    return best[1] if best else None


def step1(max_species):
    d = madin()
    if max_species:
        d = d.head(max_species)
    OUT.mkdir(parents=True, exist_ok=True)
    f = OUT / "species_accessions.tsv"
    done = pd.read_csv(f, sep="\t") if f.exists() else pd.DataFrame(columns=["species_tax_id", "genome_id"])
    todo = d[~d.species_tax_id.isin(done.species_tax_id)]
    print(f"species with pathways: {len(d)}; accessions already mapped: {len(done)}; to map: {len(todo)}", flush=True)
    rows = []
    with ThreadPoolExecutor(4) as ex:
        for k, (tid, acc) in enumerate(zip(todo.species_tax_id, ex.map(best_accession, todo.species_tax_id))):
            rows.append((tid, acc))
            if k % 200 == 0:
                print(f"  {k}/{len(todo)} mapped", flush=True)
                pd.concat([done, pd.DataFrame(rows, columns=["species_tax_id", "genome_id"])]).to_csv(f, sep="\t",
                                                                                                     index=False)
    out = pd.concat([done, pd.DataFrame(rows, columns=["species_tax_id", "genome_id"])]).drop_duplicates(
        "species_tax_id")
    out.to_csv(f, sep="\t", index=False)
    print("accessions found:", out.genome_id.notna().sum(), "/", len(out))


def step2():
    togo = importlib.import_module("20_togo_ko_profiles")
    acc = pd.read_csv(OUT / "species_accessions.tsv", sep="\t").dropna(subset=["genome_id"])
    ids = sorted(set(acc.genome_id))
    f = OUT / "togo_kofam.parquet"
    have = set(pd.read_parquet(f, columns=["genome_id"]).genome_id.unique()) if f.exists() else set()
    todo = [a for a in ids if a not in have]
    print(f"genomes: {len(ids)}; profiles cached: {len(have)}; to fetch: {len(todo)}", flush=True)
    rows, log = [], []
    with ThreadPoolExecutor(3) as ex:
        for k, (a, out, status) in enumerate(ex.map(togo.fetch, todo)):
            rows += out
            log.append((a, status, len(out)))
            if k % 250 == 0:
                print(f"  {k}/{len(todo)} fetched, {sum(s == 'ok' for _, s, _ in log)} with data", flush=True)
                if rows:
                    df = pd.DataFrame(rows, columns=["genome_id", "protein", "KO", "score", "threshold", "evalue",
                                                     "significant"])
                    (pd.concat([pd.read_parquet(f), df]) if f.exists() else df).to_parquet(f, index=False)
                    have |= set(df.genome_id)
                    rows = []
    if rows:
        df = pd.DataFrame(rows, columns=["genome_id", "protein", "KO", "score", "threshold", "evalue", "significant"])
        (pd.concat([pd.read_parquet(f), df]) if f.exists() else df).to_parquet(f, index=False)
    pd.DataFrame(log, columns=["genome_id", "status", "n"]).to_csv(OUT / "togo_fetch_log.tsv", sep="\t", index=False)
    hits = pd.read_parquet(f)
    sig = hits[hits.significant].drop_duplicates(["genome_id", "KO"])
    mat = pd.crosstab(sig.genome_id, sig.KO).clip(upper=1).astype("int8")
    mat.index.name = "genome_id"
    mat.to_csv(OUT / "ko_matrix.tsv", sep="\t")
    print("KO matrix:", mat.shape)


def step3():
    ls = importlib.import_module("05_labels_splits")
    d = madin()
    acc = pd.read_csv(OUT / "species_accessions.tsv", sep="\t").dropna(subset=["genome_id"])
    have = set(pd.read_csv(OUT / "ko_matrix.tsv", sep="\t", usecols=["genome_id"]).genome_id)
    d = d.merge(acc, on="species_tax_id")
    d = d[d.genome_id.isin(have)].drop_duplicates("genome_id").reset_index(drop=True)
    paths = d.pathways.fillna("").str.split(",").apply(lambda xs: [x.strip() for x in xs])
    for t, keys in TRAITS.items():
        d[t] = paths.apply(lambda ps: float(any(p == k or p.startswith(k + "_") for p in ps for k in keys)))
    traits = [t for t in TRAITS if d[t].sum() >= MIN_POS and (d[t] == 0).sum() >= MIN_POS]
    d["faprotax_source"] = d.data_source.fillna("").str.contains("faprotax")

    meta_cols = ["accession", "gtdb_taxonomy"]
    gm = pd.concat([pd.read_csv(ROOT / "data/raw/gtdb" / f, sep="\t", usecols=meta_cols, low_memory=False)
                    for f in ("ar53_metadata.tsv.gz", "bac120_metadata.tsv.gz")])
    gm["acc"] = gm.accession.str[3:]
    for rank, code in (("genus", "g"), ("family", "f")):
        gmap = dict(zip(gm.acc, gm.gtdb_taxonomy.str.extract(rf"{code}__([^;]*)")[0]))
        d[f"gtdb_{rank}"] = d.genome_id.map(gmap)
    # Split groups: GTDB rank where available, else the NCBI rank, else the genome itself. A genome with no
    # resolvable genus/family becomes its own singleton group, which is the conservative choice: it can never share
    # a group with another genome, so it is always fully held out. Groups must never be NaN - a missing group would
    # make the grouped splitter compare float to str and, worse, silently merge unrelated genomes.
    gen_ok = d.gtdb_genus.notna() & (d.gtdb_genus.astype(str) != "")
    d["group_genus"] = np.where(gen_ok, "g:" + d.gtdb_genus.astype(str),
                                np.where(d.genus.notna(), "ncbi:" + d.genus.astype(str),
                                         "self:" + d.genome_id.astype(str)))
    fam_ok = d.gtdb_family.notna() & (d.gtdb_family.astype(str) != "")
    d["group_family"] = np.where(fam_ok, "f:" + d.gtdb_family.astype(str),
                                 np.where(d.family.notna(),
                                          "ncbi:" + d.family.astype(str) + "|" + d.group_genus,
                                          "self:" + d.genome_id.astype(str)))
    assert not (pd.isna(d.group_genus).any() or pd.isna(d.group_family).any()), "NaN split group"
    Y = d[traits].to_numpy()
    d["fold_genus"] = ls.greedy_group_folds(d.group_genus.to_numpy(), Y, N_FOLDS, SEED)
    d["fold_family"] = ls.greedy_group_folds(d.group_family.to_numpy(), Y, N_FOLDS, SEED)
    keep = (["species_tax_id", "species", "genus", "family", "order", "class", "phylum", "superkingdom", "genome_id",
             "gtdb_genus", "gtdb_family", "faprotax_source", "data_source", "pathways", "group_genus", "group_family",
             "fold_genus", "fold_family"] + traits)
    d[keep].to_csv(OUT / "labels.tsv", sep="\t", index=False)
    digest = hashlib.sha256((OUT / "labels.tsv").read_bytes()).hexdigest()
    (OUT / "labels.sha256").write_text(f"{digest}  labels.tsv\n")
    counts = pd.DataFrame({"trait": traits, "n_pos": [int(d[t].sum()) for t in traits],
                           "n_neg": [int((d[t] == 0).sum()) for t in traits],
                           "prevalence": [float(d[t].mean()) for t in traits],
                           "n_genera_with_pos": [int(d.loc[d[t] == 1, "group_genus"].nunique()) for t in traits],
                           "frac_pos_faprotax": [float(d.loc[d[t] == 1, "faprotax_source"].mean()) for t in traits]})
    counts.to_csv(ROOT / "results/pathway_trait_counts.csv", index=False, float_format="%.3f")
    print(f"genomes: {len(d)}, evaluable traits: {len(traits)}, genus groups: {d.group_genus.nunique()}, "
          f"family groups: {d.group_family.nunique()}")
    print(counts.to_string(index=False))
    print("labels sha256", digest)


if __name__ == "__main__":
    step = sys.argv[sys.argv.index("--step") + 1] if "--step" in sys.argv else "all"
    mx = int(sys.argv[sys.argv.index("--max-species") + 1]) if "--max-species" in sys.argv else None
    if step in ("1", "all"):
        step1(mx)
    if step in ("2", "all"):
        step2()
    if step in ("3", "all"):
        step3()
