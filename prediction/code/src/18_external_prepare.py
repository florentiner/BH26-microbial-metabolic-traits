"""External validation set (docs/external_validation.md): labels, genome selection, download, protein tables.

Output: data/external/labels.tsv (species, genome_id, traits with NaN = not evaluable),
        data/external/genomes/<acc>/..., data/external/proteins/<acc>.parquet,
        data/external/unique_sequences_new.parquet (sequences not already in the internal set)
"""
import importlib
import json
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data/external"
SEED = 20260917
LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
MAP = {
    "acceptor_sulfate": ["sulfate_reduction"],
    "acceptor_sulfur": ["sulfur_reduction", "thiosulfate_reduction"],
    "acceptor_CO2_methanogen_acetogen": ["methanogenesis", "acetogenesis"],
    "acceptor_iron": ["iron_reduction"],
    "acceptor_nitrate": ["denitrification", "nitrate_reduction", "nitrate_reduction_to_ammonia"],
    "donor_hydrogen_oxidation": ["hydrogen_oxidation"],
    "donor_sulfur_oxidation": ["sulfur_oxidation", "thiosulfate_oxidation", "sulfide_oxidation", "sulfite_oxidation",
                               "anoxygenic_photoautotrophy_sulfur_oxidation"],
    "donor_iron_oxidation": ["iron_oxidation"],
    "nitrification": ["nitrification"],
}
AEROBIC = {"aerobic", "obligate aerobic", "facultative", "microaerophilic"}
ANAEROBIC = {"anaerobic", "obligate anaerobic"}


def label_species():
    d = pd.read_csv(ROOT / "data/raw/madin/condensed_species_NCBI.csv", low_memory=False)
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    d = d[d.pathways.notna()].copy()
    d = d[~d.species.str.lower().isin(set(lab.species.str.lower())) & ~d.genus.str.lower().isin(set(lab.genus.str.lower()))]
    paths = d.pathways.str.split(",").apply(lambda xs: [x.strip() for x in xs])
    for t, keys in MAP.items():  # prefix match: methanogenesis_H2CO2, hydrogen_oxidation_dark, ...
        d[t] = paths.apply(lambda ps: float(any(p == k or p.startswith(k + "_") for p in ps for k in keys)))
    met = d.metabolism.fillna("")
    d["acceptor_oxygen"] = np.where(met.isin(AEROBIC), 1.0, np.where(met.isin(ANAEROBIC), 0.0, np.nan))
    d["faprotax_source"] = d.data_source.fillna("").str.contains("faprotax")
    return d


def ncbi_best_genome(taxid):
    url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{int(taxid)}/dataset_report"
           "?filters.assembly_source=refseq&filters.has_annotation=true&page_size=100")
    for i in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                reps = json.load(r).get("reports", [])
            break
        except Exception:
            time.sleep(5 * (i + 1))
    else:
        return None
    best = None
    for r in reps:
        a = r.get("assembly_info", {})
        key = (0 if a.get("refseq_category") == "reference genome" else 1,
               LEVEL_RANK.get(a.get("assembly_level"), 9),
               -int((a.get("release_date") or "0000-00-00").replace("-", "")[:8] or 0))  # most recent first
        if best is None or key < best[0]:
            best = (key, r["accession"])
    time.sleep(0.25)
    return best[1] if best else None


def main():
    EXT.mkdir(parents=True, exist_ok=True)
    d = label_species()
    trait_cols = list(MAP) + ["acceptor_oxygen"]
    pos_any = d[[t for t in MAP if t != "acceptor_nitrate"]].max(axis=1) == 1
    pos = d[pos_any]
    neg_pool = d[~pos_any].sample(frac=1.0, random_state=SEED)
    print("novel-genus pathway species:", len(d), "| positives for a mapped non-nitrate pathway:", len(pos))
    with ThreadPoolExecutor(3) as ex:
        pos_acc = list(ex.map(ncbi_best_genome, pos.species_tax_id))
    pos = pos.assign(genome_id=pos_acc).dropna(subset=["genome_id"])
    negs, k = [], 0
    while len(negs) < len(pos) and k < len(neg_pool):
        chunk = neg_pool.iloc[k:k + 60]
        k += 60
        with ThreadPoolExecutor(3) as ex:
            accs = list(ex.map(ncbi_best_genome, chunk.species_tax_id))
        negs.append(chunk.assign(genome_id=accs).dropna(subset=["genome_id"]))
        negs = [pd.concat(negs)]
    neg = negs[0].iloc[:len(pos)]
    sel = pd.concat([pos.assign(selection="positive_enriched"), neg.assign(selection="background")])
    sel = sel.drop_duplicates("genome_id")
    out = sel[["species_tax_id", "species", "genus", "family", "phylum", "superkingdom", "genome_id", "selection",
               "faprotax_source", "data_source", "pathways", "metabolism"] + trait_cols]
    out.to_csv(EXT / "labels.tsv", sep="\t", index=False)
    print("selected genomes:", len(out), out.selection.value_counts().to_dict())
    print(out[trait_cols].sum().to_string())

    dl = importlib.import_module("01_download_genomes")
    pt = importlib.import_module("02_build_protein_tables")
    dl.OUT = pt.RAW = EXT / "genomes"
    pt.OUT = EXT / "proteins"
    pt.OUT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(4) as ex:
        status = list(ex.map(dl.download_one, out.genome_id))
    st = pd.DataFrame(status, columns=["genome_id", "status", "n_proteins", "faa_sha256"])
    st.to_csv(EXT / "download_manifest.tsv", sep="\t", index=False)
    print(st.status.value_counts().to_dict())
    internal = set(pd.read_parquet(ROOT / "data/processed/unique_sequences.parquet", columns=["seq_hash"]).seq_hash)
    uniq, qc = {}, []
    for acc in st.loc[st.status == "ok", "genome_id"]:
        qc.append(pt.build(acc, uniq))
    pd.DataFrame(qc).to_csv(EXT / "protein_table_qc.tsv", sep="\t", index=False)
    new = pd.DataFrame({"seq_hash": [h for h in uniq if h not in internal]})
    new["sequence"] = new.seq_hash.map(uniq)
    new["length"] = new.sequence.str.len()
    new.to_parquet(EXT / "unique_sequences_new.parquet", index=False)
    print("external unique sequences:", len(uniq), "new (not in internal set):", len(new))


if __name__ == "__main__":
    main()
