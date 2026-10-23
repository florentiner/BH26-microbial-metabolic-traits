"""Current NCBI lineage + GTDB genus for external genomes, and a lineage-based novelty filter (audit round 2, M3/M4).

Madin et al. (2020) species/genus names predate NCBI renames (e.g. Desulfotomaculum -> Desulfofundulus), so name
matching let training genera and even training genomes into the "novel genus" set. Here every external genome gets
its current NCBI lineage (NCBI Datasets taxonomy API, from the assembly taxid) and its GTDB R232 genus. A genome is
novel only if (a) its accession is not among the 586 internal genomes, (b) its current NCBI genus is not a training
genus, and (c) its GTDB genus (when available) is not the GTDB genus of any labelled training genome.
Output: data/external/lineage.tsv (genome_id, taxid, ncbi_domain..ncbi_species, gtdb_genus, novel_*)
"""
import json
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data/external"
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def taxid_of(acc):
    rep = json.loads((EXT / "genomes" / acc / "assembly_report.json").read_text().splitlines()[0])
    return rep.get("organism", {}).get("taxId")


def lineage(taxid):
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{taxid}/dataset_report"
    for i in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                c = json.load(r)["reports"][0]["taxonomy"].get("classification", {})
            time.sleep(0.25)
            return {f"ncbi_{k}": c.get(k, {}).get("name") for k in RANKS}
        except Exception:
            time.sleep(5 * (i + 1))
    return {f"ncbi_{k}": None for k in RANKS}


def main():
    st = pd.read_csv(EXT / "download_manifest.tsv", sep="\t")
    ids = st.loc[st.status == "ok", "genome_id"].tolist()
    taxids = [taxid_of(a) for a in ids]
    with ThreadPoolExecutor(3) as ex:
        lin = list(ex.map(lineage, taxids))
    df = pd.DataFrame(lin)
    df.insert(0, "taxid", taxids)
    df.insert(0, "genome_id", ids)

    cols = ["accession", "gtdb_taxonomy"]
    meta = pd.concat([pd.read_csv(ROOT / "data/raw/gtdb" / f, sep="\t", usecols=cols, low_memory=False)
                      for f in ("ar53_metadata.tsv.gz", "bac120_metadata.tsv.gz")])
    meta["acc"] = meta.accession.str[3:]
    meta["gtdb_genus"] = meta.gtdb_taxonomy.str.extract(r"g__([^;]*)")[0]
    gg = dict(zip(meta.acc, meta.gtdb_genus))
    df["gtdb_genus"] = df.genome_id.map(gg)

    internal_ids = set(pd.read_csv(ROOT / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig").genome_id)
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    train_genera = set(lab.loc[lab.labelled, "genus"].str.lower())
    gt = pd.read_csv(ROOT / "data/processed/gtdb_taxonomy.tsv", sep="\t")
    train_gtdb_genera = set(gt.loc[gt.genome_id.isin(lab.loc[lab.labelled, "genome_id"]), "gtdb_genus"].dropna())
    df["in_internal_586"] = df.genome_id.isin(internal_ids)
    df["ncbi_genus_in_training"] = df.ncbi_genus.fillna("").str.lower().isin(train_genera)
    df["gtdb_genus_in_training"] = df.gtdb_genus.isin(train_gtdb_genera)
    df["novel"] = ~(df.in_internal_586 | df.ncbi_genus_in_training | df.gtdb_genus_in_training)
    df.to_csv(EXT / "lineage.tsv", sep="\t", index=False)
    print("lineage resolved:", int(df.ncbi_genus.notna().sum()), "/", len(df))
    print("in internal 586:", int(df.in_internal_586.sum()), "| NCBI genus in training:", int(df.ncbi_genus_in_training.sum()),
          "| GTDB genus in training:", int(df.gtdb_genus_in_training.sum()), "| novel:", int(df.novel.sum()))


if __name__ == "__main__":
    main()
