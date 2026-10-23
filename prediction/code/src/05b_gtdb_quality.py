"""GTDB R232 taxonomy and CheckM2 genome quality for every genome (sensitivity analyses S2 and S4).

Matching: GTDB accession RS_<GCF> or GB_<paired GCA>, version-insensitive fallback on the accession stem.
Output: data/processed/gtdb_taxonomy.tsv, data/processed/genome_qc.tsv
  quality_ok = CheckM2 completeness >= 90 and contamination <= 5
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/gtdb"
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def main():
    cols = ["accession", "gtdb_taxonomy", "checkm2_completeness", "checkm2_contamination",
            "ncbi_genbank_assembly_accession"]
    meta = pd.concat([pd.read_csv(RAW / f, sep="\t", usecols=cols, low_memory=False)
                      for f in ("ar53_metadata.tsv.gz", "bac120_metadata.tsv.gz")], ignore_index=True)
    meta["stem"] = meta.accession.str[3:].str.split(".").str[0].str[4:]  # numeric part shared by GCF/GCA
    meta["acc_nover"] = meta.accession.str[3:].str.split(".").str[0]
    by_acc = meta.set_index(meta.accession.str[3:])
    by_nover = meta.drop_duplicates("acc_nover").set_index("acc_nover")
    by_gca = meta.drop_duplicates("ncbi_genbank_assembly_accession").set_index("ncbi_genbank_assembly_accession")

    tsv = pd.read_csv(ROOT / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig")
    rows = []
    for acc in tsv.genome_id:
        paired = None
        rep = ROOT / "data/raw/genomes" / acc / "assembly_report.json"
        if rep.exists():
            r = json.loads(rep.read_text().splitlines()[0])
            paired = r.get("pairedAccession") or r.get("paired_accession")
        hit, how = None, None
        for key, table, label in ((acc, by_acc, "exact"), (acc.split(".")[0], by_nover, "version"),
                                  (paired, by_gca, "paired_gca"),
                                  (paired.split(".")[0] if paired else None, by_nover, "paired_gca_version")):
            if key is not None and key in table.index:
                hit, how = table.loc[key], label
                break
        row = {"genome_id": acc, "gtdb_match": how}
        if hit is not None:
            tax = dict(x.split("__", 1) for x in hit.gtdb_taxonomy.split(";"))
            for rank, code in zip(RANKS, "d p c o f g s".split()):
                row[f"gtdb_{rank}"] = tax.get(code) or None
            row["checkm2_completeness"] = hit.checkm2_completeness
            row["checkm2_contamination"] = hit.checkm2_contamination
        rows.append(row)
    g = pd.DataFrame(rows)
    g.to_csv(ROOT / "data/processed/gtdb_taxonomy.tsv", sep="\t", index=False)
    qc = g[["genome_id", "checkm2_completeness", "checkm2_contamination"]].copy()
    qc["quality_ok"] = (qc.checkm2_completeness >= 90) & (qc.checkm2_contamination <= 5)
    level = {}
    for acc in tsv.genome_id:
        rep = ROOT / "data/raw/genomes" / acc / "assembly_report.json"
        r = json.loads(rep.read_text().splitlines()[0]) if rep.exists() else {}
        level[acc] = r.get("assemblyInfo", r.get("assembly_info", {})).get("assemblyLevel",
                     r.get("assembly_info", {}).get("assembly_level"))
    qc["assembly_level"] = qc.genome_id.map(level)
    qc.to_csv(ROOT / "data/processed/genome_qc.tsv", sep="\t", index=False)
    print("GTDB match:", g.gtdb_match.value_counts(dropna=False).to_dict())
    print("quality_ok:", qc.quality_ok.value_counts().to_dict(), "assembly levels:",
          qc.assembly_level.value_counts(dropna=False).to_dict())
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    m = g.merge(lab[["genome_id", "genus", "family", "labelled"]], on="genome_id")
    ml = m[m.labelled]
    print("labelled genomes without GTDB genus:", ml.gtdb_genus.isna().sum())
    print("labelled: NCBI genera", ml.genus.nunique(), "GTDB genera", ml.gtdb_genus.nunique(),
          "| NCBI families", ml.family.nunique(), "GTDB families", ml.gtdb_family.nunique())
    print("NCBI genera split across >1 GTDB genus:",
          int((ml.groupby("genus").gtdb_genus.nunique() > 1).sum()))


if __name__ == "__main__":
    main()
