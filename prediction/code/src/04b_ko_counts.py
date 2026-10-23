"""KO copy number per genome for the energy panel (research loop R1).

count = number of protein loci in the genome assigned the KO (KofamScan rule, data/processed/kofam_hits.parquet).
Output: data/processed/ko_counts_panel.tsv (genome x KO integer counts)
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data/processed"


def main():
    hits = pd.read_parquet(PROC / "kofam_hits.parquet", columns=["seq_hash", "KO"]).drop_duplicates()
    panel = pd.read_csv(PROC / "ko_panel.tsv", sep="\t")
    ids = pd.read_csv(ROOT / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig").genome_id
    rows = {}
    for acc in ids:
        prot = pd.read_parquet(PROC / "proteins" / f"{acc}.parquet", columns=["seq_hash"])
        rows[acc] = prot.merge(hits, on="seq_hash").KO.value_counts()
    mat = pd.DataFrame(rows).T.reindex(columns=sorted(panel.KO)).fillna(0).astype(int)
    mat.index.name = "genome_id"
    pres = pd.read_csv(PROC / "ko_matrix_panel.tsv", sep="\t", index_col=0)
    assert ((mat > 0).astype(int).loc[pres.index, pres.columns].values == pres.values).all(), "counts != presence"
    mat.to_csv(PROC / "ko_counts_panel.tsv", sep="\t")
    print(mat.shape, "max copy number", int(mat.values.max()), "cells >1:", int((mat.values > 1).sum()))


if __name__ == "__main__":
    main()
