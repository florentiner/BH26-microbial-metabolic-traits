"""Build genome-coordinate-ordered protein tables (protocol gate G4).

RefSeq protein.faa is non-redundant (one record per WP_ accession) and sorted by
accession, so its record order is NOT genomic order and identical paralogs are
collapsed. We therefore rebuild each proteome from the GFF3: one protein per CDS
locus (Parent gene), ordered by (contig order in GFF, start coordinate).

Output:
  data/processed/proteins/<acc>.parquet  columns: pos, contig, contig_idx, start, end,
                                          strand, locus, protein_id, seq_hash, length
  data/processed/unique_sequences.parquet  seq_hash, sequence, length  (content-addressed)
  data/processed/protein_table_qc.tsv      per-genome QC counts
"""
import hashlib
from pathlib import Path
from urllib.parse import unquote

import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/genomes"
OUT = ROOT / "data/processed/proteins"


def parse_attrs(s):
    return dict(kv.split("=", 1) for kv in s.strip().split(";") if "=" in kv)


def build(acc, uniq):
    seqs = {r.id: str(r.seq).rstrip("*") for r in SeqIO.parse(RAW / acc / "protein.faa", "fasta")}
    loci, contig_order = {}, {}
    n_cds_rows = n_no_pid = 0
    with open(RAW / acc / "genomic.gff") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "CDS":
                continue
            n_cds_rows += 1
            a = parse_attrs(f[8])
            pid = a.get("protein_id")
            if pid is None:  # pseudogene CDS without translation
                n_no_pid += 1
                continue
            contig_order.setdefault(f[0], len(contig_order))
            key = (f[0], a.get("Parent", a.get("ID")))  # multi-segment CDS share one Parent
            start, end = int(f[3]), int(f[4])
            if key in loci:
                loci[key]["start"] = min(loci[key]["start"], start)
                loci[key]["end"] = max(loci[key]["end"], end)
            else:
                loci[key] = dict(contig=f[0], start=start, end=end, strand=f[6],
                                 locus=unquote(key[1]), protein_id=pid)
    df = pd.DataFrame(loci.values())
    missing = sorted(set(df.protein_id) - set(seqs))
    df = df[df.protein_id.isin(seqs)].copy()
    df["contig_idx"] = df.contig.map(contig_order)
    df = df.sort_values(["contig_idx", "start", "end"]).reset_index(drop=True)
    df["pos"] = range(len(df))
    hashes, lengths = [], []
    for pid in df.protein_id:
        s = seqs[pid]
        h = hashlib.md5(s.encode()).hexdigest()
        uniq.setdefault(h, s)
        hashes.append(h)
        lengths.append(len(s))
    df["seq_hash"], df["length"] = hashes, lengths
    df = df[["pos", "contig", "contig_idx", "start", "end", "strand", "locus", "protein_id", "seq_hash", "length"]]
    df.to_parquet(OUT / f"{acc}.parquet", index=False)
    return dict(genome_id=acc, n_faa_records=len(seqs), n_cds_rows=n_cds_rows, n_cds_without_protein=n_no_pid,
                n_loci=len(df), n_missing_protein_ids=len(missing), n_contigs=len(contig_order),
                n_faa_unused=len(set(seqs) - set(df.protein_id)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ids = pd.read_csv(ROOT / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig")["genome_id"]
    uniq, qc = {}, []
    for i, acc in enumerate(ids):
        qc.append(build(acc, uniq))
        if i % 100 == 0:
            print(i, qc[-1], flush=True)
    qc = pd.DataFrame(qc)
    qc.to_csv(ROOT / "data/processed/protein_table_qc.tsv", sep="\t", index=False)
    u = pd.DataFrame({"seq_hash": list(uniq), "sequence": list(uniq.values())})
    u["length"] = u.sequence.str.len()
    u.to_parquet(ROOT / "data/processed/unique_sequences.parquet", index=False)
    print(qc.describe().T)
    print("total loci", qc.n_loci.sum(), "unique sequences", len(u), "residues", u.length.sum())
    assert (qc.n_missing_protein_ids == 0).all(), "GFF protein_id without FASTA record"
    assert (qc.n_faa_unused == 0).all(), "FASTA record not placed on a GFF locus"


if __name__ == "__main__":
    main()
