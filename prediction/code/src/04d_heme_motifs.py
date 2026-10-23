"""Multiheme c-type cytochrome features (research loop R11): CXXCH heme-binding motif counts per protein.

Per genome: number of proteins with >=2, >=3, >=5, >=10 CXXCH motifs (log1p), their fraction of the proteome,
total motif count per 1000 proteins, and the maximum motifs in a single protein. Label-free, sequence-only.
Output: data/processed/genome_features_heme.parquet (and data/external/ for external genomes if present)
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MOTIF = re.compile(r"(?=C..CH)")


def genome_features(prot, motif_count):
    m = prot.seq_hash.map(motif_count).fillna(0).to_numpy()
    n = len(m)
    f = {}
    for k in (2, 3, 5, 10):
        f[f"heme_n_ge{k}"] = np.log1p((m >= k).sum())
        f[f"heme_frac_ge{k}"] = (m >= k).mean() * 1000
    f["heme_motifs_per_1000_prot"] = m.sum() / n * 1000
    f["heme_max_in_protein"] = np.log1p(m.max())
    return f


def motif_counts(seq_file):
    u = pd.read_parquet(seq_file, columns=["seq_hash", "sequence"])
    counts = dict(zip(u.seq_hash, (len(MOTIF.findall(s)) for s in u.sequence)))
    return counts


def main():
    proc = ROOT / "data/processed"
    counts = motif_counts(proc / "unique_sequences.parquet")
    ext = ROOT / "data/external"
    if (ext / "unique_sequences_new.parquet").exists():
        counts.update(motif_counts(ext / "unique_sequences_new.parquet"))
    targets = [(proc / "proteins", pd.read_csv(ROOT / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig")
                .genome_id, proc / "genome_features_heme.parquet")]
    if "--external" in sys.argv:
        st = pd.read_csv(ext / "download_manifest.tsv", sep="\t")
        targets = [(ext / "proteins", st.loc[st.status == "ok", "genome_id"], ext / "genome_features_heme.parquet")]
    for pdir, ids, out in targets:
        rows = {a: genome_features(pd.read_parquet(pdir / f"{a}.parquet", columns=["seq_hash"]), counts) for a in ids}
        df = pd.DataFrame(rows).T.astype(np.float32)
        df.index.name = "genome_id"
        df.to_parquet(out)
        print(out.name, df.shape, df.describe().loc[["mean", "max"]].round(2).to_string())


if __name__ == "__main__":
    main()
