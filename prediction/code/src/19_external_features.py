"""Annotate external genomes with exactly the internal feature pipelines (KOfam energy panel, MagicLamp LithoGenie).

Sequences already scanned internally reuse their hits (identical sequence -> identical scores); only new
sequences are searched, with the same HMMs and thresholds.
Output: data/external/ko_matrix_panel.tsv, data/external/genome_features_lithogenie.parquet
"""
import importlib

import numpy as np
import pandas as pd

from redox_common import PROC, ROOT

EXT = ROOT / "data/external"
ko = importlib.import_module("04_kofam_scan")
ml = importlib.import_module("04c_magiclamp_scan")


def main():
    st = pd.read_csv(EXT / "download_manifest.tsv", sep="\t")
    ids = st.loc[st.status == "ok", "genome_id"].tolist()
    new = pd.read_parquet(EXT / "unique_sequences_new.parquet")
    print("external genomes", len(ids), "new sequences to scan", len(new), flush=True)
    genome_hashes = {a: pd.read_parquet(EXT / "proteins" / f"{a}.parquet", columns=["seq_hash"]).seq_hash for a in ids}

    _, panel = ko.build_panel()
    hmms, thr, stype = ko.load_hmms(panel)
    ko_new = ko.scan(hmms, thr, stype, new.seq_hash.tolist(), new.sequence.tolist())
    ko_hits = pd.concat([pd.read_parquet(PROC / "kofam_hits.parquet"), ko_new], ignore_index=True)
    mat = ko.presence_matrix(ko_hits, genome_hashes, panel.KO)
    mat.to_csv(EXT / "ko_matrix_panel.tsv", sep="\t")
    print("KO matrix", mat.shape, flush=True)

    lib = ml.load_library()
    import pyhmmer
    hm = []
    for name, p in zip(lib.hmm, lib.path):
        with pyhmmer.plan7.HMMFile(p) as fh:
            h = fh.read()
        h.name = name.encode()
        hm.append(h)
    cutoff = dict(zip(lib.hmm, lib.bitscore_cutoff))
    seqs = pyhmmer.easel.TextSequenceBlock(
        pyhmmer.easel.TextSequence(name=h.encode(), sequence=s) for h, s in zip(new.seq_hash, new.sequence)
    ).digitize(pyhmmer.easel.Alphabet.amino())
    rows = []
    for hits in pyhmmer.hmmer.hmmsearch(hm, seqs, cpus=8, E=10.0):
        q = hits.query.name.decode() if isinstance(hits.query.name, bytes) else hits.query.name
        for h in hits:
            if (h.score >= cutoff[q]) if cutoff[q] > 0 else (h.evalue <= ml.FALLBACK_E):
                rows.append((h.name.decode() if isinstance(h.name, bytes) else h.name, q, float(h.score)))
    lg_hits = pd.concat([pd.read_parquet(PROC / "magiclamp_hits.parquet"),
                         pd.DataFrame(rows, columns=["seq_hash", "hmm", "score"])], ignore_index=True)
    hh = lg_hits[["seq_hash", "hmm"]].drop_duplicates()
    counts = {a: pd.DataFrame({"seq_hash": hs}).merge(hh, on="seq_hash").hmm.value_counts()
              for a, hs in genome_hashes.items()}
    lgm = pd.DataFrame(counts).T.reindex(index=ids, columns=lib.hmm).fillna(0)
    feat = np.log1p(lgm).astype(np.float32)
    feat.columns = [f"lg_{c}" for c in feat.columns]
    feat.index.name = "genome_id"
    feat.to_parquet(EXT / "genome_features_lithogenie.parquet")
    print("LithoGenie features", feat.shape)


if __name__ == "__main__":
    main()
