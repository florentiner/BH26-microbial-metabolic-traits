"""Curated lithotrophy / iron redox HMMs (research loop R3): MagicLamp LithoGenie + FeGenie.

HMM sets (MagicLamp commit in data/raw/magiclamp/COMMIT):
  hmms/litho/*.hmm                         (LithoGenie; cutoffs hmms/litho/hmm-meta.txt)
  hmms/iron/iron_oxidation, iron_reduction  (FeGenie; cutoffs hmms/iron/HMM-bitcutoffs.txt)
Duplicated HMM names are kept once (LithoGenie copy). Assignment follows the tools: full-sequence bit score >=
the HMM's cutoff (hmmsearch -T). HMMs with no or zero cutoff use a conservative E-value <= 1e-10 instead of
accepting every reported hit (LithoGenie's default). All proteins are searched at once (bit scores do not depend
on database size; E-values do, which only affects the fallback HMMs).

Output: data/processed/magiclamp_hits.parquet, data/processed/genome_features_lithogenie.parquet
        (log1p count of loci per HMM per genome), data/processed/magiclamp_hmms.tsv
"""
import importlib
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyhmmer

ROOT = Path(__file__).resolve().parents[1]
ML = ROOT / "data/raw/magiclamp/MagicLamp-master/hmms"
PROC = ROOT / "data/processed"
FALLBACK_E = 1e-10


def load_library():
    meta = pd.read_csv(ML / "litho/hmm-meta.txt", sep="\t")
    cut = dict(zip(meta.hmm, pd.to_numeric(meta.bitscore_cutoff, errors="coerce")))
    react = dict(zip(meta.hmm, meta.reaction))
    fe = pd.read_csv(ML / "iron/HMM-bitcutoffs.txt", sep="\t")
    fe_cut = dict(zip(fe.HMM, pd.to_numeric(fe.soft_bitscore_cutoff, errors="coerce")))
    rows = []
    for f in sorted((ML / "litho").glob("*.hmm")):
        rows.append((f.stem, "LithoGenie", react.get(f.stem, "unknown"), cut.get(f.stem, np.nan), f))
    for sub in ("iron_oxidation", "iron_reduction"):
        for f in sorted((ML / "iron" / sub).glob("*.hmm")):
            rows.append((f.stem, "FeGenie", sub.replace("_", "-"), fe_cut.get(f.stem, np.nan), f))
    lib = pd.DataFrame(rows, columns=["hmm", "source", "reaction", "bitscore_cutoff", "path"])
    lib = lib.drop_duplicates("hmm", keep="first").reset_index(drop=True)
    lib["rule"] = np.where(lib.bitscore_cutoff > 0, "bitscore>=cutoff", f"evalue<={FALLBACK_E}")
    return lib


def main():
    lib = load_library()
    lib.drop(columns="path").to_csv(PROC / "magiclamp_hmms.tsv", sep="\t", index=False)
    print(lib.source.value_counts().to_dict(), "fallback-E HMMs:", int((lib.rule != "bitscore>=cutoff").sum()))
    hmms = []
    for name, p in zip(lib.hmm, lib.path):
        with pyhmmer.plan7.HMMFile(p) as fh:
            h = fh.read()
        h.name = name.encode() if isinstance(h.name, bytes) else name  # file NAME fields are not unique
        hmms.append(h)
    cutoff = dict(zip(lib.hmm, lib.bitscore_cutoff))
    u = pd.read_parquet(PROC / "unique_sequences.parquet", columns=["seq_hash", "sequence"])
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = pyhmmer.easel.TextSequenceBlock(
        pyhmmer.easel.TextSequence(name=h.encode(), sequence=s) for h, s in zip(u.seq_hash, u.sequence)
    ).digitize(alphabet)
    del u
    rows, t0 = [], time.time()
    for k, hits in enumerate(pyhmmer.hmmer.hmmsearch(hmms, seqs, cpus=8, E=10.0)):
        q = hits.query.name
        q = q.decode() if isinstance(q, bytes) else q
        c = cutoff[q]
        for h in hits:
            ok = h.score >= c if c > 0 else h.evalue <= FALLBACK_E
            if ok:
                n = h.name.decode() if isinstance(h.name, bytes) else h.name
                rows.append((n, q, float(h.score)))
        if k % 25 == 0:
            print(f"{k + 1}/{len(hmms)} {q} {time.time() - t0:.0f}s hits={len(rows)}", flush=True)
    hits = pd.DataFrame(rows, columns=["seq_hash", "hmm", "score"])
    hits.to_parquet(PROC / "magiclamp_hits.parquet", index=False)

    ids = pd.read_csv(ROOT / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig").genome_id
    hh = hits[["seq_hash", "hmm"]].drop_duplicates()
    counts = {}
    for acc in ids:
        prot = pd.read_parquet(PROC / "proteins" / f"{acc}.parquet", columns=["seq_hash"])
        counts[acc] = prot.merge(hh, on="seq_hash").hmm.value_counts()
    mat = pd.DataFrame(counts).T.reindex(columns=lib.hmm).fillna(0)
    mat.index.name = "genome_id"
    feat = np.log1p(mat).astype(np.float32)
    feat.columns = [f"lg_{c}" for c in feat.columns]
    feat.to_parquet(PROC / "genome_features_lithogenie.parquet")
    print("done", feat.shape, "HMMs present in >=1 genome:", int((mat.sum() > 0).sum()), f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    importlib.invalidate_caches()
    main()
