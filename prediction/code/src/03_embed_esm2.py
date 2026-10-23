"""Per-protein ESM-2 t12 35M embeddings (the Bacformer-base protein encoder), content-addressed.

Pooling replicates bacformer.pp.embed_prot_seqs.generate_protein_embeddings (HF path):
mean of last_hidden_state over attention_mask (includes <cls>/<eos>), truncation to
1024 tokens. Sequences are length-bucketed and padded to fixed bucket lengths to avoid
MPS kernel recompilation; padding is excluded by the attention mask.

Sequences from labelled genomes are embedded first so the evaluation can start early.
Output: data/processed/esm2_t12_35M/shard_XXXX.npz  (hashes: str[N], emb: float16[N,480])
Resumable: hashes present in completed shards are skipped.
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from transformers import AutoModel, AutoTokenizer

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/processed/esm2_t12_35M"
MODEL = "facebook/esm2_t12_35M_UR50D"
REVISION = "6fbf070e65b0b7291e7bbcd451118c216cff79d8"
MAX_LEN = 1024
BUCKETS = [64, 128, 192, 256, 320, 384, 512, 640, 768, 1024]
TOKEN_BUDGET = 16384
SHARD_SIZE = 50_000


def done_hashes():
    s = set()
    for f in sorted(OUT.glob("shard_*.npz")):
        s.update(np.load(f, allow_pickle=False)["hashes"].tolist())
    return s


def load_model(dev):
    tok = AutoTokenizer.from_pretrained(MODEL, revision=REVISION)
    model = AutoModel.from_pretrained(MODEL, revision=REVISION, add_pooling_layer=False).eval().to(dev)
    return tok, (model.half() if dev != "cpu" else model)


def embed_batches(seqs, tok, model, dev):
    """Yield (start, end, float16 embeddings) over length-sorted `seqs` in fixed-bucket batches."""
    i = 0
    while i < len(seqs):
        L = min(len(seqs[i]) + 2, MAX_LEN)
        B = next(b for b in BUCKETS if b >= L)
        bs = max(1, TOKEN_BUDGET // B)
        j = i  # keep the batch inside one bucket
        while j < len(seqs) and j - i < bs and min(len(seqs[j]) + 2, MAX_LEN) <= B:
            j += 1
        inp = tok(seqs[i:j], return_tensors="pt", padding="max_length", max_length=B, truncation=True)
        inp = {k: v.to(dev) for k, v in inp.items()}
        with torch.no_grad():
            h = model(**inp).last_hidden_state.float()
        m = inp["attention_mask"].float()
        emb = (torch.einsum("ijk,ij->ik", h, m) / m.sum(1, keepdim=True)).cpu().numpy().astype(np.float16)
        if not np.isfinite(emb).all():
            sys.exit(f"non-finite embedding in batch starting at {i}")
        yield i, j, emb
        i = j


def embed_sequences(sequences, dev):
    """Embed arbitrary sequences (any order); returns float16 array aligned to input order."""
    tok, model = load_model(dev)
    order = np.argsort([len(s) for s in sequences], kind="stable")
    sorted_seqs = [sequences[k] for k in order]
    out = np.zeros((len(sequences), model.config.hidden_size), dtype=np.float16)
    for i, j, emb in embed_batches(sorted_seqs, tok, model, dev):
        out[order[i:j]] = emb
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    dev = "mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu")
    tok, model = load_model(dev)

    labels = pd.read_csv(ROOT / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig")
    unlabeled = labels["Donor categories: Not in the electron donor / acceptor database"] == 1
    prio = set()
    for acc in labels.loc[~unlabeled, "genome_id"]:
        prio.update(pd.read_parquet(ROOT / f"data/processed/proteins/{acc}.parquet", columns=["seq_hash"]).seq_hash)
    u = pd.read_parquet(ROOT / "data/processed/unique_sequences.parquet")
    done = done_hashes()
    u = u[~u.seq_hash.isin(done)].copy()
    u["prio"] = ~u.seq_hash.isin(prio)  # False (labelled) sorts first
    u = u.sort_values(["prio", "length"]).reset_index(drop=True)
    print(f"device={dev} done={len(done)} todo={len(u)} (labelled-first {int((~u.prio).sum())})", flush=True)

    shard_id = len(list(OUT.glob("shard_*.npz")))
    buf_h, buf_e = [], []
    seqs, hashes = u.sequence.tolist(), u.seq_hash.tolist()
    t0, n = time.time(), 0
    for i, j, emb in embed_batches(seqs, tok, model, dev):
        buf_h += hashes[i:j]
        buf_e.append(emb)
        n += j - i
        if len(buf_h) >= SHARD_SIZE or j >= len(seqs):
            np.savez(OUT / f"shard_{shard_id:04d}.npz", hashes=np.array(buf_h), emb=np.concatenate(buf_e))
            shard_id += 1
            buf_h, buf_e = [], []
            rate = n / (time.time() - t0)
            print(f"shard {shard_id} saved; {n}/{len(seqs)} ; {rate:.1f} prot/s ; eta {(len(seqs) - n) / rate / 3600:.2f} h", flush=True)


if __name__ == "__main__":
    main()
