"""Bacformer-large genome vectors computed locally on the Apple GPU (MPS), for genomes the Kaggle quota cannot cover.

Same computation as kaggle/esmc_bacformer_large/run_template.py: ESM-C 300M (Synthyra/ESMplusplus_small @ f840ce88)
mask-weighted mean per protein, then Bacformer-large (macwiatrak/bacformer-large-masked-complete-genomes @ 84f85c37)
over windows of 6 000 proteins in genome order with contig ids clipped relative to the window start, mean over
positions. Equivalence to the GPU pipeline was checked earlier on a shared genome (cosine 0.9999985).

Proteins come from the local predict cache, so this runs offline once the proteomes are downloaded.

usage: python src/74_local_embed.py <accessions.txt> <out.parquet> [--device mps|cpu]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
PLM, PLM_REV = "Synthyra/ESMplusplus_small", "f840ce8824a3f8cbd025e6a0d0462c3d70fa4ac7"
BFL, BFL_REV = "macwiatrak/bacformer-large-masked-complete-genomes", "84f85c37a05e55559142b842c1569e265183b554"
BUCKETS = [128, 256, 384, 512, 768, 1024]
TOKEN_BUDGET = 16384


def proteins(acc):
    """(contig_idx, seq_hash, sequence) per protein in genome order, from the predict cache."""
    sys.path.insert(0, str(ROOT / "src"))
    import importlib
    pt = importlib.import_module("02_build_protein_tables")
    pt.RAW = ROOT / "data/predict_cache/genomes"
    pt.OUT = ROOT / "data/predict_cache/proteins"
    uniq = {}
    pt.build(acc, uniq)
    df = pd.read_parquet(pt.OUT / f"{acc}.parquet")
    return [(int(r.contig_idx), r.seq_hash) for r in df.itertuples()], uniq


def esmc_embed(seqs, device):
    from transformers import AutoModel
    m = AutoModel.from_pretrained(PLM, revision=PLM_REV, trust_remote_code=True).eval().to(device)
    tok = m.tokenizer
    items = sorted(seqs.items(), key=lambda kv: len(kv[1]))
    out, i = {}, 0
    while i < len(items):
        L = min(len(items[i][1]) + 2, 1024)
        B = next(b for b in BUCKETS if b >= L)
        j = i
        while j < len(items) and j - i < max(1, TOKEN_BUDGET // B) and min(len(items[j][1]) + 2, 1024) <= B:
            j += 1
        inp = tok([s for _, s in items[i:j]], return_tensors="pt", padding="max_length", max_length=B,
                  truncation=True)
        inp = {k: v.to(device) for k, v in inp.items()}
        with torch.no_grad():
            h = m(**inp).last_hidden_state.float()
        msk = inp["attention_mask"].float()
        e = (torch.einsum("ijk,ij->ik", h, msk) / msk.sum(1, keepdim=True)).cpu().numpy().astype(np.float16)
        assert np.isfinite(e).all()
        for (hsh, _), v in zip(items[i:j], e):
            out[hsh] = v
        i = j
        if len(out) % 2000 < (j - i):
            print(f"  ESM-C {len(out)}/{len(items)}", flush=True)
    del m
    return out


def bacformer_vectors(order, emb, device):
    from transformers import AutoModel
    m = AutoModel.from_pretrained(BFL, revision=BFL_REV, trust_remote_code=True).eval().to(device)
    rows = {}
    for acc, prots in order.items():
        contig = np.array([c for c, _ in prots])
        x = np.stack([emb[h] for _, h in prots]).astype(np.float32)
        tot, n = np.zeros(x.shape[1], dtype=np.float64), 0
        for s in range(0, len(x), 6000):
            pe = torch.tensor(x[s:s + 6000]).unsqueeze(0).to(device)
            cid = torch.tensor(np.clip(contig[s:s + 6000] - contig[s], 0, 999)).unsqueeze(0).to(device)
            with torch.no_grad():
                h = m(protein_embeddings=pe, contig_ids=cid, attention_mask=torch.ones_like(cid),
                      return_dict=True).last_hidden_state[0].float()
            assert torch.isfinite(h).all(), f"{acc}: non-finite"
            tot += h.sum(0).cpu().numpy()
            n += h.shape[0]
        rows[acc] = (tot / n).astype(np.float32)
        print(f"  {acc}: {len(prots)} proteins -> vector", flush=True)
    return rows


def main():
    accs = [a.strip() for a in open(sys.argv[1]) if a.strip()]
    out_path = Path(sys.argv[2])
    device = sys.argv[sys.argv.index("--device") + 1] if "--device" in sys.argv else \
        ("mps" if torch.backends.mps.is_available() else "cpu")
    print(f"{len(accs)} genomes on {device}", flush=True)
    order, seqs = {}, {}
    for a in accs:
        order[a], uq = proteins(a)
        seqs.update(uq)
    uniq = {h: seqs[h] for _, h in {p for v in order.values() for p in v}}
    print(f"{sum(len(v) for v in order.values())} proteins, {len(uniq)} unique", flush=True)
    emb = esmc_embed(uniq, device)
    rows = bacformer_vectors(order, emb, device)
    df = pd.DataFrame.from_dict(rows, orient="index")
    df.columns = [f"bacformer_large_mean_{i}" for i in range(df.shape[1])]
    df.to_parquet(out_path)
    print(f"written {out_path} {df.shape}")


if __name__ == "__main__":
    main()
