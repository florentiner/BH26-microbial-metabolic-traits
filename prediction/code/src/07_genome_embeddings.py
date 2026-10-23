"""Genome-level representations from cached protein embeddings.

  esm2_mean / esm2_max : mean / max over all protein loci of ESM-2 t12 35M protein embeddings (arm C)
  bacformer_mean       : Bacformer-base (macwiatrak/bacformer-masked-complete-genomes, frozen) contextual
                         protein embeddings, mean over protein tokens (special tokens excluded) (arm E)

Bacformer input = proteins in genome-coordinate order, one list per contig (contig order from GFF),
exactly as bacformer.pp.protein_embeddings_to_inputs expects. The model has 6000 positions; genomes
whose token count (proteins + CLS + one SEP per contig + END) exceeds this are processed in consecutive
windows of whole contigs (contigs longer than a window are split) and the protein-token outputs pooled.

Output: data/processed/genome_features_{esm2_mean,esm2_max,bacformer_mean}.parquet (index genome_id)
        data/processed/bacformer_windows.tsv (n_windows per genome, for the methods section)
"""
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from bacformer.modeling import SPECIAL_TOKENS_DICT
from bacformer.pp.embed_prot_seqs import protein_embeddings_to_inputs
from transformers import AutoModel

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data/processed"
BACFORMER = "macwiatrak/bacformer-masked-complete-genomes"
BACFORMER_REV = "d1c461f1db75902e6502e5ba7588ea1a4ea96b29"
BACFORMER_MAG = "macwiatrak/bacformer-masked-MAG"  # sensitivity arm E_MAG (174/358 labelled assemblies are drafts)
BACFORMER_MAG_REV = "340a08cf5d5cd7d9f1ffac9075cde733d108a82b"
MAX_TOKENS = 6000


def load_protein_embeddings():
    hashes, embs = [], []
    for f in sorted((PROC / "esm2_t12_35M").glob("shard_*.npz")):
        z = np.load(f)
        hashes.append(z["hashes"])
        embs.append(z["emb"])
    hashes = np.concatenate(hashes)
    assert len(hashes) == len(set(hashes.tolist())), "duplicate hashes across shards"
    return pd.Index(hashes), np.concatenate(embs)


def windows(contig_lists, max_tokens=MAX_TOKENS):
    """Split a list of contigs (each a list of row indices) into windows fitting the position limit."""
    budget = max_tokens - 2  # CLS + END
    out, cur, used = [], [], 0
    for c in contig_lists:
        while c:
            room = budget - used - 1  # one SEP per contig
            if room <= 0:
                out.append(cur)
                cur, used = [], 0
                continue
            take = c[:room]
            cur.append(take)
            used += len(take) + 1
            c = c[room:]
    if cur:
        out.append(cur)
    return out


def load_bacformer(dev, mag=False):
    name, rev = (BACFORMER_MAG, BACFORMER_MAG_REV) if mag else (BACFORMER, BACFORMER_REV)
    return AutoModel.from_pretrained(name, revision=rev, trust_remote_code=True).eval().to(dev).half()


def bacformer_mean(prot, emb, model, dev):
    """Mean over protein tokens of frozen Bacformer-base contextual embeddings (windowed at 6000 tokens)."""
    assert (prot.pos.to_numpy() == np.arange(len(prot))).all()
    prot_tok = SPECIAL_TOKENS_DICT["PROT_EMB"]
    contigs = [list(g.index) for _, g in prot.groupby("contig_idx", sort=True)]
    wins = windows(contigs)
    acc_sum, acc_n = np.zeros(emb.shape[1], dtype=np.float64), 0
    for w in wins:
        inp = protein_embeddings_to_inputs([[emb[i] for i in c] for c in w], max_n_proteins=MAX_TOKENS,
                                           bacformer_model_type="base")
        assert inp["special_tokens_mask"].shape[1] <= MAX_TOKENS
        assert int((inp["special_tokens_mask"] == prot_tok).sum()) == sum(map(len, w))
        x = {k: (v.to(dev).half() if v.is_floating_point() else v.to(dev)) for k, v in inp.items()}
        with torch.no_grad():
            h = model(**x, return_dict=True).last_hidden_state[0].float()
        hp = h[x["special_tokens_mask"][0] == prot_tok].cpu().numpy()
        assert np.isfinite(hp).all(), "non-finite Bacformer output"
        acc_sum += hp.sum(0)
        acc_n += len(hp)
    assert acc_n == len(prot)
    return (acc_sum / acc_n).astype(np.float32), len(wins)


def genome_representations(prot, emb, model, dev, mag_model=None):
    """prot: ordered protein table of one genome; emb: float32 (n_loci, 480) aligned to prot rows.
    Returns dict of genome vectors and the number of Bacformer windows."""
    bf, n_win = bacformer_mean(prot, emb, model, dev)
    out = {"esm2_mean": emb.mean(0), "esm2_max": emb.max(0), "bacformer_mean": bf}
    if mag_model is not None:
        out["bacformer_mag_mean"] = bacformer_mean(prot, emb, mag_model, dev)[0]
    return out, n_win


def main():
    import sys
    ids = pd.read_csv(ROOT / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig")["genome_id"]
    if "--labelled-only" in sys.argv:  # evaluation needs only labelled genomes; unlabelled ones are added later
        lab = pd.read_csv(PROC / "labels.tsv", sep="\t")
        ids = lab.loc[lab.labelled, "genome_id"].reset_index(drop=True)
    index, E = load_protein_embeddings()
    dev = "mps" if torch.backends.mps.is_available() else "cpu"
    model, mag_model = load_bacformer(dev), load_bacformer(dev, mag=True)
    feats, win_rows = {}, []
    for n, acc in enumerate(ids):
        prot = pd.read_parquet(PROC / "proteins" / f"{acc}.parquet")
        pos = index.get_indexer(prot.seq_hash)
        assert (pos >= 0).all(), f"{acc}: {(pos < 0).sum()} proteins without embeddings"
        reps, n_win = genome_representations(prot, E[pos].astype(np.float32), model, dev, mag_model)
        for k, v in reps.items():
            feats.setdefault(k, []).append(v)
        win_rows.append((acc, len(prot), prot.contig_idx.nunique(), n_win))
        if n % 50 == 0:
            print(n, acc, len(prot), "proteins", n_win, "windows", flush=True)
    for name, rows in feats.items():
        df = pd.DataFrame(np.vstack(rows), index=pd.Index(ids, name="genome_id"))
        df.columns = [f"{name}_{i}" for i in range(df.shape[1])]
        df.to_parquet(PROC / f"genome_features_{name}.parquet")
    pd.DataFrame(win_rows, columns=["genome_id", "n_proteins", "n_contigs", "n_windows"]).to_csv(
        PROC / "bacformer_windows.tsv", sep="\t", index=False)
    print("done")


if __name__ == "__main__":
    main()
