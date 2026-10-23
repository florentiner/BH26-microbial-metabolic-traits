"""Reproducibility checks requested by audit round 1 (item M3). Run after 03, 04 (re-run) and 07.

1. KOfam: re-run with the frozen code (data/processed/kofam_hits.parquet) vs first run (*.run1.*): identical.
2. ESM-2: 500 random stored protein embeddings recomputed with the current 03_embed_esm2.py code: cosine ~ 1.
3. Protein order: every genome's protein table is sorted by (contig order, start) and contains each CDS locus once.
4. Frozen split hash unchanged.
Output: results/qc/reproducibility_checks.json (script exits non-zero if any check fails)
"""
import importlib
import json
import sys

import numpy as np
import pandas as pd
import torch

from redox_common import PROC, RES, sha256

emb_mod = importlib.import_module("03_embed_esm2")
ge = importlib.import_module("07_genome_embeddings")


def main():
    res = {}
    h1 = pd.read_parquet(PROC / "kofam_hits.run1.parquet").sort_values(["seq_hash", "KO"]).reset_index(drop=True)
    h2 = pd.read_parquet(PROC / "kofam_hits.parquet").sort_values(["seq_hash", "KO"]).reset_index(drop=True)
    m1 = pd.read_csv(PROC / "ko_matrix_panel.run1.tsv", sep="\t", index_col=0)
    m2 = pd.read_csv(PROC / "ko_matrix_panel.tsv", sep="\t", index_col=0)
    res["kofam_rerun_identical"] = bool(h1.equals(h2) and m1.equals(m2))
    res["kofam_hits_n"] = [len(h1), len(h2)]

    index, E = ge.load_protein_embeddings()
    u = pd.read_parquet(PROC / "unique_sequences.parquet").set_index("seq_hash")
    res["esm2_all_unique_sequences_embedded"] = bool(len(index) == len(u) and index.isin(u.index).all())
    rng = np.random.default_rng(0)
    pick = rng.choice(len(index), 500, replace=False)
    dev = "mps" if torch.backends.mps.is_available() else "cpu"
    new = emb_mod.embed_sequences(u.loc[index[pick], "sequence"].tolist(), dev).astype(np.float32)
    old = E[pick].astype(np.float32)
    cos = (new * old).sum(1) / np.linalg.norm(new, axis=1) / np.linalg.norm(old, axis=1)
    res["esm2_recompute_cosine_min"] = float(cos.min())
    res["esm2_recompute_ok"] = bool(cos.min() > 0.999)

    bad = []
    qc = pd.read_csv(PROC / "protein_table_qc.tsv", sep="\t")
    for acc in qc.genome_id:
        t = pd.read_parquet(PROC / "proteins" / f"{acc}.parquet")
        order_ok = (t.sort_values(["contig_idx", "start", "end"], kind="stable").pos.to_numpy() == t.pos.to_numpy()).all()
        if not (order_ok and t.locus.is_unique and (t.pos.to_numpy() == np.arange(len(t))).all()):
            bad.append(acc)
    res["protein_order_bad_genomes"] = bad
    expected = (PROC / "splits.sha256").read_text().split()[0]
    res["splits_hash_ok"] = sha256(PROC / "splits.tsv") == expected
    ok = res["kofam_rerun_identical"] and res["esm2_recompute_ok"] and not bad and res["splits_hash_ok"] and \
        res["esm2_all_unique_sequences_embedded"]
    res["all_ok"] = bool(ok)
    (RES / "qc").mkdir(parents=True, exist_ok=True)
    (RES / "qc" / "reproducibility_checks.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
