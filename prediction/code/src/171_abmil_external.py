"""External macro AUPRC for attention-based MIL (arm D, Table 1), the one Table 1 row without it.

As for the other external rows (src/22_external_eval.py, src/138_external_missing_rows.py): the model is trained on all
358 labelled genomes with the procedure of src/06b_arm_D_abmil.py and applied unchanged to the 562 external genomes.
The bag prefilter (all proteins or only proteins with an energy-panel KOfam hit) is chosen by the same inner grouped
CV as in each outer fold of arm D; the test model is the same 5-seed ensemble.
Per-protein ESM-2 t12 35M embeddings of the external proteins come from the Kaggle run that produced the external
ESM-2 genome vectors (your-kaggle-user/redox-esm2-bacformer-base-external, output downloaded to
data/kaggle/base_external_protein_esm2/). If the prefilter is selected, the external proteins are scanned with the
energy-panel HMMs and KofamScan thresholds (src/04_kofam_scan.scan), reusing the training hits for shared sequences.
Output: results/external_abmil.csv (per trait), results/external_abmil_predictions.csv
"""
import glob
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import PROC, RES, inner_folds, load_labelled, safe_ap  # noqa: E402

D = importlib.import_module("06b_arm_D_abmil")
EXT = ROOT / "data/external"
SHARDS = ROOT / "data/kaggle/base_external_protein_esm2"


def external_bags(ids, prefilter):
    emb = {}
    for f in sorted(glob.glob(str(SHARDS / "**/esm2_shard_*.npz"), recursive=True)):
        z = np.load(f)
        emb.update(zip(z["hashes"].tolist(), z["emb"]))
    prots = {a: pd.read_parquet(EXT / "proteins" / f"{a}.parquet", columns=["seq_hash"]).seq_hash for a in ids}
    need = {h for s in prots.values() for h in s}
    miss = need - set(emb)
    print(f"external proteins {len(need)}, embedded {len(need) - len(miss)}, missing {len(miss)}", flush=True)
    if miss:
        # sequences shared with the training genomes were embedded locally
        ge = importlib.import_module("07_genome_embeddings")
        index, E = ge.load_protein_embeddings()
        pos = index.get_indexer(list(miss))
        for h, p in zip(miss, pos):
            if p >= 0:
                emb[h] = E[p]
        miss = need - set(emb)
        assert not miss, f"{len(miss)} external proteins without an embedding"
    keep_hash = None
    if prefilter:
        f = EXT / "kofam_panel_hits_external.parquet"
        if not f.exists():
            ko = importlib.import_module("04_kofam_scan")
            new = pd.read_parquet(EXT / "unique_sequences_new.parquet")
            _, panel = ko.build_panel()
            hmms, thr, stype = ko.load_hmms(panel)
            hits = ko.scan(hmms, thr, stype, new.seq_hash.tolist(), new.sequence.tolist())
            hits.to_parquet(f)
        keep_hash = set(pd.read_parquet(f, columns=["seq_hash"]).seq_hash) | set(
            pd.read_parquet(PROC / "kofam_hits.parquet", columns=["seq_hash"]).seq_hash)
    bags = []
    for a in ids:
        x = np.stack([emb[h] for h in prots[a]]).astype(np.float32)
        if keep_hash is not None:
            k = prots[a].isin(keep_hash).to_numpy()
            x = x[k] if k.any() else x[:1]
        bags.append(x)
    return bags


def main():
    lab, traits, Y = load_labelled()
    bags = D.load_bags(lab)
    tr = np.arange(len(lab))
    groups = lab["group_genus"].to_numpy()
    strat = Y[:, int(np.argmin(Y.sum(0)))]
    folds = inner_folds(strat, groups)
    score = {}
    for ci, cfg in enumerate(D.CONFIGS):
        P = np.zeros((len(tr), len(traits)))
        for itr, iva in folds:
            P[iva], _ = D.train_predict(bags[cfg["prefilter"]], Y, tr[itr], tr[iva], seed=ci)
        score[ci] = float(np.nanmean([safe_ap(Y[:, j], P[:, j]) for j in range(len(traits))]))
        print("inner CV", cfg, round(score[ci], 4), flush=True)
    bi = max(score, key=lambda c: (score[c], -c))
    cfg = D.CONFIGS[bi]
    models = [D.train_model(bags[cfg["prefilter"]], Y, tr, seed=100 + s) for s in range(D.N_SEEDS_TEST)]

    ext = pd.read_csv(RES / "external_validation_predictions.csv")
    ext = ext[ext.arm == "B2"][["genome_id", "trait", "y"]]
    ids = sorted(ext.genome_id.unique())
    xb = external_bags(ids, cfg["prefilter"])
    P = np.mean([D.predict(m, xb, np.arange(len(ids)))[0] for m in models], axis=0)
    S = pd.DataFrame(P, index=ids, columns=traits)
    ext["score"] = [S.at[g, t] for g, t in zip(ext.genome_id, ext.trait)]
    ext.assign(arm="Attention-based MIL", config=json.dumps(cfg)).to_csv(RES / "external_abmil_predictions.csv", index=False)
    per = ext.groupby("trait").apply(lambda g: safe_ap(g.y.to_numpy(), g.score.to_numpy())).rename("auprc").reset_index()
    per.to_csv(RES / "external_abmil.csv", index=False, float_format="%.4f")
    print(per.round(3).to_string(index=False))
    print("selected", cfg, "external macro AUPRC", round(per.auprc.mean(), 4))


if __name__ == "__main__":
    main()
