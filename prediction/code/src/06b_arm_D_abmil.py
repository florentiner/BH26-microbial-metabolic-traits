"""Arm D: gated-attention multiple-instance learning (ABMIL) over ESM-2 protein embeddings.

One bag per genome (instances = protein loci). Shared trunk, one gated-attention pooling and one
linear head per trait, class-weighted BCE over the 14 traits.

Declared budget (reduced from protocol §05 for local CPU compute):
  configs  : bag prefilter {off, on} x trunk width 128 x instance dropout 0.2 x entropy bonus 0
             prefilter "on" = keep only proteins with a KOfam hit in the energy-metabolism panel
             (label-free, fixed before training). Note: this injects curated annotation knowledge, so
             "D with prefilter" vs C is not a pure attention-vs-mean-pooling contrast.
  selection: inner grouped CV (4 folds, 1 seed, stratified on the rarest trait) by macro AUPRC; per-trait MCC
             thresholds from inner out-of-fold predictions of the selected config's 5-seed ensemble
  test     : 5-seed ensemble (mean probability) refitted on the outer training fold
  training : fixed 30 epochs, AdamW lr 1e-3, wd 1e-2, batch 16 bags; no early stopping (no validation
             data is touched).
Output: results/oof/D.parquet ; results/attention_collapse_D.csv (max attention weight diagnostics)
"""
import json
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

from redox_common import LEVELS, PROC, RES, best_mcc_threshold, inner_folds, load_labelled, safe_ap

torch.set_num_threads(6)
EPOCHS, LR, WD, BATCH, HIDDEN, ATT, INST_DROP = 30, 1e-3, 1e-2, 16, 128, 64, 0.2
CONFIGS = [{"model": "abmil", "prefilter": False}, {"model": "abmil", "prefilter": True}]
# audit round 4 (C-1): the selected config was prefilter=on in 7/10 fold-levels, which injects curated KO knowledge,
# so arm D as selected is not a clean attention-vs-mean-pooling contrast. `--prefilter off` runs the pure MIL arm
# (arm name D_nopf) for the protocol's Q3.
ARM = "D"
if "--prefilter" in sys.argv and sys.argv[sys.argv.index("--prefilter") + 1] == "off":
    CONFIGS = [{"model": "abmil", "prefilter": False}]
    ARM = "D_nopf"
N_SEEDS_TEST = 5


class ABMIL(nn.Module):
    def __init__(self, d_in, n_traits):
        super().__init__()
        self.trunk = nn.Sequential(nn.LayerNorm(d_in), nn.Linear(d_in, HIDDEN), nn.GELU(), nn.Dropout(0.1))
        self.V = nn.Linear(HIDDEN, ATT)
        self.U = nn.Linear(HIDDEN, ATT)
        self.w = nn.Linear(ATT, n_traits)
        self.head_w = nn.Parameter(torch.randn(n_traits, HIDDEN) * 0.01)
        self.head_b = nn.Parameter(torch.zeros(n_traits))

    def forward(self, x, mask):  # x (B,N,D) mask (B,N) bool
        h = self.trunk(x)
        a = self.w(torch.tanh(self.V(h)) * torch.sigmoid(self.U(h)))  # (B,N,T)
        a = a.masked_fill(~mask[..., None], float("-inf")).softmax(dim=1)
        z = torch.einsum("bnt,bnh->bth", a, h)
        return (z * self.head_w).sum(-1) + self.head_b, a


def make_batch(bags, idx, train, rng):
    xs = []
    for i in idx:
        b = bags[i]
        if train and len(b) > 1:
            keep = rng.random(len(b)) >= INST_DROP
            if keep.sum() == 0:
                keep[rng.integers(len(b))] = True
            b = b[keep]
        xs.append(b)
    n = max(len(b) for b in xs)
    X = np.zeros((len(xs), n, xs[0].shape[1]), dtype=np.float32)
    M = np.zeros((len(xs), n), dtype=bool)
    for k, b in enumerate(xs):
        X[k, :len(b)] = b
        M[k, :len(b)] = True
    return torch.from_numpy(X), torch.from_numpy(M)


def train_model(bags, Y, tr, seed):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = ABMIL(bags[0].shape[1], Y.shape[1])
    pos = Y[tr].sum(0)
    pw = torch.tensor((len(tr) - pos) / np.maximum(pos, 1), dtype=torch.float32)
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=pw)
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
    Yt = torch.tensor(Y, dtype=torch.float32)
    for _ in range(EPOCHS):
        model.train()
        order = rng.permutation(tr)
        for s in range(0, len(order), BATCH):
            idx = order[s:s + BATCH]
            X, M = make_batch(bags, idx, True, rng)
            loss = loss_fn(model(X, M)[0], Yt[idx])
            opt.zero_grad()
            loss.backward()
            opt.step()
    return model.eval()


def predict(model, bags, te):
    probs, maxatt = [], []
    with torch.no_grad():
        for s in range(0, len(te), BATCH):
            X, M = make_batch(bags, te[s:s + BATCH], False, None)
            logits, a = model(X, M)
            probs.append(torch.sigmoid(logits).numpy())
            maxatt.append(a.amax(dim=1).numpy())
    return np.vstack(probs), np.vstack(maxatt)


def train_predict(bags, Y, tr, te, seed):
    return predict(train_model(bags, Y, tr, seed), bags, te)


def ensemble_thresholds(bags, Y, tr, folds):
    """Per-trait MCC thresholds from inner out-of-fold predictions of the same 5-seed ensemble used at test time."""
    P = np.zeros((len(tr), Y.shape[1]))
    for itr, iva in folds:
        P[iva] = np.mean([train_predict(bags, Y, tr[itr], tr[iva], seed=100 + s)[0] for s in range(N_SEEDS_TEST)],
                         axis=0)
    return [best_mcc_threshold(Y[tr, j], P[:, j]) for j in range(Y.shape[1])]


def load_bags(lab):
    from importlib import import_module
    ge = import_module("07_genome_embeddings")
    index, E = ge.load_protein_embeddings()
    hits = pd.read_parquet(PROC / "kofam_hits.parquet", columns=["seq_hash"])
    panel_hashes = set(hits.seq_hash)
    full, filt = [], []
    for acc in lab.genome_id:
        prot = pd.read_parquet(PROC / "proteins" / f"{acc}.parquet", columns=["seq_hash"])
        pos = index.get_indexer(prot.seq_hash)
        assert (pos >= 0).all()
        emb = E[pos].astype(np.float32)
        full.append(emb)
        keep = prot.seq_hash.isin(panel_hashes).to_numpy()
        filt.append(emb[keep] if keep.any() else emb[:1])
    return {False: full, True: filt}


def main():
    lab, traits, Y = load_labelled()
    bags = load_bags(lab)
    print("median bag size full/prefiltered:", np.median([len(b) for b in bags[False]]),
          np.median([len(b) for b in bags[True]]), flush=True)
    rows, diag = [], []
    for level in LEVELS:
        fold = lab[f"fold_{level}"].to_numpy()
        groups = lab[f"group_{level}"].to_numpy()
        for k in sorted(np.unique(fold)):
            t0 = time.time()
            tr, te = np.where(fold != k)[0], np.where(fold == k)[0]
            # inner folds stratified on the rarest trait present in the training fold
            strat = Y[tr, int(np.argmin(Y[tr].sum(0)))]
            folds = inner_folds(strat, groups[tr])
            oof = {}
            for ci, cfg in enumerate(CONFIGS):
                P = np.zeros((len(tr), len(traits)))
                for itr, iva in folds:
                    P[iva], _ = train_predict(bags[cfg["prefilter"]], Y, tr[itr], tr[iva], seed=ci)
                oof[ci] = P
            score = {ci: np.nanmean([safe_ap(Y[tr, j], oof[ci][:, j]) for j in range(len(traits))]) for ci in oof}
            bi = max(score, key=lambda c: (score[c], -c))
            thr = ensemble_thresholds(bags[CONFIGS[bi]["prefilter"]], Y, tr, folds)
            preds, atts = zip(*[train_predict(bags[CONFIGS[bi]["prefilter"]], Y, tr, te, seed=100 + s)
                                for s in range(N_SEEDS_TEST)])
            p_te = np.mean(preds, axis=0)
            att = np.mean(atts, axis=0)
            for j, t in enumerate(traits):
                for r, i in enumerate(te):
                    rows.append((ARM, level, t, lab.genome_id[i], int(k), Y[i, j], float(p_te[r, j]), thr[j],
                                 json.dumps(CONFIGS[bi], sort_keys=True), score[bi]))
                diag.append(dict(level=level, fold=int(k), trait=t, median_max_attention=float(np.median(att[:, j])),
                                 frac_bags_max_att_gt_0_9=float(np.mean(att[:, j] > 0.9)),
                                 seed_sd_mean=float(np.mean(np.std([p[:, j] for p in preds], axis=0)))))
            print(f"[D/{level}] fold {k}: inner macro AUPRC {score} -> {CONFIGS[bi]} ({time.time() - t0:.0f}s)",
                  flush=True)
    df = pd.DataFrame(rows, columns=["arm", "level", "trait", "genome_id", "fold", "y", "prob", "threshold",
                                     "config", "inner_auprc"])
    df.to_parquet(RES / "oof" / f"{ARM}.parquet", index=False)
    pd.DataFrame(diag).to_csv(RES / f"attention_collapse_{ARM}.csv", index=False)


if __name__ == "__main__":
    main()
