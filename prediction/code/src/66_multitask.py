"""Multi-task transfer from the pathway task to donors/acceptors (preregistration_trees_xai.md, Addendum R).

One network over the KO columns shared by both tasks, with 34 outputs (14 donor/acceptor + 20 pathway) and a masked
loss so every genome contributes only the labels it has. Arm ST is the same network trained on the 14 donor and
acceptor labels alone, which isolates the transfer from the architecture. Evaluation uses the frozen donor/acceptor
folds, and pathway genomes from a test fold's genus (or family) are removed from training.

usage: python src/66_multitask.py [--arms MT ST] [--seeds 3]
Output: results/trees_xai/da/oof/{MT,ST}.parquet
"""
import argparse
import importlib
import json

import numpy as np
import pandas as pd

a32 = importlib.import_module("32_pathway_arms")
tx = importlib.import_module("47_trees_xai")
from pathway_common import ko_matrix, load_pathways  # noqa: E402
from redox_common import LEVELS, best_mcc_threshold, inner_folds, load_features, load_labelled  # noqa: E402


class MaskedMLP(a32.MultiLabelMLP):
    """PB4's network with a masked loss: unobserved labels contribute nothing to the gradient."""

    def fit(self, X, Y, M):
        import torch
        from torch import nn
        torch.manual_seed(self.seed)
        self.dev = self.device()
        self.mu, sd = X.mean(0), X.std(0)
        self.sd = np.where(sd > 1e-8, sd, 1.0)
        Xt = torch.tensor((X - self.mu) / self.sd, dtype=torch.float32, device=self.dev)
        Yt = torch.tensor(Y, dtype=torch.float32, device=self.dev)
        Mt = torch.tensor(M, dtype=torch.float32, device=self.dev)
        pos = (Yt * Mt).sum(0).clamp(min=1)
        obs = Mt.sum(0).clamp(min=1)
        w = ((obs - pos) / pos).clamp(max=50)
        self.net = nn.Sequential(nn.Linear(X.shape[1], self.h), nn.ReLU(), nn.Dropout(self.p),
                                 nn.Linear(self.h, Y.shape[1])).to(self.dev)
        opt = torch.optim.AdamW(self.net.parameters(), lr=self.lr, weight_decay=1e-4)
        lossf = nn.BCEWithLogitsLoss(pos_weight=w, reduction="none")
        g = torch.Generator().manual_seed(self.seed)
        for _ in range(self.epochs):
            self.net.train()
            for idx in torch.randperm(len(Xt), generator=g).to(self.dev).split(128):
                opt.zero_grad()
                loss = lossf(self.net(Xt[idx]), Yt[idx]) * Mt[idx]
                (loss.sum() / Mt[idx].sum().clamp(min=1)).backward()
                opt.step()
        return self


def build():
    """Shared KO matrix and masked label matrix for the union of both labelled sets."""
    da, dtr, dY = load_labelled()
    pw, ptr, pY = load_pathways()
    Xd = load_features("energy_panel", da.genome_id)
    Xp = ko_matrix(pw.genome_id)
    cols = [c for c in Xd.columns if c in Xp.columns]
    Xd, Xp = Xd[cols], Xp[cols]
    traits = list(dtr) + list(ptr)
    nd, npw = len(da), len(pw)
    Y = np.zeros((nd + npw, len(traits)), dtype=np.float32)
    M = np.zeros_like(Y)
    Y[:nd, :len(dtr)], M[:nd, :len(dtr)] = dY, 1.0
    Y[nd:, len(dtr):], M[nd:, len(dtr):] = pY, 1.0
    X = np.vstack([Xd.to_numpy(), Xp.to_numpy()]).astype(np.float32)
    keep = X.std(axis=0) > 0
    return da, pw, dtr, traits, X[:, keep], Y, M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="+", default=["MT", "ST"])
    ap.add_argument("--seeds", type=int, default=3)
    args = ap.parse_args()
    da, pw, dtr, traits, X, Y, M = build()
    nd, nda = len(da), len(dtr)
    print(f"{nd} donor/acceptor + {len(pw)} pathway genomes, {X.shape[1]} shared KO features, {len(traits)} labels",
          flush=True)
    out = tx.OUT / "da" / "oof"
    out.mkdir(parents=True, exist_ok=True)
    for arm in args.arms:
        cols = slice(0, len(traits)) if arm == "MT" else slice(0, nda)
        rows = []
        for level in LEVELS:
            fold = da[f"fold_{level}"].to_numpy()
            groups = da[f"group_{level}"].to_numpy()
            rank = "genus" if level == "genus" else "family"
            for k in sorted(np.unique(fold)):
                tr, te = np.where(fold != k)[0], np.where(fold == k)[0]
                itr, iva = inner_folds(Y[tr, int(np.argmin(Y[tr, :nda].sum(0)))], groups[tr])[0]
                if arm == "MT":
                    # Block the auxiliary pathway genomes against BOTH the outer test fold and the inner
                    # validation fold, by genome id and by lineage name under either taxonomy (the two tables
                    # come from different vintages: 16 genera and 30 families disagree, and 174 genomes are
                    # shared outright). Audit of 2026-09-21, finding 3.
                    def allowed(block_rows):
                        names = set()
                        for col in (rank, f"gtdb_{rank}", "group_genus", "group_family"):
                            if col in da.columns:
                                names |= set(da.loc[block_rows, col].dropna())
                        ids = set(da.loc[block_rows, "genome_id"])
                        keep = np.ones(len(pw), dtype=bool)
                        for col in (rank, f"gtdb_{rank}", "group_genus", "group_family"):
                            if col in pw.columns:
                                keep &= ~pw[col].isin(names).to_numpy()
                        keep &= ~pw.genome_id.isin(ids).to_numpy()
                        return nd + np.where(keep)[0]
                    idx = np.concatenate([tr, allowed(te)])
                    inner_idx = np.concatenate([tr[itr], allowed(np.concatenate([te, tr[iva]]))])
                else:
                    idx = tr
                    inner_idx = tr[itr]
                P_in = np.mean([MaskedMLP(hidden=256, seed=s).fit(X[inner_idx], Y[inner_idx, cols],
                                                                  M[inner_idx, cols]).predict_proba(X[tr[iva]])
                                for s in range(args.seeds)], axis=0)
                P_te = np.mean([MaskedMLP(hidden=256, seed=100 + s).fit(X[idx], Y[idx, cols], M[idx, cols])
                                .predict_proba(X[te]) for s in range(args.seeds)], axis=0)
                for j, t in enumerate(dtr):
                    thr = best_mcc_threshold(Y[tr[iva], j], P_in[:, j]) if 0 < Y[tr[iva], j].sum() < len(iva) else 0.5
                    for i, ix in enumerate(te):
                        rows.append((arm, level, t, da.genome_id[ix], int(k), Y[ix, j], float(P_te[i, j]), float(thr),
                                     json.dumps({"model": arm, "hidden": 256, "n_train": int(len(idx)), "device": MaskedMLP.device()}), np.nan))
                print(f"[{arm}/{level}] fold {k}: {len(idx)} training genomes", flush=True)
        df = pd.DataFrame(rows, columns=["arm", "level", "trait", "genome_id", "fold", "y", "prob", "threshold",
                                         "config", "inner_auprc"])
        df.to_parquet(out / f"{arm}.parquet", index=False)
        print(out / f"{arm}.parquet", df.groupby(["arm", "level"]).size().to_dict(), flush=True)


if __name__ == "__main__":
    main()
