"""Self-supervised pretraining on unlabelled KO profiles, then supervised fine-tuning (Addendum T).

The one deep-learning lever not yet tested: use the ~3 400 pathway genomes' KO profiles WITHOUT their labels to
pretrain an encoder (masked-feature reconstruction, a denoising autoencoder), then fine-tune encoder + head on the
358 donor/acceptor genomes. Arm SSL. Its control is ST (src/66_multitask.py): the identical network and features
trained from scratch, so SSL − ST isolates the value of the unlabelled data.

Leakage control is the one used for MT after the audit: pretraining genomes are blocked by genome id and by
lineage under both taxonomies against the outer test fold and the inner validation fold.

usage: PB4_DEVICE=mps python src/81_ssl_pretrain.py [--seeds 3] [--pre-epochs 40]
Output: results/trees_xai/da/oof/SSL.parquet
"""
import argparse
import importlib
import json

import numpy as np
import pandas as pd

mt = importlib.import_module("66_multitask")
tx = importlib.import_module("47_trees_xai")
from redox_common import LEVELS, best_mcc_threshold, inner_folds  # noqa: E402


class Pretrained:
    """Encoder pretrained by masked reconstruction, then fine-tuned with a per-trait weighted BCE head."""

    def __init__(self, hidden=256, dropout=0.3, seed=0, pre_epochs=40, ft_epochs=60, lr=1e-3):
        self.h, self.p, self.seed, self.pre_epochs, self.ft_epochs, self.lr = hidden, dropout, seed, pre_epochs, ft_epochs, lr

    def _dev(self):
        return mt.a32.MultiLabelMLP.device()

    def pretrain(self, X_unl):
        import torch
        from torch import nn
        torch.manual_seed(self.seed)
        dev = self._dev()
        self.mu, sd = X_unl.mean(0), X_unl.std(0)
        self.sd = np.where(sd > 1e-8, sd, 1.0)
        X = torch.tensor((X_unl - self.mu) / self.sd, dtype=torch.float32, device=dev)
        raw = torch.tensor(X_unl, dtype=torch.float32, device=dev)
        d = X.shape[1]
        self.enc = nn.Sequential(nn.Linear(d, self.h), nn.ReLU(), nn.Dropout(self.p)).to(dev)
        dec = nn.Linear(self.h, d).to(dev)
        opt = torch.optim.AdamW(list(self.enc.parameters()) + list(dec.parameters()), lr=self.lr, weight_decay=1e-4)
        lossf = nn.BCEWithLogitsLoss()
        g = torch.Generator().manual_seed(self.seed)
        for _ in range(self.pre_epochs):
            self.enc.train()
            for idx in torch.randperm(len(X), generator=g).to(dev).split(128):
                mask = (torch.rand(len(idx), d, generator=g) < 0.2).to(dev)         # hide 20% of the KOs
                xin = X[idx].masked_fill(mask, 0.0)
                opt.zero_grad()
                loss = lossf(dec(self.enc(xin))[mask], raw[idx][mask])          # reconstruct the hidden ones
                loss.backward()
                opt.step()
        return self

    def fit(self, X, Y):
        import torch
        from torch import nn
        dev = self._dev()
        torch.manual_seed(self.seed + 1)
        Xt = torch.tensor((X - self.mu) / self.sd, dtype=torch.float32, device=dev)
        Yt = torch.tensor(Y, dtype=torch.float32, device=dev)
        pos = Yt.sum(0).clamp(min=1)
        w = ((len(Yt) - pos) / pos).clamp(max=50)
        self.head = nn.Linear(self.h, Y.shape[1]).to(dev)
        opt = torch.optim.AdamW(list(self.enc.parameters()) + list(self.head.parameters()), lr=self.lr,
                                weight_decay=1e-4)
        lossf = nn.BCEWithLogitsLoss(pos_weight=w)
        g = torch.Generator().manual_seed(self.seed + 1)
        for _ in range(self.ft_epochs):
            self.enc.train()
            for idx in torch.randperm(len(Xt), generator=g).to(dev).split(128):
                opt.zero_grad()
                lossf(self.head(self.enc(Xt[idx])), Yt[idx]).backward()
                opt.step()
        return self

    def predict_proba(self, X):
        import torch
        self.enc.eval()
        with torch.no_grad():
            xt = torch.tensor((X - self.mu) / self.sd, dtype=torch.float32, device=self._dev())
            return torch.sigmoid(self.head(self.enc(xt))).cpu().numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--pre-epochs", type=int, default=40)
    args = ap.parse_args()
    da, pw, dtr, traits, X, Y, M = mt.build()
    nd, nda = len(da), len(dtr)
    print(f"{nd} labelled + {len(pw)} unlabelled-for-this-task genomes, {X.shape[1]} features", flush=True)
    rows = []
    for level in LEVELS:
        fold = da[f"fold_{level}"].to_numpy()
        groups = da[f"group_{level}"].to_numpy()
        rank = "genus" if level == "genus" else "family"
        for k in sorted(np.unique(fold)):
            tr, te = np.where(fold != k)[0], np.where(fold == k)[0]
            itr, iva = inner_folds(Y[tr, int(np.argmin(Y[tr, :nda].sum(0)))], groups[tr])[0]

            def allowed(block_rows):                 # same blocking as the repaired MT arm
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

            unl_outer = np.concatenate([tr, allowed(te)])                       # features only, no labels used
            unl_inner = np.concatenate([tr[itr], allowed(np.concatenate([te, tr[iva]]))])
            P_in = np.mean([Pretrained(seed=s, pre_epochs=args.pre_epochs).pretrain(X[unl_inner])
                            .fit(X[tr[itr]], Y[tr[itr], :nda]).predict_proba(X[tr[iva]]) for s in range(args.seeds)], axis=0)
            P_te = np.mean([Pretrained(seed=100 + s, pre_epochs=args.pre_epochs).pretrain(X[unl_outer])
                            .fit(X[tr], Y[tr, :nda]).predict_proba(X[te]) for s in range(args.seeds)], axis=0)
            for j, t in enumerate(dtr):
                thr = best_mcc_threshold(Y[tr[iva], j], P_in[:, j]) if 0 < Y[tr[iva], j].sum() < len(iva) else 0.5
                for i, ix in enumerate(te):
                    rows.append(("SSL", level, t, da.genome_id[ix], int(k), Y[ix, j], float(P_te[i, j]), float(thr),
                                 json.dumps({"model": "SSL", "hidden": 256, "n_pretrain": int(len(unl_outer)),
                                             "pre_epochs": args.pre_epochs, "device": mt.a32.MultiLabelMLP.device()}),
                                 np.nan))
            print(f"[SSL/{level}] fold {k}: pretrained on {len(unl_outer)} genomes", flush=True)
    df = pd.DataFrame(rows, columns=["arm", "level", "trait", "genome_id", "fold", "y", "prob", "threshold",
                                     "config", "inner_auprc"])
    out = tx.OUT / "da" / "oof"
    df.to_parquet(out / "SSL.parquet", index=False)
    print(out / "SSL.parquet", df.groupby(["arm", "level"]).size().to_dict())


if __name__ == "__main__":
    main()
