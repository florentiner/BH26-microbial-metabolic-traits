"""Nested grouped CV for the pathway arms (docs/preregistration_pathways.md).

usage: python src/32_pathway_arms.py {nulls,PA,PB1,PB2,PB3,PB4} [--jobs N]
  nulls : prevalence and taxonomic group-mean nulls
  PA    : KEGG module completeness rule (no fitting; threshold by MCC on the outer training fold)
  PB1   : 852-KO energy panel      PB2 : full significant-KO profile      PB3 : full KO profile + module completeness
  PB4   : multi-label MLP on the full KO profile - one shared representation, all traits jointly (hypothesis P6)
Output: results/pathways/oof/<arm>.parquet
"""
import argparse
import json

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import matthews_corrcoef

from pathway_common import RES, energy_panel_subset, ko_matrix, load_pathways, module_completeness, rule_scores
from redox_common import (LEVELS, best_mcc_threshold, fit_models, nested_cv, prevalence_null,
                          taxonomic_nn_null)

LR_C = [0.01, 0.1, 1.0, 10.0]
OUT = RES / "pathways/oof"


def configs():
    return [{"model": "logreg", "C": c} for c in LR_C] + [{"model": "xgb", "max_depth": d} for d in (2, 4)]


def fit_predict(cfg, Xtr, ytr, Xte, seed=0):
    return np.mean([m.predict_proba(Xte)[:, 1] for m in fit_models(cfg, Xtr, ytr)], axis=0)


def run_rules(lab, traits, Y, level, S):
    fold = lab[f"fold_{level}"].to_numpy()
    rows = []
    for j, t in enumerate(traits):
        for k in sorted(np.unique(fold)):
            tr, te = np.where(fold != k)[0], np.where(fold == k)[0]
            ytr, s_tr = Y[tr, j], S[tr, j]
            cands = {">0": 1e-6, "0.5": 0.5, "0.8": 0.8, "1.0": 1.0}
            mcc = {n: matthews_corrcoef(ytr, s_tr >= v) for n, v in cands.items()}
            best = max(mcc, key=lambda n: (mcc[n], -cands[n]))
            for i in te:
                rows.append(("PA", level, t, lab.genome_id[i], int(k), Y[i, j], float(S[i, j]), cands[best],
                             json.dumps({"model": "module_rule", "theta": best}), np.nan))
    return pd.DataFrame(rows, columns=["arm", "level", "trait", "genome_id", "fold", "y", "prob", "threshold",
                                       "config", "inner_auprc"])



class MultiLabelMLP:
    """One shared hidden layer over the KO profile, one sigmoid output per trait (hypothesis P6).

    All traits are fitted jointly, so a rare trait borrows representation from the common ones. Missing labels are
    not possible here (closed world), so the loss is a plain masked-free BCE with per-trait positive weighting.
    """

    def __init__(self, hidden=256, dropout=0.3, lr=1e-3, epochs=60, seed=0):
        self.h, self.p, self.lr, self.epochs, self.seed = hidden, dropout, lr, epochs, seed

    @staticmethod
    def device():
        """Apple-silicon GPU (MPS) when available, else CPU; override with PB4_DEVICE=cpu for bitwise reproducibility.

        The device changes where the arithmetic runs, not the model, loss, optimiser or mini-batch order (drawn on the
        CPU). Dropout masks, however, come from the device's own random generator, so an MPS run is equivalent to a
        different random seed, not a bitwise copy of a CPU run; MPS kernels are also not bitwise deterministic.
        Measured on a test problem: mean |p_CPU - p_MPS| 0.016, AUPRC difference ~0.005 - seed-level variation.
        """
        import os
        import torch
        want = os.environ.get("PB4_DEVICE", "auto")
        if want == "cpu" or not torch.backends.mps.is_available():
            return "cpu"
        return "mps"

    def fit(self, X, Y):
        import torch
        from torch import nn
        torch.manual_seed(self.seed)
        self.dev = self.device()
        self.mu, sd = X.mean(0), X.std(0)
        # zero-variance columns (a gene absent from - or present in - every training genome, common under
        # genus-blocked folds) get scale 1, as in sklearn's StandardScaler. Dividing by sd + 1e-6 turned such a
        # gene into a +-1e6 input on any test genome that differed, corrupting predictions and attributions.
        self.sd = np.where(sd > 1e-8, sd, 1.0)
        Xt = torch.tensor((X - self.mu) / self.sd, dtype=torch.float32, device=self.dev)
        Yt = torch.tensor(Y, dtype=torch.float32, device=self.dev)
        pos = Yt.sum(0).clamp(min=1)
        w = ((len(Yt) - pos) / pos).clamp(max=50)
        self.net = nn.Sequential(nn.Linear(X.shape[1], self.h), nn.ReLU(), nn.Dropout(self.p),
                                 nn.Linear(self.h, Y.shape[1])).to(self.dev)
        opt = torch.optim.AdamW(self.net.parameters(), lr=self.lr, weight_decay=1e-4)
        lossf = nn.BCEWithLogitsLoss(pos_weight=w)
        n = len(Xt)
        g = torch.Generator().manual_seed(self.seed)
        for _ in range(self.epochs):
            self.net.train()
            for idx in torch.randperm(n, generator=g).to(self.dev).split(128):  # order drawn on CPU: same on any device
                opt.zero_grad()
                lossf(self.net(Xt[idx]), Yt[idx]).backward()
                opt.step()
        return self

    def predict_proba(self, X):
        import torch
        self.net.eval()
        with torch.no_grad():
            xt = torch.tensor((X - self.mu) / self.sd, dtype=torch.float32, device=self.dev)
            return torch.sigmoid(self.net(xt)).cpu().numpy()


def run_mlp(arm, X, lab, traits, Y, level, seeds=3):
    """Nested CV for the multi-label arm: hidden size selected by inner grouped CV on macro AUPRC over all traits."""
    from redox_common import inner_folds, safe_ap
    fold = lab[f"fold_{level}"].to_numpy()
    groups = lab[f"group_{level}"].to_numpy()
    grid = [{"model": "mlp", "hidden": h, "dropout": 0.3} for h in (128, 256, 512)]
    rows = []
    for k in sorted(np.unique(fold)):
        tr, te = np.where(fold != k)[0], np.where(fold == k)[0]
        folds = inner_folds(Y[tr, int(np.argmin(Y[tr].sum(0)))], groups[tr])
        scores = {}
        for ci, cfg in enumerate(grid):
            oof = np.zeros_like(Y[tr])
            for itr, iva in folds:
                m = MultiLabelMLP(hidden=cfg["hidden"], dropout=cfg["dropout"], seed=ci).fit(X[tr][itr], Y[tr][itr])
                oof[iva] = m.predict_proba(X[tr][iva])
            scores[ci] = np.nanmean([safe_ap(Y[tr][:, j], oof[:, j]) for j in range(len(traits))])
        bi = max(scores, key=lambda c: (scores[c], -c))
        cfg = grid[bi]
        # threshold per trait from the inner OOF of the selected config, refitted to keep memory flat
        oof = np.zeros_like(Y[tr])
        for itr, iva in folds:
            oof[iva] = MultiLabelMLP(hidden=cfg["hidden"], seed=bi).fit(X[tr][itr], Y[tr][itr]).predict_proba(X[tr][iva])
        thr = [best_mcc_threshold(Y[tr][:, j], oof[:, j]) for j in range(len(traits))]
        P = np.mean([MultiLabelMLP(hidden=cfg["hidden"], seed=100 + s).fit(X[tr], Y[tr]).predict_proba(X[te])
                     for s in range(seeds)], axis=0)
        for j, t in enumerate(traits):
            for i, ix in enumerate(te):
                rows.append((arm, level, t, lab.genome_id[ix], int(k), Y[ix, j], float(P[i, j]), float(thr[j]),
                             json.dumps({**cfg, "device": MultiLabelMLP.device()}), float(scores[bi])))
    return pd.DataFrame(rows, columns=["arm", "level", "trait", "genome_id", "fold", "y", "prob", "threshold",
                                       "config", "inner_auprc"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm")
    ap.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args()
    lab, traits, Y = load_pathways()
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"{len(lab)} genomes, {len(traits)} traits, "
          f"{lab.group_genus.nunique()} genus groups, {lab.group_family.nunique()} family groups", flush=True)

    if args.arm == "nulls":
        out = [f(lab, traits, Y, lvl) for lvl in LEVELS for f in (prevalence_null, taxonomic_nn_null)]
    elif args.arm == "PA":
        S = rule_scores(lab, traits)
        out = [run_rules(lab, traits, Y, lvl, S) for lvl in LEVELS]
    elif args.arm == "PB4":
        X = ko_matrix(lab.genome_id).to_numpy()
        X = X[:, X.std(axis=0) > 0]
        print("feature matrix", X.shape, flush=True)
        out = [run_mlp(args.arm, X, lab, traits, Y, lvl) for lvl in LEVELS]
    else:
        if args.arm == "PB1":
            X = energy_panel_subset(lab.genome_id).to_numpy()
        elif args.arm == "PB2":
            X = ko_matrix(lab.genome_id).to_numpy()
        elif args.arm == "PB3":
            X = np.hstack([ko_matrix(lab.genome_id).to_numpy(), module_completeness(lab.genome_id).to_numpy()])
        else:
            raise SystemExit(f"unknown arm {args.arm}")
        X = X[:, X.std(axis=0) > 0]  # label-free constant-column drop
        print("feature matrix", X.shape, flush=True)
        out = Parallel(n_jobs=args.jobs)(
            delayed(nested_cv)(args.arm, X, lab, [traits[j]], Y[:, [j]], configs(), fit_predict, lvl)
            for lvl in LEVELS for j in range(len(traits)))
    df = pd.concat(out, ignore_index=True)
    df.to_parquet(OUT / f"{args.arm}.parquet", index=False)
    print(OUT, df.groupby(["arm", "level"]).size().to_dict())


if __name__ == "__main__":
    main()
