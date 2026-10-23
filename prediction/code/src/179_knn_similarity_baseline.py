"""Review round 3: a nearest-neighbour baseline that transfers labels from the most similar training genomes.

For each test genome, the k training genomes (other folds only, same frozen genus- or family-level folds as RedoxKO)
with the highest Jaccard similarity of KO presence give a similarity-weighted share of positives as the score.
Main variant: the 852-KO panel of RedoxKO (all 358 genomes), k = 5. Sensitivity: k = 1 and 10, and the full
KofamScan profile (306 genomes with a full profile). Paired difference RedoxKO minus kNN with a cluster bootstrap
over genera or families.
Output: results/additional/rf_paper/knn_baseline.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/additional/rf_paper"
SEED, B = 20261022, 2000


def jaccard(X):
    X = X.astype(np.float32)
    inter = X @ X.T
    n = X.sum(1)
    union = n[:, None] + n[None, :] - inter
    return np.where(union > 0, inter / np.maximum(union, 1), 0.0)


def knn_oof(S, ids, folds, Y, k):
    rows = []
    idx = {g: i for i, g in enumerate(ids)}
    for g, f in folds.items():
        i = idx[g]
        tr = np.array([idx[h] for h, fh in folds.items() if fh != f])
        s = S[i, tr]
        top = tr[np.argsort(-s)[:k]]
        w = S[i, top] + 1e-9
        for t in Y.columns:
            rows.append((g, t, float(np.dot(w, Y.iloc[top][t].to_numpy()) / w.sum())))
    return pd.DataFrame(rows, columns=["genome_id", "trait", "prob"])


def macro(d):
    return float(np.mean([average_precision_score(g.y, g.prob) for _, g in d.groupby("trait") if 0 < g.y.sum() < len(g)]))


def paired(a, b, cl):
    m = a.merge(b, on=["genome_id", "trait"], suffixes=("_a", "_b"))
    m["cl"] = cl.reindex(m.genome_id).astype(str).to_numpy()
    rng = np.random.default_rng(SEED)
    cls = m.cl.unique()
    ix = {c: np.where(m.cl.to_numpy() == c)[0] for c in cls}
    out = []
    for _ in range(B):
        x = m.iloc[np.concatenate([ix[c] for c in rng.choice(cls, len(cls))])]
        try:
            out.append(macro(x.rename(columns={"y_a": "y", "prob_a": "prob"})) -
                       macro(x.rename(columns={"y_a": "y", "prob_b": "prob"}).drop(columns="prob_a")))
        except Exception:
            pass
    return np.nanpercentile(out, [2.5, 97.5]).tolist()


def main():
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    lab = lab[lab.labelled == True].set_index("genome_id")  # noqa: E712
    rf = pd.read_parquet(ROOT / "results/trees_xai/da/oof/RF.parquet")
    traits = sorted(rf.trait.unique())
    Y = lab[traits].astype(int)
    panel = pd.read_csv(ROOT / "data/processed/ko_matrix_panel.tsv", sep="\t", index_col=0).reindex(lab.index).fillna(0)
    full = pd.read_parquet(ROOT / "data/processed/togo_kofam_internal.parquet", columns=["genome_id", "KO", "significant"])
    full = full[full.significant & full.genome_id.isin(lab.index)]
    fullm = pd.crosstab(full.genome_id, full.KO).clip(upper=1)
    res = {}
    for level in ("genus", "family"):
        r = rf[rf.level == level]
        folds = r.drop_duplicates("genome_id").set_index("genome_id").fold
        cl = lab[level]
        a = r[["genome_id", "trait", "y", "prob"]]
        out = {"RedoxKO": macro(a)}
        S = jaccard(panel.loc[folds.index].to_numpy())
        for k in (1, 5, 10):
            kn = knn_oof(S, list(folds.index), folds.to_dict(), Y.loc[folds.index], k).merge(a[["genome_id", "trait", "y"]])
            out[f"kNN_panel_k{k}"] = macro(kn)
            if k == 5:
                out["diff_k5"] = out["RedoxKO"] - out["kNN_panel_k5"]
                out["diff_k5_ci"] = paired(a, kn[["genome_id", "trait", "prob"]], cl)
        sub = folds[folds.index.isin(fullm.index)]
        Sf = jaccard(fullm.loc[sub.index].to_numpy())
        kf = knn_oof(Sf, list(sub.index), sub.to_dict(), Y.loc[sub.index], 5).merge(a[["genome_id", "trait", "y"]])
        out["full_profile_genomes"] = int(len(sub))
        out["kNN_full_k5"] = macro(kf)
        out["RedoxKO_same_genomes"] = macro(a[a.genome_id.isin(sub.index)])
        res[level] = out
    (OUT / "knn_baseline.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=lambda v: round(float(v), 3)))




def extra():
    """Own bootstrap intervals in cross-validation, and the external genomes (same genome-category pairs as RedoxKO)."""
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    lab = lab[lab.labelled == True].set_index("genome_id")  # noqa: E712
    rf = pd.read_parquet(ROOT / "results/trees_xai/da/oof/RF.parquet")
    traits = sorted(rf.trait.unique())
    Y = lab[traits].astype(int)
    panel = pd.read_csv(ROOT / "data/processed/ko_matrix_panel.tsv", sep="\t", index_col=0).reindex(lab.index).fillna(0)
    res = json.loads((OUT / "knn_baseline.json").read_text())
    for level in ("genus", "family"):
        r = rf[rf.level == level]
        folds = r.drop_duplicates("genome_id").set_index("genome_id").fold
        S = jaccard(panel.loc[folds.index].to_numpy())
        kn = knn_oof(S, list(folds.index), folds.to_dict(), Y.loc[folds.index], 5).merge(r[["genome_id", "trait", "y"]])
        kn["cl"] = lab[level].reindex(kn.genome_id).astype(str).to_numpy()
        rng = np.random.default_rng(SEED)
        cls = kn.cl.unique()
        ix = {c: np.where(kn.cl.to_numpy() == c)[0] for c in cls}
        bs = [macro(kn.iloc[np.concatenate([ix[c] for c in rng.choice(cls, len(cls))])]) for _ in range(B)]
        res[level]["kNN_panel_k5_ci"] = np.nanpercentile(bs, [2.5, 97.5]).tolist()
    # external genomes: neighbours among all 358 training genomes
    p = pd.read_csv(ROOT / "results/external_validation_missing_rows_predictions.csv").query("arm == 'Random forest'")
    yl = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t").set_index("genome_id")
    p = p[p.trait.isin(yl.columns)].copy()
    p["y"] = [yl.at[g, t] for g, t in zip(p.genome_id, p.trait)]
    p = p.dropna(subset=["y"])
    ref = pd.read_csv(ROOT / "results/external_validation_predictions.csv").query("arm == 'B2'")
    p = p.merge(ref[["genome_id", "trait"]], on=["genome_id", "trait"])   # same genome-category pairs as the paper
    ext = pd.read_csv(ROOT / "data/external/ko_matrix_panel.tsv", sep="\t", index_col=0)
    ext = ext.reindex(columns=panel.columns).fillna(0)
    Xt, Xe = panel.to_numpy(np.float32), ext.to_numpy(np.float32)
    inter = Xe @ Xt.T
    union = Xe.sum(1)[:, None] + Xt.sum(1)[None, :] - inter
    Sx = np.where(union > 0, inter / np.maximum(union, 1), 0)
    pos = {g: i for i, g in enumerate(ext.index)}
    sc = []
    for g, t in zip(p.genome_id, p.trait):
        s = Sx[pos[g]]
        top = np.argsort(-s)[:5]
        w = s[top] + 1e-9
        sc.append(float(np.dot(w, Y[t].to_numpy()[top]) / w.sum()))
    q = p.assign(prob=sc)
    res["external"] = {"kNN_panel_k5": macro(q), "RedoxKO_check": macro(p.rename(columns={"score": "prob"}))}
    (OUT / "knn_baseline.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: res[k] for k in ("external",)}, indent=1), res["genus"]["kNN_panel_k5_ci"], res["family"]["kNN_panel_k5_ci"])


if __name__ == "__main__":
    main() if len(sys.argv) == 1 else extra()
