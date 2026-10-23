"""Robustness to incomplete genomes: RedoxKO versus the random forest on the MAG-isolate pairs of Table 5.

Proteins of both genomes of each pair come from data/additional/mags/faa (the Pyrodigal calls of src/88_real_mags.py).
Both are annotated with the 852-KO panel (KofamScan rules of training); RedoxKO is the released model
(results/models/B2), the forest is fitted on all 358 labeled genomes with min_samples_leaf and threshold chosen by
inner grouped CV. Agreement of calls between MAG and isolate is reported per completeness bin, with a paired
species-cluster bootstrap of the difference.
usage: python src/141_mag_rf_vs_redoxko.py
Output: results/additional/mag_rf_vs_redoxko.csv (+ data/additional/mags/ko_panel_pairs.parquet)
"""
import hashlib
import importlib
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox_common import best_mcc_threshold, inner_folds, load_features, load_labelled, safe_ap  # noqa: E402

pr = __import__("predict")
tx = importlib.import_module("47_trees_xai")
tx.TREE_JOBS = 8
MAG = ROOT / "data/additional/mags"


def read_faa(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf).rstrip("*")
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).rstrip("*")
    return [s for s in seqs.values() if s]


def main():
    pairs = pd.read_csv(MAG / "pairs.tsv", sep="\t")
    kp = MAG / "ko_panel_pairs.parquet"
    if kp.exists():
        X = pd.read_parquet(kp)
    else:
        tables, seqs = {}, {}
        for g in sorted(set(pairs.isolate) | set(pairs.mag)):
            f = MAG / "faa" / f"{g}.faa"
            if not f.exists():
                continue
            ss = read_faa(f)
            hs = [hashlib.md5(s.encode()).hexdigest() for s in ss]
            seqs.update(zip(hs, ss))
            tables[g] = pd.DataFrame({"seq_hash": hs})
        print(f"{len(tables)} genomes, {len(seqs)} unique proteins", flush=True)
        X = pr._isolated("ko_features", tables, seqs)
        X.to_parquet(kp)
    lab, traits, Y = load_labelled()
    Xt = load_features("energy_panel", lab.genome_id)
    keep = Xt.columns[Xt.std(axis=0) > 0]
    Xi, Xp = Xt[keep].to_numpy(), X.reindex(columns=keep).fillna(0).to_numpy(dtype=np.float32)
    bundle = joblib.load(ROOT / "results/models/B2/model.joblib")
    groups = lab.group_genus.to_numpy()
    calls = {"RedoxKO": {}, "RF": {}}
    for j, t in enumerate(traits):
        b = bundle[t]
        calls["RedoxKO"][t] = pd.Series((np.asarray(b["model"].predict_proba(X.reindex(columns=b["model"].columns).fillna(0)))
                                         >= b["threshold"]).astype(int), index=X.index)
        y, m = Y[:, j], ~np.isnan(Y[:, j])
        folds = inner_folds(y[m], groups[m])
        best = None
        for cfg in tx.CONFIGS["rf"]:
            oof = np.zeros(m.sum())
            for itr, iva in folds:
                oof[iva] = tx.fit_predict(cfg, Xi[m][itr], y[m][itr], Xi[m][iva])
            ap = safe_ap(y[m], oof)
            if best is None or ap > best[0]:
                best = (ap, cfg, oof)
        thr = best_mcc_threshold(y[m], best[2])
        calls["RF"][t] = pd.Series((tx.fit_predict(best[1], Xi[m], y[m], Xp) >= thr).astype(int), index=X.index)
        print("fitted", t, flush=True)
    pairs = pairs[pairs.isolate.isin(X.index) & pairs.mag.isin(X.index)].reset_index(drop=True)
    agree = {n: np.array([[int(c[t][r.isolate] == c[t][r.mag]) for t in traits] for r in pairs.itertuples()]).mean(axis=1)
             for n, c in calls.items()}
    rows = []
    rng = np.random.default_rng(20261021)
    for b in ("50-70", "70-90", "90-100"):
        s = pairs.bin.to_numpy() == b
        sp = pairs.species[s].to_numpy()
        a_b, a_r = agree["RedoxKO"][s], agree["RF"][s]
        u = np.unique(sp)
        idx = {v: np.where(sp == v)[0] for v in u}
        d = []
        for _ in range(2000):
            k = np.concatenate([idx[v] for v in rng.choice(u, len(u))])
            d.append(a_r[k].mean() - a_b[k].mean())
        rows.append(dict(bin=b, pairs=int(s.sum()), RedoxKO=a_b.mean(), RF=a_r.mean(), rf_minus_redoxko=a_r.mean() - a_b.mean(),
                         ci_low=np.percentile(d, 2.5), ci_high=np.percentile(d, 97.5)))
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results/additional/mag_rf_vs_redoxko.csv", index=False, float_format="%.4f")
    print(out.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
