"""Rules extracted from the random forest itself, SIRUS-style (docs/preregistration_trees_xai.md, Addendum I).

For each trait the forest (CV-selected configuration, 3 seeds x 500 trees) is fitted on training genera; every tree
contributes its paths of depth 1 and 2 from the root; the 10 most frequent paths are kept. A rule states the
probability of the trait among genomes the path covers. Nested over genus folds, plus the external genomes (DA).

usage: python src/53_forest_rules.py --task {da,ox,pw}
Output: results/rulebook/<task>_forest_rules.csv
"""
import argparse
import importlib
import json
import re
from collections import Counter

import numpy as np
import pandas as pd

tx = importlib.import_module("47_trees_xai")
rb = importlib.import_module("51_rulebook")
OUT = tx.RES / "rulebook"
TOP = 10


def forest_paths(forests, cols):
    """Frequency of every depth-1 and depth-2 root path across all trees of all seed forests."""
    cnt, n_trees = Counter(), 0
    for f in forests:
        for est in f.estimators_:
            t = est.tree_
            n_trees += 1
            if t.children_left[0] == -1:
                continue
            f0 = cols[t.feature[0]]
            for v0, child in ((0, t.children_left[0]), (1, t.children_right[0])):   # left: x <= 0.5 -> absent
                cnt[((f0, v0),)] += 1
                if t.children_left[child] != -1:
                    f1 = cols[t.feature[child]]
                    if f1 != f0:
                        for v1 in (0, 1):
                            cnt[tuple(sorted(((f0, v0), (f1, v1))))] += 1
    return {p: c / n_trees for p, c in cnt.items()}


def cover(Xdf, path):
    m = np.ones(len(Xdf), bool)
    for g, v in path:
        m &= (Xdf[g].to_numpy() > 0) == bool(v)
    return m


def path_text(path, trait, p, base):
    cond = " AND ".join(f"{rb.MECH.get(g, (g, ''))[0]} ({g}) {'present' if v else 'absent'}" for g, v in path)
    direction = "raises" if p > base else "lowers"
    return f"IF {cond} THEN P({trait}) = {p:.0%} ({direction} it from {base:.0%})"


def rulebook_genes(task, trait):
    f = OUT / f"{task}_rules.csv"
    if not f.exists():
        return set()
    r = pd.read_csv(f)
    return set(re.findall(r"K\d{5}", " ".join(r[r.trait == trait].rule)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", choices=tx.TASKS, required=True)
    task = ap.parse_args().task
    lab, traits, Y, Xdf = tx.load_task(task)
    rb.register_task(Xdf, traits, Y)
    cols = list(Xdf.columns)
    X = Xdf.to_numpy(dtype=np.float32)
    fold = lab.fold_genus.to_numpy()
    oof = pd.read_parquet(tx.OUT / task / "oof/RF.parquet")
    oof = oof[oof.level == "genus"]
    ext = None
    if task == "da":
        ext = rb.load_external()
    rows = []
    for j, t in enumerate(traits):
        y = Y[:, j].astype(int)
        cfg = json.loads(oof[oof.trait == t].config.value_counts().index[0])
        base = float(y.mean())
        heldout, top_by_fold = {}, {}
        for k in sorted(np.unique(fold)):
            tr, te = fold != k, fold == k
            fr = forest_paths(tx.fit_all(cfg, X[tr], y[tr]), cols)
            top = [p for p, _ in sorted(fr.items(), key=lambda kv: -kv[1])[:TOP]]
            top_by_fold[k] = set(top)
            for p in top:
                m = cover(Xdf[te], p)
                a = heldout.setdefault(p, [0, 0])
                a[0] += int((m & (y[te] == 1)).sum())
                a[1] += int(m.sum())
        fr_all = forest_paths(tx.fit_all(cfg, X, y), cols)
        final = sorted(fr_all.items(), key=lambda kv: -kv[1])[:TOP]
        known = rulebook_genes(task, t)
        for rank, (p, freq) in enumerate(final, 1):
            m = cover(Xdf, p)
            p_train = float(y[m].mean()) if m.any() else np.nan
            kh, nh = heldout.get(p, [0, 0])
            ph = kh / nh if nh else np.nan
            lo, hi = rb.wilson(kh, nh) if nh else (np.nan, np.nan)
            genes = {g for g, _ in p}
            _, verdict, level = rb.interpret(p, 1 if p_train > base else 0, t)
            mech, _, _ = rb.interpret(p, 1, t)
            row = dict(task=task, trait=t, rank=rank, rule=path_text(p, t, p_train, base), genes=", ".join(sorted(genes)),
                       depth=len(p), tree_frequency=freq, stability=f"{sum(p in s for s in top_by_fold.values())}/5",
                       base_rate=base, p_in_sample=p_train, coverage_genomes=int(m.sum()),
                       p_heldout=ph, heldout_support=nh, heldout_ci_low=lo, heldout_ci_high=hi,
                       direction="raises" if p_train > base else "lowers",
                       informative_heldout=bool(nh and (lo > base or hi < base)),
                       shares_gene_with_rulebook=bool(genes & known), mechanism=mech, interpretation=verdict,
                       interpretation_level=level)
            if ext is not None:
                ex, el = ext
                col = rb.EXTERNAL_LABEL.get(t, t)
                if col in el.columns and genes <= set(ex.columns):
                    le = el[col].reindex(ex.index)
                    ok = le.notna().to_numpy()
                    me = cover(ex, p) & ok
                    ke, ne = int((me & (le.fillna(-1).to_numpy() == 1)).sum()), int(me.sum())
                    be = float(le[ok].mean())
                    elo, ehi = rb.wilson(ke, ne) if ne else (np.nan, np.nan)
                    row.update(external_base_rate=be, p_external=ke / ne if ne else np.nan, external_support=ne,
                               external_ci_low=elo, external_ci_high=ehi,
                               informative_external=bool(ne and (elo > be or ehi < be)))
            rows.append(row)
        print(f"[{task}] {t}: top path freq {final[0][1]:.2f}; "
              f"{sum(r['shares_gene_with_rulebook'] for r in rows if r['trait'] == t)}/{TOP} share a gene with the rulebook")
    d = pd.DataFrame(rows)
    d.to_csv(OUT / f"{task}_forest_rules.csv", index=False, float_format="%.4f")
    first = d[d["rank"] == 1]
    print(f"\n[{task}] I-P1: most frequent forest path shares a gene with the rulebook for "
          f"{int(first.shares_gene_with_rulebook.sum())} of {len(first)} traits")
    show = ["trait", "rule", "tree_frequency", "stability", "p_heldout"] + (["p_external"] if "p_external" in d else [])
    print(d[d["rank"] <= 2][show].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
