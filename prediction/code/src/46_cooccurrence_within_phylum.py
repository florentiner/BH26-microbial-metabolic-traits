"""Re-test group co-occurrence (Test D of src/41) under the within-phylum null.

Test D shuffled whole genera, which removes within-genus replication but not confounding by phylum: two groups that
are both common in, say, Pseudomonadota will co-occur across genera for that reason alone. Here each genus is one
data point and one group's presence is permuted only among genera of the same phylum.

Output: results/jaccard/cooccurrence_within_phylum.csv
usage: python src/46_cooccurrence_within_phylum.py
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/jaccard"
SEED = 20260918


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    jc = load("jc", "src/40_jaccard_correlations.py")
    env = load("env", "src/45_environment_scan.py")
    rng = np.random.default_rng(SEED)
    X, groups, role, ranks, info = jc.load_html(ROOT / "jaccard_similarity.html")
    G = pd.DataFrame({n: X[k].max(axis=1) for n, k in groups.items()})
    genus = ranks["genus"].fillna("?")
    gl = G.groupby(genus).mean().ge(0.5).astype(float)          # one row per genus
    phylum = pd.DataFrame({"g": genus, "p": ranks["phylum"].fillna("?")}).groupby("g").p \
               .agg(lambda v: v.value_counts().index[0]).reindex(gl.index).to_numpy()
    names = list(G.columns)
    rows = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            xa, xb = gl[a].to_numpy(), gl[b].to_numpy()
            if min(xa.sum(), (1 - xa).sum(), xb.sum(), (1 - xb).sum()) < 5:
                continue
            _, p_gen = env.perm_stat(xa, xb, "binary", rng, n=5000)
            phi, p_phy = env.perm_within(xa, xb, "binary", phylum, rng, n=5000)
            rows.append(dict(group_a=a, group_b=b, n_genera=len(xa), n_genera_both=int((xa * xb).sum()),
                             phi_genus_level=phi, p_across_genera=p_gen, p_within_phylum=p_phy,
                             shares_kos=bool(set(groups[a]) & set(groups[b]))))
    d = pd.DataFrame(rows)
    d["q_across_genera"] = jc.bh(d.p_across_genera.to_numpy())
    d["q_within_phylum"] = jc.bh(d.p_within_phylum.to_numpy())
    d["verdict"] = np.select(
        [d.shares_kos & (d.q_across_genera < 0.05),
         (d.q_across_genera < 0.05) & (d.q_within_phylum < 0.05),
         d.q_across_genera < 0.05],
        ["artefact of shared KOs", "co-occur beyond phylum identity", "lineage-bound"], "not significant")
    d = d.sort_values(["verdict", "phi_genus_level"], ascending=[True, False])
    d.to_csv(OUT / "cooccurrence_within_phylum.csv", index=False, float_format="%.4f")
    for v in ("co-occur beyond phylum identity", "lineage-bound", "artefact of shared KOs"):
        s = d[d.verdict == v]
        print(f"\n=== {v}: {len(s)} ===")
        if len(s):
            print(s[["group_a", "group_b", "n_genera_both", "phi_genus_level", "q_across_genera",
                     "q_within_phylum"]].round(3).head(12).to_string(index=False))


if __name__ == "__main__":
    main()
