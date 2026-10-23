"""Are the test sets autotrophs, like the training collection? And does RedoxKO work equally well on both kinds?

Autotrophy proxy (genome-based, label-free): a genome is autotroph-like when at least one carbon-fixation KEGG module
is at least 80% complete (Calvin cycle M00165, reductive citrate cycle M00173, Wood-Ljungdahl M00377,
3-hydroxypropionate bi-cycle M00376, hydroxypropionate-hydroxybutyrate cycle M00375, dicarboxylate-hydroxybutyrate
cycle M00374), or CO2 reduction to methane (M00567, the archaeal route of hydrogenotrophic methanogens), computed
from the 852-KO panel. Carbon-fixation genes alone do not prove autotrophic growth, so the
proxy is checked on the training collection, which was chosen as autotrophs.

Scores: RedoxKO (arm B2) on the 562 external genomes and on the 255 growth experiments of the 100 primary new species.
usage: python src/133_autotroph_subset.py
Output: results/additional/autotroph_subset.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pathway_common import completeness_from_matrix  # noqa: E402
from redox_common import PROC, load_labelled, safe_ap  # noqa: E402

CFIX = ["M00165", "M00173", "M00377", "M00376", "M00375", "M00374", "M00567"]
SEED, B = 20261021, 2000


def autotroph_like(X):
    md = pd.read_csv(PROC / "kegg_modules_energy.tsv", sep="\t")
    md = md[md.module.isin(CFIX)]
    C = completeness_from_matrix(X.to_numpy(), X.columns, md)
    return pd.Series(C.max(axis=1) >= 0.8, index=X.index), pd.DataFrame(C, index=X.index, columns=md.module)


def macro(d, col):
    au, ap = [], []
    for _, g in d.groupby("trait"):
        if 0 < g.y.sum() < len(g):
            au.append(roc_auc_score(g.y, g[col]))
            ap.append(safe_ap(g.y.to_numpy(), g[col].to_numpy()))
    return (np.mean(au) if au else np.nan), (np.mean(ap) if ap else np.nan), len(au)


def boot_diff(d, cluster):
    """Macro AUROC autotroph-like minus other, bootstrap over clusters (genera or papers)."""
    rng = np.random.default_rng(SEED)
    cl = d[cluster].unique()
    idx = {c: np.where(d[cluster].to_numpy() == c)[0] for c in cl}
    out = []
    for _ in range(B):
        s = d.iloc[np.concatenate([idx[c] for c in rng.choice(cl, len(cl))])]
        a, o = macro(s[s.auto], "score")[0], macro(s[~s.auto], "score")[0]
        if np.isfinite(a) and np.isfinite(o):
            out.append(a - o)
    return np.percentile(out, [2.5, 97.5])


def main():
    lab, traits, Y = load_labelled()
    train = pd.read_csv(PROC / "ko_matrix_panel.tsv", sep="\t", index_col=0).loc[lab.genome_id]
    ext_X = pd.read_csv(ROOT / "data/external/ko_matrix_panel.tsv", sep="\t", index_col=0)
    new_X = pd.read_parquet(ROOT / "validation_new_organisms/ko_panel_validation.parquet")
    cols = train.columns
    a_tr, _ = autotroph_like(train)
    a_ex, _ = autotroph_like(ext_X.reindex(columns=cols).fillna(0))
    a_nw, _ = autotroph_like(new_X.reindex(columns=cols).fillna(0))

    ext = pd.read_csv(ROOT / "results/external_validation_predictions.csv")
    ext = ext[ext.arm == "B2"].dropna(subset=["y"]).copy()
    lin = pd.read_csv(ROOT / "data/external/lineage.tsv", sep="\t")
    ext = ext.merge(lin[["genome_id", "ncbi_genus"]], on="genome_id", how="left")
    ext["auto"] = ext.genome_id.map(a_ex).fillna(False).astype(bool)

    new = pd.read_csv(ROOT / "validation_new_organisms/scored178_B2_balanced100.csv")
    new = new[~new.category.isin(["donor_other", "acceptor_other"])].rename(columns={"category": "trait", "expected": "y"})
    new["auto"] = new.accession.map(a_nw).fillna(False).astype(bool)

    rows = [dict(set="training (358)", genomes=len(a_tr), autotroph_like=float(a_tr.mean()))]
    for name, d, gid, cl in [("external genomes", ext, "genome_id", "ncbi_genus"), ("new species", new, "accession", "paper")]:
        share = d.drop_duplicates(gid).auto.mean()
        for sub, s in [("all", d), ("autotroph-like", d[d.auto]), ("other", d[~d.auto])]:
            au, ap, nt = macro(s, "score")
            rows.append(dict(set=name, subset=sub, genomes=s[gid].nunique(), tests=len(s), autotroph_like=share,
                             macro_auroc=au, macro_auprc=ap, categories=nt))
        lo, hi = boot_diff(d, cl)
        rows[-1].update(auroc_diff_auto_minus_other=macro(d[d.auto], "score")[0] - macro(d[~d.auto], "score")[0],
                        ci_low=lo, ci_high=hi)
    out = pd.DataFrame(rows)
    (ROOT / "results/additional").mkdir(parents=True, exist_ok=True)
    out.to_csv(ROOT / "results/additional/autotroph_subset.csv", index=False, float_format="%.4f")
    pd.set_option("display.width", 220)
    print(out.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
