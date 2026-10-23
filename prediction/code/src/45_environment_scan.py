"""Reference-style scan: which enzyme groups track which environmental conditions?

The preliminary observation "aerobic-respiration enzymes correlate with oxygen level" is one instance of a general
question: for each functional group x each recorded condition, is there an association beyond shared ancestry?
This script asks it for every condition Madin et al. 2020 record with usable coverage for these genomes.

Environmental conditions (primary family):  oxygen class, growth / optimum temperature, thermophily, optimum pH,
                                            doubling time
Cell traits (secondary family, reported separately): Gram stain, sporulation, motility

Method - analysis at the genus level, because genomes of one genus are not independent observations:
  * each genus contributes one data point: marker present if >= half its genomes carry it, trait = genus mean
    (numeric) or majority (binary);
  * numeric trait: Spearman rho between marker presence and trait (the rank-biserial correlation), plus the median
    trait value in marker+ vs marker- genera for interpretation;
  * binary trait: phi (MCC) between marker presence and trait;
  * p (across genera): trait permuted across all genera. Removes within-genus replication only - NOT confounding by
    deeper lineage: a marker found only in archaea will "correlate" with anything archaea share.
  * p (within phylum): trait permuted only among genera of the same phylum, so phylum and domain structure is kept
    intact. An association that survives this is not explained by which phylum an organism belongs to. This is the
    test that decides whether a correlation is a functional link rather than lineage identity.
  BH q within each family, separately for both p-values.

Output: results/jaccard/environment_scan.csv
usage: python src/45_environment_scan.py
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/jaccard"
SEED = 20260918
N_PERM = 10000
MIN_GENERA = 25          # genera with the trait recorded
MIN_MINORITY = 5         # genera on each side of the marker split

spec = importlib.util.spec_from_file_location("jc", ROOT / "src/40_jaccard_correlations.py")
jc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(jc)


# Refined marker sets from the definition tests (src/41-43): the diagnostic gene or the physiologically correct
# subset of a group. Scanned alongside the groups as defined, and reported as a separate marker family.
REFINED = {
    "aa3 + cbb3 oxidases (respiratory O2 reduction)": ["K02274", "K02275", "K02276", "K00404", "K00405", "K00406"],
    "aa3 cytochrome c oxidase only": ["K02274", "K02275", "K02276"],
    "cytochrome bd only (O2 scavenging)": ["K00425", "K00426"],
    "amoB only (ammonia monooxygenase)": ["K10945"],
    "hao only (hydroxylamine oxidoreductase)": ["K10535"],
    "mcrA only (methyl-CoM reductase)": ["K00399"],
    "dsrB only (dissimilatory sulfite reductase)": ["K11181"],
    "narGH / nxrAB only": ["K00370", "K00371"],
    "nosZ only (N2O reductase)": ["K00376"],
    "nirK or nirS (NO-forming nitrite reductase)": ["K00368", "K15864"],
}


def conditions(sp, mad):
    """Condition name -> (family, kind, per-genome Series). Madin columns are mapped by species name."""
    m = lambda c: sp.map(mad[c])  # noqa: E731
    num = lambda c: pd.to_numeric(m(c), errors="coerce")  # noqa: E731
    ox = m("metabolism").map(lambda v: jc.OXYGEN.get(str(v).strip().lower(), np.nan))
    rt = m("range_tmp").astype("string").str.lower()
    gs = m("gram_stain").astype("string").str.lower()
    spo = m("sporulation").astype("string").str.lower()
    mot = m("motility").astype("string").str.lower()
    return {
        "oxygen: uses O2 (aerobe/facultative/microaerophile)":
            ("environment", "binary", ox.map({"aerobe": 1, "facultative": 1, "microaerophile": 1, "anaerobe": 0})),
        "oxygen: strict anaerobe": ("environment", "binary", ox.map(lambda v: np.nan if pd.isna(v) else
                                                                     float(v == "anaerobe"))),
        "growth temperature (C)": ("environment", "numeric", num("growth_tmp")),
        "optimum temperature (C)": ("environment", "numeric", num("optimum_tmp")),
        "thermophile (thermophilic / extreme)": ("environment", "binary",
                                                 rt.map(lambda v: np.nan if pd.isna(v) else
                                                        float("thermophil" in v))),
        "optimum pH": ("environment", "numeric", num("optimum_ph")),
        "doubling time (h, log10)": ("environment", "numeric", np.log10(num("doubling_h"))),
        "Gram-positive": ("cell trait", "binary", gs.map(lambda v: np.nan if pd.isna(v) else float(v == "positive"))),
        "spore-forming": ("cell trait", "binary", spo.map(lambda v: np.nan if pd.isna(v) else float(v == "yes"))),
        "motile": ("cell trait", "binary", mot.map(lambda v: np.nan if pd.isna(v) else float(v in ("yes", "flagella")))),
    }


def perm_stat(x, y, kind, rng, n=N_PERM):
    """Observed statistic and two-sided permutation p, trait shuffled across genera (vectorised)."""
    if kind == "numeric":
        ry = rankdata(y)
        ry = (ry - ry.mean()) / ry.std()
        rx = rankdata(x)
        rx = (rx - rx.mean()) / rx.std()
        obs = float(np.mean(rx * ry))
        P = np.array([rng.permutation(ry) for _ in range(n)])
        null = (P * rx).mean(axis=1)
    else:
        xs = (x - x.mean()) / x.std()
        ys = (y - y.mean()) / y.std()
        obs = float(np.mean(xs * ys))
        P = np.array([rng.permutation(ys) for _ in range(n)])
        null = (P * xs).mean(axis=1)
    p = (np.sum(np.abs(null) >= abs(obs) - 1e-12) + 1) / (n + 1)
    return obs, float(p)


def perm_within(x, y, kind, strata, rng, n=5000):
    """Two-sided p with the trait permuted only within strata (phyla). Strata with one genus stay fixed."""
    if kind == "numeric":
        yy = rankdata(y)
        xx = rankdata(x)
    else:
        yy, xx = y.astype(float), x.astype(float)
    xx = (xx - xx.mean()) / xx.std()
    yy = (yy - yy.mean()) / yy.std()
    obs = float(np.mean(xx * yy))
    idx = [np.where(strata == s_)[0] for s_ in np.unique(strata)]
    idx = [i for i in idx if len(i) > 1]
    if not idx:
        return obs, np.nan
    count = 0
    yp = yy.copy()
    for _ in range(n):
        for i in idx:
            yp[i] = yy[rng.permutation(i)]
        if abs(np.mean(xx * yp)) >= abs(obs) - 1e-12:
            count += 1
    return obs, (count + 1) / (n + 1)


def main():
    rng = np.random.default_rng(SEED)
    X, groups, role, ranks, info = jc.load_html(ROOT / "jaccard_similarity.html")
    G = pd.DataFrame({name: X[kos].max(axis=1) for name, kos in groups.items()})
    R = pd.DataFrame({name: X[[k for k in kos if k in X.columns]].max(axis=1) for name, kos in REFINED.items()})
    marker_set = {**{c: "group as defined" for c in G.columns}, **{c: "refined marker" for c in R.columns}}
    role = {**role, **{c: "refined" for c in R.columns}}
    G = pd.concat([G, R], axis=1)
    mad = pd.read_csv(ROOT / "data/raw/madin/condensed_species_NCBI.csv", low_memory=False) \
            .drop_duplicates("species").set_index("species")
    sp = info["species"].str.replace("Candidatus ", "", regex=False).str.split().str[:2].str.join(" ")
    genus = ranks["genus"].fillna("?")
    phylum_of_genus = pd.DataFrame({"genus": genus, "phylum": ranks["phylum"].fillna("?")}) \
        .groupby("genus").phylum.agg(lambda v: v.value_counts().index[0])

    rows = []
    for cname, (family, kind, series) in conditions(sp, mad).items():
        have = series.notna()
        if have.sum() == 0:
            continue
        df = pd.DataFrame({"genus": genus[have], "y": series[have].astype(float)})
        for gname in G.columns:
            df["x"] = G.loc[have, gname].astype(float)
            gl = df.groupby("genus").agg(x=("x", "mean"), y=("y", "mean"), n=("x", "size"))
            gl["x"] = (gl.x >= 0.5).astype(float)
            if kind == "binary":
                gl["y"] = (gl.y >= 0.5).astype(float)
            n_pos, n_neg = int(gl.x.sum()), int((gl.x == 0).sum())
            if len(gl) < MIN_GENERA or min(n_pos, n_neg) < MIN_MINORITY:
                continue
            if kind == "binary" and min(int(gl.y.sum()), int((gl.y == 0).sum())) < MIN_MINORITY:
                continue
            stat, p = perm_stat(gl.x.to_numpy(), gl.y.to_numpy(), kind, rng)
            strata = phylum_of_genus.reindex(gl.index).fillna("?").to_numpy()
            _, p_phy = perm_within(gl.x.to_numpy(), gl.y.to_numpy(), kind, strata, rng)
            row = dict(family=family, condition=cname, kind=kind, functional_group=gname,
                       marker_set=marker_set[gname], role=role[gname],
                       n_genomes=int(have.sum()), n_genera=len(gl), n_genera_marker_pos=n_pos,
                       n_genera_marker_neg=n_neg, statistic=stat,
                       statistic_name="Spearman rho" if kind == "numeric" else "phi (MCC)", p_permutation=p,
                       p_within_phylum=p_phy, n_phyla=int(len(set(strata))))
            if kind == "numeric":
                row["median_marker_pos"] = float(gl.y[gl.x == 1].median())
                row["median_marker_neg"] = float(gl.y[gl.x == 0].median())
            else:
                row["trait_rate_marker_pos"] = float(gl.y[gl.x == 1].mean())
                row["trait_rate_marker_neg"] = float(gl.y[gl.x == 0].mean())
            rows.append(row)
    out = pd.DataFrame(rows)
    out["q_bh"] = np.nan
    out["q_bh_within_phylum"] = np.nan
    for fam, ms in out[["family", "marker_set"]].drop_duplicates().itertuples(index=False):
        m = (out.family == fam) & (out.marker_set == ms)
        out.loc[m, "q_bh"] = jc.bh(out.loc[m, "p_permutation"].to_numpy())
        out.loc[m, "q_bh_within_phylum"] = jc.bh(out.loc[m, "p_within_phylum"].to_numpy())
    out["verdict"] = np.select(
        [(out.q_bh < 0.05) & (out.q_bh_within_phylum < 0.05), out.q_bh < 0.05],
        ["functional link (survives within-phylum test)",
         "lineage-bound (significant, but not separable from phylum identity in this data)"],
        "not significant")
    out["abs_stat"] = out.statistic.abs()
    out = out.sort_values(["family", "q_bh", "abs_stat"], ascending=[False, True, False]).drop(columns="abs_stat")
    out.to_csv(OUT / "environment_scan.csv", index=False, float_format="%.4f")

    for fam, ms in (("environment", "group as defined"), ("environment", "refined marker"),
                    ("cell trait", "group as defined")):
        d = out[(out.family == fam) & (out.marker_set == ms)]
        real = d[d.verdict.str.startswith("functional")]
        conf = d[d.verdict.str.startswith("lineage")]
        print(f"\n=== {fam} / {ms}: {len(d)} tests | {len(real)} functional links | {len(conf)} lineage-bound ===")
        cols = ["condition", "functional_group", "n_genera", "statistic", "median_marker_pos", "median_marker_neg",
                "trait_rate_marker_pos", "trait_rate_marker_neg", "q_bh", "q_bh_within_phylum"]
        print("-- functional links:")
        print(real[[c for c in cols if c in real.columns]].round(3).to_string(index=False))
        print("-- lineage-bound:")
        print(conf[[c for c in cols if c in conf.columns]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
