"""X4 (docs/preregistration_additional_experiments.md): habitat case study on uncultured lineages.

GTDB R232 species representatives that are MAGs (CheckM2 completeness >= 90, contamination <= 5) whose NCBI
isolation source matches exactly one habitat class (keyword lists fixed in the pre-registration), at most one per
GTDB genus and 60 per class (seed 20261014). B2 predictions from nucleotide FASTA via src/predict.py; the
pre-registered expectations X4-E1..E5 are tested with one-sided Fisher / Mann-Whitney tests and BH correction.

usage: python src/90_habitat_case.py select|download|predict|analyze
"""
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
mags_mod = __import__("88_real_mags")
WORK = ROOT / "data/additional/habitat"
FNA = ROOT / "data/additional/mags/fna"  # shared FASTA cache
PRED = WORK / "pred"
OUT = ROOT / "results/additional"
SEED = 20261014
KEYS = {
    "hydrothermal": ["hydrothermal", "vent", "black smoker", "chimney"],
    "acid_mine": ["acid mine", "mine drainage", r"\bAMD\b"],
    "human_gut": ["human gut", "human feces", "human faeces", "human stool"],
    "digester": ["anaerobic digester", "biogas", "digestor"],
    "marine_surface": ["seawater", "surface water", "ocean water"],
}


def classify(src):
    if not isinstance(src, str) or not src.strip() or src.lower() in ("none", "missing", "not applicable"):
        return None
    hits = [c for c, ks in KEYS.items() if any(re.search(k if k.startswith(r"\b") else re.escape(k), src,
                                                         flags=0 if k.startswith(r"\b") else re.I) for k in ks)]
    if "marine_surface" in hits and "hydrothermal" in hits:  # marine_surface excludes any hydrothermal keyword
        hits.remove("marine_surface")
    return hits[0] if len(hits) == 1 else None


def select():
    g = pd.read_parquet(ROOT / "data/additional/gtdb/gtdb232_slim.parquet")
    g = g[(g.gtdb_representative.astype(str).str.lower().isin(["t", "true"]))
          & (g.ncbi_genome_category == "derived from metagenome")
          & (g.checkm2_completeness >= 90) & (g.checkm2_contamination <= 5)].copy()
    g["habitat"] = g.ncbi_isolation_source.map(classify)
    g = g[g.habitat.notna()]
    g["genus"] = g.gtdb_taxonomy.str.split(";").str[-2]
    g["mag"] = g.ncbi_genbank_assembly_accession.where(g.ncbi_genbank_assembly_accession.notna(), g.accession.str[3:])
    rng = np.random.default_rng(SEED)
    rows = []
    for h in KEYS:
        d = g[g.habitat == h].sort_values("mag")
        d = d.iloc[rng.permutation(len(d))].drop_duplicates("genus")
        rows.append(d.head(60))
    sel = pd.concat(rows)
    WORK.mkdir(parents=True, exist_ok=True)
    sel[["mag", "habitat", "gtdb_taxonomy", "checkm2_completeness", "checkm2_contamination",
         "ncbi_isolation_source"]].to_csv(WORK / "selection.tsv", sep="\t", index=False)
    print(g.habitat.value_counts().to_dict(), "\nselected", sel.habitat.value_counts().to_dict())


def download():
    sel = pd.read_csv(WORK / "selection.tsv", sep="\t")
    FNA.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(4) as ex:
        res = list(ex.map(mags_mod.fetch_fna, sorted(sel.mag)))
    st = pd.DataFrame(res, columns=["accession", "status"])
    st.to_csv(WORK / "download_status.tsv", sep="\t", index=False)
    print(st.status.str.split(":").str[0].value_counts())


def predict(batch=20):
    PRED.mkdir(parents=True, exist_ok=True)
    st = pd.read_csv(WORK / "download_status.tsv", sep="\t")
    accs = sorted(st[st.status == "ok"].accession)
    done = {p.stem[2:] for p in PRED.glob("b_*.csv")}
    todo = [a for a in accs if a not in done]
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        out = PRED / f"batch_{i:05d}_{chunk[0]}.csv"
        cmd = [sys.executable, str(ROOT / "src/predict.py"), "--arm", "B2", "--pathways", "none", "--out", str(out)]
        for a in chunk:
            cmd += ["--fasta", str(FNA / f"{a}.fna")]
        t0 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print("FAILED", chunk[0], r.stderr[-800:], flush=True)
            continue
        for a in chunk:
            (PRED / f"b_{a}.csv").write_text("done\n")
        print(f"{i + len(chunk)}/{len(todo)} in {time.time() - t0:.0f}s", flush=True)


def analyze():
    sel = pd.read_csv(WORK / "selection.tsv", sep="\t").set_index("mag")
    P = pd.concat([pd.read_csv(p) for p in sorted(PRED.glob("batch_*.csv"))], ignore_index=True)
    gcol = "genome" if "genome" in P else P.columns[0]
    pcol = "score"
    sc = P.pivot_table(index=gcol, columns="trait", values=pcol)
    cl = P.pivot_table(index=gcol, columns="trait", values="predicted")
    sc["habitat"] = cl["habitat"] = sel.habitat.reindex(sc.index)
    frac = cl.groupby("habitat").mean(numeric_only=True).T
    frac.loc["n_genomes"] = cl.groupby("habitat").size()
    frac.to_csv(OUT / "x4_predicted_fraction_by_habitat.csv")

    def fisher(trait, a, b):  # one-sided: trait more frequent in a than in b
        x, y = cl[cl.habitat == a][trait], cl[cl.habitat == b][trait]
        t = [[int(x.sum()), int(len(x) - x.sum())], [int(y.sum()), int(len(y) - y.sum())]]
        return fisher_exact(t, alternative="greater")[1], x.mean(), y.mean()

    tests = []
    for name, trait, a, b in [("X4-E1", "donor_hydrogen_oxidation", "hydrothermal", "marine_surface"),
                              ("X4-E2", "donor_sulfur_oxidation", "hydrothermal", "human_gut"),
                              ("X4-E4", "acceptor_CO2_methanogen_acetogen", "digester", "marine_surface"),
                              ("X4-E5a", "acceptor_oxygen", "marine_surface", "human_gut"),
                              ("X4-E5b", "acceptor_oxygen", "marine_surface", "digester")]:
        p, fa, fb = fisher(trait, a, b)
        tests.append(dict(test=name, trait=trait, higher_in=a, lower_in=b, frac_a=fa, frac_b=fb, p=p,
                          kind="Fisher one-sided on released calls"))
    for other in [h for h in KEYS if h != "acid_mine"]:
        x = sc[sc.habitat == "acid_mine"]["donor_iron_oxidation"]
        y = sc[sc.habitat == other]["donor_iron_oxidation"]
        tests.append(dict(test=f"X4-E3 vs {other}", trait="donor_iron_oxidation", higher_in="acid_mine",
                          lower_in=other, frac_a=x.median(), frac_b=y.median(),
                          p=mannwhitneyu(x, y, alternative="greater").pvalue,
                          kind="Mann-Whitney one-sided on raw scores (medians shown)"))
    T = pd.DataFrame(tests)
    T["q_bh"] = multipletests(T.p, method="fdr_bh")[1]
    T["supported"] = T.q_bh < 0.05
    T.to_csv(OUT / "x4_tests.csv", index=False)
    pd.set_option("display.width", 220)
    print(frac.round(2).to_string())
    print(T.round(4).to_string(index=False))


if __name__ == "__main__":
    {"select": select, "download": download, "predict": predict, "analyze": analyze}[sys.argv[1]]()
