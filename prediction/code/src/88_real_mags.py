"""X2 (docs/preregistration_additional_experiments.md): predictions from real metagenome-assembled genomes versus the
isolate genome of the same GTDB R232 species.

Pairs: our evaluation isolates (external 562 set, 143 blind newly described species) whose GTDB species cluster
contains genomes "derived from metagenome" (CheckM2 completeness >= 50, contamination <= 10). Per species one MAG is
drawn at random from each completeness bin [50,70), [70,90), [90,100] that is populated (seed 20261014, so up to 3).
Both the MAG and the isolate are predicted from nucleotide FASTA through the same Pyrodigal path of src/predict.py
(deviation D3: the isolate FASTA path was added so that gene calling is identical on both sides). Model: B2.

usage: python src/88_real_mags.py select|download|predict|analyze
"""
import io
import json
import subprocess
import sys
import time
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "data/additional/mags"
FNA = WORK / "fna"
PRED = WORK / "pred"
OUT = ROOT / "results/additional"
SEED = 20261014
BINS = [50, 70, 90, 100.01]
BIN_NAMES = ["50-70", "70-90", "90-100"]
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA"


def base(a):
    return a.split("_", 1)[1].split(".")[0] if a[:3] in ("RS_", "GB_") else a.split("_", 1)[1].split(".")[0]


def select():
    g = pd.read_parquet(ROOT / "data/additional/gtdb/gtdb232_slim.parquet")
    g["base"] = g.accession.str[3:].map(base)
    g["species"] = g.gtdb_taxonomy.str.split(";").str[-1]
    ext = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t")
    new = pd.read_csv(ROOT / "validation_new_organisms/scored178_B2.csv")
    new = sorted(new[new.blind == True].accession.unique())  # noqa: E712
    iso = pd.DataFrame({"isolate": list(ext.genome_id) + new, "set": ["external"] * len(ext) + ["new"] * len(new)})
    iso["base"] = iso.isolate.map(base)
    iso = iso.merge(g[["base", "species"]], on="base", how="inner").drop_duplicates("isolate")
    mags = g[(g.ncbi_genome_category == "derived from metagenome") & (g.checkm2_completeness >= 50)
             & (g.checkm2_contamination <= 10)].copy()
    mags["bin"] = pd.cut(mags.checkm2_completeness, BINS, right=False, labels=BIN_NAMES).astype(str)
    mags["mag"] = mags.ncbi_genbank_assembly_accession.where(mags.ncbi_genbank_assembly_accession.notna(),
                                                             mags.accession.str[3:])
    rng = np.random.default_rng(SEED)
    rows = []
    for _, r in iso.sort_values("isolate").iterrows():
        cand = mags[mags.species == r.species]
        for b in BIN_NAMES:
            c = cand[cand.bin == b].sort_values("mag")
            if len(c):
                m = c.iloc[int(rng.integers(len(c)))]
                rows.append(dict(isolate=r.isolate, set=r.set, species=r.species, mag=m.mag, bin=b,
                                 completeness=m.checkm2_completeness, contamination=m.checkm2_contamination,
                                 isolation_source=m.ncbi_isolation_source))
    pairs = pd.DataFrame(rows)
    WORK.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(WORK / "pairs.tsv", sep="\t", index=False)
    print(pairs.groupby(["set", "bin"]).size(), "\nMAGs", len(pairs), "isolates", pairs.isolate.nunique())


def fetch_fna(acc):
    target = FNA / f"{acc}.fna"
    if target.exists() and target.stat().st_size > 0:
        return acc, "ok"
    for i in range(6):
        try:
            with urllib.request.urlopen(API.format(acc=acc), timeout=600) as r:
                z = zipfile.ZipFile(io.BytesIO(r.read()))
            fn = [n for n in z.namelist() if n.endswith((".fna", ".fa", ".fasta"))]
            if not fn:
                return acc, "no_fasta"
            tmp = target.with_suffix(".tmp")
            with open(tmp, "wb") as f:
                for n in fn:
                    f.write(z.read(n))
            tmp.rename(target)
            time.sleep(0.3)
            return acc, "ok"
        except Exception as e:  # noqa: BLE001
            time.sleep(10 * (i + 1))
            err = str(e)
    return acc, f"error:{err}"


def download():
    pairs = pd.read_csv(WORK / "pairs.tsv", sep="\t")
    FNA.mkdir(parents=True, exist_ok=True)
    accs = sorted(set(pairs.mag) | set(pairs.isolate))
    with ThreadPoolExecutor(4) as ex:
        res = list(ex.map(fetch_fna, accs))
    st = pd.DataFrame(res, columns=["accession", "status"])
    st.to_csv(WORK / "download_status.tsv", sep="\t", index=False)
    print(st.status.str.split(":").str[0].value_counts())


def predict(batch=20):
    PRED.mkdir(parents=True, exist_ok=True)
    st = pd.read_csv(WORK / "download_status.tsv", sep="\t")
    accs = sorted(st[st.status == "ok"].accession)
    done = {p.stem for p in PRED.glob("*.csv")}
    todo = [a for a in accs if f"b_{a}" not in done]
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


def load_predictions():
    frames = [pd.read_csv(p) for p in sorted(PRED.glob("batch_*.csv"))]
    return pd.concat(frames, ignore_index=True)


def analyze():
    pairs = pd.read_csv(WORK / "pairs.tsv", sep="\t")
    P = load_predictions()
    gcol = "genome" if "genome" in P else P.columns[0]
    P = P.rename(columns={gcol: "genome"})
    pcol = "score"
    traits = sorted(P.trait.unique())
    sc = P.pivot_table(index="genome", columns="trait", values=pcol)
    cl = P.pivot_table(index="genome", columns="trait", values="predicted")
    have = set(sc.index)
    pairs = pairs[pairs.mag.isin(have) & pairs.isolate.isin(have)].copy()
    rows = []
    for _, r in pairs.iterrows():
        for t in traits:
            rows.append(dict(isolate=r.isolate, mag=r.mag, set=r.set, bin=r.bin, completeness=r.completeness,
                             trait=t, s_iso=sc.at[r.isolate, t], s_mag=sc.at[r.mag, t],
                             c_iso=int(cl.at[r.isolate, t]), c_mag=int(cl.at[r.mag, t])))
    D = pd.DataFrame(rows)
    D["agree"] = (D.c_iso == D.c_mag).astype(int)
    D["mag_only_neg"] = ((D.c_iso == 1) & (D.c_mag == 0)).astype(int)
    D["mag_only_pos"] = ((D.c_iso == 0) & (D.c_mag == 1)).astype(int)
    D["abs_diff"] = (D.s_iso - D.s_mag).abs()
    D.to_csv(OUT / "x2_pair_calls.csv", index=False)

    def boot_ci(df, col, n=1000):
        rng = np.random.default_rng(SEED)
        sp = df.isolate.unique()
        idx = {s: np.where(df.isolate.to_numpy() == s)[0] for s in sp}
        v = [df[col].to_numpy()[np.concatenate([idx[s] for s in rng.choice(sp, len(sp))])].mean()
             for _ in range(n)]
        return np.percentile(v, [2.5, 97.5])

    summ = []
    for b, d in D.groupby("bin"):
        lo, hi = boot_ci(d, "agree")
        summ.append(dict(bin=b, mags=d.mag.nunique(), species=d.isolate.nunique(), agreement=d.agree.mean(),
                         agreement_ci=f"[{lo:.3f}, {hi:.3f}]", mean_abs_score_diff=d.abs_diff.mean(),
                         mag_only_negative=int(d.mag_only_neg.sum()), mag_only_positive=int(d.mag_only_pos.sum())))
    summ = pd.DataFrame(summ)
    summ.to_csv(OUT / "x2_agreement_by_bin.csv", index=False)
    per_trait = D.groupby(["bin", "trait"]).agg(agreement=("agree", "mean"), n=("agree", "size")).reset_index()
    per_trait.to_csv(OUT / "x2_agreement_per_trait.csv", index=False)

    # accuracy against labels where the isolate has them
    ext = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t").set_index("genome_id")
    new = pd.read_csv(ROOT / "validation_new_organisms/scored178_B2.csv")
    new = new[(new.blind == True) & ~new.category.isin(["donor_other", "acceptor_other"])]  # noqa: E712
    lab = []
    for _, r in D.iterrows():
        y = np.nan
        if r.set == "external" and r.trait in ext.columns:
            y = ext.at[r.isolate, r.trait]
        elif r.set == "new":
            m = new[(new.accession == r.isolate) & (new.category == r.trait)]
            if len(m):
                y = m.expected.iloc[0]
        lab.append(y)
    D["y"] = lab
    L = D[D.y.notna()].copy()
    acc = []
    for b, d in L.groupby("bin"):
        def macro_auc(df, col):
            v = [roc_auc_score(g.y, g[col]) for _, g in df.groupby("trait") if g.y.nunique() == 2]
            return float(np.mean(v)) if v else np.nan
        rng = np.random.default_rng(SEED)
        sp = d.isolate.unique()
        idx = {s: np.where(d.isolate.to_numpy() == s)[0] for s in sp}
        deltas = []
        for _ in range(1000):
            dd = d.iloc[np.concatenate([idx[s] for s in rng.choice(sp, len(sp))])]
            deltas.append(macro_auc(dd, "s_mag") - macro_auc(dd, "s_iso"))
        lo, hi = np.nanpercentile(deltas, [2.5, 97.5])
        acc.append(dict(bin=b, labelled_calls=len(d), species=d.isolate.nunique(),
                        auroc_isolate=macro_auc(d, "s_iso"), auroc_mag=macro_auc(d, "s_mag"),
                        auroc_mag_minus_isolate=macro_auc(d, "s_mag") - macro_auc(d, "s_iso"),
                        delta_ci=f"[{lo:+.3f}, {hi:+.3f}]",
                        mcc_isolate=matthews_corrcoef(d.y, d.c_iso), mcc_mag=matthews_corrcoef(d.y, d.c_mag)))
    acc = pd.DataFrame(acc)
    acc.to_csv(OUT / "x2_accuracy_by_bin.csv", index=False)
    pd.set_option("display.width", 200)
    print(summ.round(3).to_string(index=False))
    print(acc.round(3).to_string(index=False))


if __name__ == "__main__":
    {"select": select, "download": download, "predict": predict, "analyze": analyze}[sys.argv[1]]()
