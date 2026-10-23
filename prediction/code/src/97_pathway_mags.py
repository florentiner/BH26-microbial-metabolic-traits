"""X2 extension: the released pathway model (PB3) on the real MAG-isolate pairs, both predicted from nucleotide FASTA
through src/predict.py (Pyrodigal; full 10 517-KO profile with the training KOfam release 2026-08-25).

usage: python src/97_pathway_mags.py predict | analyze
"""
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
FNA = ROOT / "data/additional/mags/fna"
PRED = ROOT / "data/additional/mags/pred_pathways"
OUT = ROOT / "results/additional"
SEED = 20261014


def predict(batch=10):
    PRED.mkdir(parents=True, exist_ok=True)
    st = pd.read_csv(ROOT / "data/additional/mags/download_status.tsv", sep="\t")
    accs = sorted(st[st.status == "ok"].accession)
    done = {p.stem[2:] for p in PRED.glob("b_*.csv")}
    todo = [a for a in accs if a not in done]
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        out = PRED / f"batch_{i:05d}_{chunk[0]}.csv"
        cmd = [sys.executable, str(ROOT / "src/predict.py"), "--arm", "B2", "--pathways", "model", "--out", str(out)]
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
    pairs = pd.read_csv(ROOT / "data/additional/mags/pairs.tsv", sep="\t")
    P = pd.concat([pd.read_csv(p) for p in sorted(PRED.glob("batch_*_pathways_predicted.csv"))], ignore_index=True)
    g = "genome_id" if "genome_id" in P else P.columns[0]
    sc = P.pivot_table(index=g, columns="pathway", values="score")
    cl = P.pivot_table(index=g, columns="pathway", values="predicted")
    pairs = pairs[pairs.mag.isin(sc.index) & pairs.isolate.isin(sc.index)]
    rows = []
    for _, r in pairs.iterrows():
        for t in sc.columns:
            rows.append(dict(isolate=r.isolate, mag=r.mag, bin=r.bin, trait=t, s_iso=sc.at[r.isolate, t],
                             s_mag=sc.at[r.mag, t], c_iso=int(cl.at[r.isolate, t]), c_mag=int(cl.at[r.mag, t])))
    D = pd.DataFrame(rows)
    D["agree"] = (D.c_iso == D.c_mag).astype(int)
    D.to_csv(OUT / "x2_pathway_pair_calls.csv", index=False)
    rng = np.random.default_rng(SEED)
    summ = []
    for b, d in D.groupby("bin"):
        sp = d.isolate.unique()
        idx = {s: np.where(d.isolate.to_numpy() == s)[0] for s in sp}
        v = [d.agree.to_numpy()[np.concatenate([idx[s] for s in rng.choice(sp, len(sp))])].mean()
             for _ in range(1000)]
        lo, hi = np.percentile(v, [2.5, 97.5])
        summ.append(dict(bin=b, mags=d.mag.nunique(), agreement=d.agree.mean(), agreement_ci=f"[{lo:.3f}, {hi:.3f}]",
                         mean_abs_score_diff=(d.s_iso - d.s_mag).abs().mean()))
    S = pd.DataFrame(summ)
    S.to_csv(OUT / "x2_pathway_agreement_by_bin.csv", index=False)
    print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    {"predict": predict, "analyze": analyze}[sys.argv[1]]()
