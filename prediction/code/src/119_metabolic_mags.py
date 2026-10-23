"""Experiment M3 (docs/preregistration_M_metabolic_benchmarks.md): METABOLIC calls on MAGs and their paired isolates,
with the same Pyrodigal gene calls used for RedoxKO (src/88_real_mags.py, src/predict.py proteins_from_fasta), and
call agreement by MAG completeness compared with RedoxKO (results/additional/x2_pair_calls.csv).
Nucleotide FASTA are downloaded again from NCBI and deleted after gene calling; proteins are kept.
usage: python src/119_metabolic_mags.py proteins | scan | analyze
Outputs: data/additional/mags/faa/*.faa, results/additional/metabolic_bench/mag_rows.parquet, m3.json
"""
import importlib
import json
import re
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
m88 = importlib.import_module("88_real_mags")
m103 = importlib.import_module("103_metabolic_features")
from metabolic_calls import calls_from_rows  # noqa: E402

FAA = ROOT / "data/additional/mags/faa"
OUT = ROOT / "results/additional/metabolic_bench"
SEED, B = 20261017, 2000
pairs = pd.read_csv(ROOT / "results/additional/x2_pair_calls.csv")
ACCS = sorted(set(pairs.isolate) | set(pairs.mag))


def one_protein(acc):
    out = FAA / f"{acc}.faa"
    if out.exists() and out.stat().st_size > 0:
        return acc, "ok"
    a, st = m88.fetch_fna(acc)
    if st != "ok":
        return acc, st
    import predict
    fna = m88.FNA / f"{acc}.fna"
    _, seqs = predict.proteins_from_fasta(fna)
    tmp = out.with_suffix(".tmp")
    with open(tmp, "w") as f:
        for h, s in seqs.items():
            f.write(f">{h}\n{s.rstrip('*')}\n")
    tmp.rename(out)
    fna.unlink()
    return acc, "ok"


def proteins():
    FAA.mkdir(parents=True, exist_ok=True)
    m88.FNA.mkdir(parents=True, exist_ok=True)
    done = 0
    with ThreadPoolExecutor(4) as ex:
        res = []
        for r in ex.map(one_protein, ACCS):
            res.append(r)
            done += 1
            if done % 25 == 0:
                print(f"proteins {done}/{len(ACCS)}", flush=True)
    st = pd.DataFrame(res, columns=["accession", "status"])
    st.to_csv(OUT / "mag_protein_status.tsv", sep="\t", index=False)
    print(st.status.value_counts(), flush=True)


def scan():
    t = m103.table()
    names = sorted({h.strip() for hm in t["Hmm file"].astype(str) for h in re.split(r",\s*", hm)
                    if h.strip() and h != "nan"})
    todo = [(a, FAA / f"{a}.faa", names) for a in ACCS if (FAA / f"{a}.faa").exists()]
    rows = {}
    with ProcessPoolExecutor(10) as ex:
        for i, (a, r) in enumerate(ex.map(m103.scan_genome, todo, chunksize=2), 1):
            rows[a] = r
            if i % 25 == 0:
                print(f"scan {i}/{len(todo)}", flush=True)
    df = pd.DataFrame.from_dict(rows, orient="index").astype(np.float32)
    df.columns = [f"metabolic_{c}" for c in df.columns]
    df.to_parquet(OUT / "mag_rows.parquet")
    print("scan done", df.shape, flush=True)


def analyze():
    C = calls_from_rows(pd.read_parquet(OUT / "mag_rows.parquet"))
    cats = list(C.columns)
    d = pairs[pairs.trait.isin(cats) & pairs.isolate.isin(C.index) & pairs.mag.isin(C.index)].copy()
    d["m_iso"] = [C.at[a, t] for a, t in zip(d.isolate, d.trait)]
    d["m_mag"] = [C.at[a, t] for a, t in zip(d.mag, d.trait)]
    d["agree_R"] = (d.c_iso == d.c_mag).astype(int)
    d["agree_M"] = (d.m_iso == d.m_mag).astype(int)
    d["lost_R"] = ((d.c_iso == 1) & (d.c_mag == 0)).astype(int)
    d["lost_M"] = ((d.m_iso == 1) & (d.m_mag == 0)).astype(int)
    res = {"n_pairs": int(d[["isolate", "mag"]].drop_duplicates().shape[0]), "categories": cats, "bins": {}}
    for b, g in d.groupby("bin"):
        groups = {k: x for k, x in g.groupby("isolate")}
        keys = np.array(list(groups))
        rng = np.random.default_rng(SEED)

        def f(s):
            pr, pm = s[s.c_iso == 1], s[s.m_iso == 1]
            return [s.agree_R.mean() - s.agree_M.mean(),
                    (pr.lost_R.mean() if len(pr) else np.nan) - (pm.lost_M.mean() if len(pm) else np.nan)]
        bs = np.array([f(pd.concat([groups[k] for k in rng.choice(keys, len(keys))])) for _ in range(B)])
        pr, pm = g[g.c_iso == 1], g[g.m_iso == 1]
        res["bins"][b] = {"pairs": int(g[["isolate", "mag"]].drop_duplicates().shape[0]),
                          "agree_RedoxKO": float(g.agree_R.mean()), "agree_METABOLIC": float(g.agree_M.mean()),
                          "d_agree": float(g.agree_R.mean() - g.agree_M.mean()),
                          "d_agree_ci": [float(x) for x in np.nanpercentile(bs[:, 0], [2.5, 97.5])],
                          "lost_pos_RedoxKO": float(pr.lost_R.mean()), "lost_pos_METABOLIC": float(pm.lost_M.mean()),
                          "d_lost_ci": [float(x) for x in np.nanpercentile(bs[:, 1], [2.5, 97.5])]}
    lo = res["bins"]["50-70"]["d_agree_ci"][0]
    res["M3_primary_success"] = bool(lo > 0)
    (OUT / "m3.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    {"proteins": proteins, "scan": scan, "analyze": analyze}[sys.argv[1]]()
