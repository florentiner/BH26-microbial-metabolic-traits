"""Run the released pathway model (PB3) on the prospective genomes, from the local cache only.

src/predict.py re-downloads every proteome through NCBI before annotating, which hung when the network came back
and is pure waste here: all 167 proteomes are already in data/predict_cache. This script reads the cache, scans
only the KO profiles the released model needs, uses every core, and logs progress.

usage: python src/79_pathway_validation_predict.py <accessions.txt> <out.csv> [--cpus 11]
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
CACHE = ROOT / "data/predict_cache"


def proteins(accs):
    """seq_hash lists per genome and the unique sequences, from the local cache (no network)."""
    import importlib
    pt = importlib.import_module("02_build_protein_tables")
    pt.RAW, pt.OUT = CACHE / "genomes", CACHE / "proteins"
    tables, uniq = {}, {}
    for a in accs:
        pt.build(a, uniq)
        tables[a] = pd.read_parquet(pt.OUT / f"{a}.parquet")
    return tables, uniq


def main():
    import importlib
    import joblib
    accs = [a.strip() for a in open(sys.argv[1]) if a.strip()]
    out_csv = Path(sys.argv[2])
    cpus = int(sys.argv[sys.argv.index("--cpus") + 1]) if "--cpus" in sys.argv else 11
    import json
    rel = json.loads((ROOT / "results/models/pathways/RELEASED.json").read_text())
    bundle = joblib.load(ROOT / f"results/models/pathways/{rel['arm']}/model.joblib")
    cols = {c for v in bundle.values() for c in v["model"].columns}
    kos = sorted(c for c in cols if str(c).startswith("K"))
    mods = sorted(c for c in cols if str(c).startswith("mod_"))
    print(f"{rel['arm']}: {len(bundle)} traits, {len(kos)} KO profiles, {len(mods)} module features", flush=True)

    t0 = time.time()
    tables, uniq = proteins(accs)
    print(f"{len(tables)} genomes, {len(uniq)} unique proteins from cache ({time.time() - t0:.0f}s)", flush=True)

    ko_mod = importlib.import_module("04_kofam_scan")
    kl = pd.read_csv(ko_mod.KOFAM / "ko_list.gz", sep="\t")
    kl = kl[kl.knum.isin(set(kos))].copy()
    kl["threshold"] = pd.to_numeric(kl.threshold, errors="coerce")
    panel = kl[["knum", "threshold", "score_type"]].rename(columns={"knum": "KO"})
    hmms, thr, stype = ko_mod.load_hmms(panel)
    print(f"loaded {len(hmms)} HMM profiles ({time.time() - t0:.0f}s); scanning on {cpus} cores", flush=True)
    hashes = list(uniq)
    hits = ko_mod.scan(hmms, thr, stype, hashes, [uniq[h] for h in hashes], cpus=cpus, log_every=200)
    X = ko_mod.presence_matrix(hits, {a: t.seq_hash for a, t in tables.items()}, panel.KO)
    X = X.reindex(columns=kos).fillna(0.0)
    print(f"KO matrix {X.shape} ({time.time() - t0:.0f}s)", flush=True)
    X.to_parquet(out_csv.with_name(out_csv.stem + "_ko_matrix.parquet"))   # keep the scan: every pathway arm can reuse it

    if mods:
        pm = importlib.import_module("29_pathway_modules")
        md = pd.read_csv(ROOT / "data/processed/kegg_modules_energy.tsv", sep="\t")
        M = pd.DataFrame(index=X.index, columns=mods, dtype=float)
        for acc in X.index:
            kset = set(X.columns[X.loc[acc] > 0])
            for m, defn in zip(md.module, md.definition):
                if f"mod_{m}" in M.columns:
                    M.loc[acc, f"mod_{m}"] = float(np.mean([pm.step_present(st, kset) for st in pm.parse_steps(defn)]))
        X = X.join(M.fillna(0.0))
    rows = []
    for trait, b in bundle.items():
        p = np.asarray(b["model"].predict_proba(X.reindex(columns=b["model"].columns).fillna(0.0)))
        p = p[:, 1] if p.ndim > 1 and p.shape[1] > 1 else np.ravel(p)
        thr_t = float(b["threshold"])
        for acc, s in zip(X.index, p):
            rows.append(dict(genome_id=acc, trait=trait, arm=rel["arm"], score=float(s), threshold=thr_t,
                             predicted=int(s >= thr_t)))
    d = pd.DataFrame(rows)
    d.to_csv(out_csv, index=False)
    print(f"written {out_csv}: {d.genome_id.nunique()} genomes x {d.trait.nunique()} traits "
          f"({time.time() - t0:.0f}s total)")


if __name__ == "__main__":
    main()
