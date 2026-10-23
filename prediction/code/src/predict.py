"""Predict electron donors, electron acceptors and energy-metabolism pathways for a genome.

Input is the genome and nothing else: no taxonomy, no organism name, no literature. Two ways in:
  * an NCBI accession (all RefSeq GCF_ are annotated; GCA_ only if submitter-annotated) - proteins and coordinates
    are taken from the assembly's annotation;
  * a nucleotide FASTA of an unannotated genome (`--fasta file.fna [--fasta ...]`) - genes are called with Pyrodigal
    (Prodigal algorithm, meta mode for short inputs) and proteins are ordered by contig and start coordinate exactly
    as the annotated path does.

usage:
  python src/predict.py GCF_000008665.1 [GCF_... ...] [--arm released|A|B1r|B2|B2XL|B2E|C|D|E|E_MAG] [--out preds.csv]
  python src/predict.py --file accessions.txt
  python src/predict.py --fasta my_genome.fna [--fasta other.fna]

Three outputs from the one input, all derived from the predicted proteome and nothing else:
  <out>.csv                      electron-donor and electron-acceptor categories (trained model)
  <out>_pathways_predicted.csv   energy-metabolism pathways (trained model, results/models/pathways/)
  <out>_pathways.csv             KEGG module completeness, the deterministic annotation readout (evidence, not a
                                 prediction) -- `--pathways modules|model|both|none` selects which to write.
The pathway model scans the KO profiles its features need, so it is slower than the 852-KO donor/acceptor panel.

Pipeline (identical to training): NCBI Datasets download (proteins + GFF3) -> genome-ordered protein
table -> KOfam energy-panel annotation (KO arms) and/or ESM-2 / Bacformer embeddings (embedding arms)
-> saved per-trait models from results/models/<arm>/ -> probability + binary call at the stored
threshold. `in_training_set` flags genomes whose labels were used to fit the model (their predictions are
not independent evaluations); `genus_in_training_set` flags genomes from a genus (NCBI name) that contributed
training labels -- held-out-genus CV estimates are the relevant accuracy only when this flag is False.
Protein embeddings are computed in float16 on MPS/CUDA (as in training) and float32 on CPU (cosine similarity to
float16 > 0.99999 in our check).
New genomes are cached in data/predict_cache/ and never written into the training data directories.
"""
import argparse
import json
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
dl = __import__("01_download_genomes")
pt = __import__("02_build_protein_tables")
ko_mod = __import__("04_kofam_scan")
# XGBoost, PyTorch and pyhmmer each bring their own OpenMP runtime; in one macOS process some orders segfault and
# others deadlock. Each native stage therefore runs in a fresh spawned subprocess (ko_features_isolated,
# embedding_features_isolated); the main process only loads the saved models and predicts.
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
import multiprocessing as mp  # noqa: E402


def _isolated(fn_name, *args, **kw):
    with ProcessPoolExecutor(max_workers=1, mp_context=mp.get_context("spawn")) as ex:
        return ex.submit(_call, fn_name, args, kw).result()


def _call(fn_name, args, kw):
    import sys as _sys
    _sys.path.insert(0, str(ROOT / "src"))
    mod = __import__("predict")
    return getattr(mod, fn_name)(*args, **kw)
from redox_common import PROC, RES, RULE_GROUPS, apply_calibrator  # noqa: E402

KO_ARMS, EMB_ARMS, FUSION_ARMS = {"A", "B1r", "B2", "B2XL", "RF"}, {"C", "D", "E", "E_MAG"}, {"B2E"}


CACHE = ROOT / "data/predict_cache"


def proteins_from_fasta(path):
    """Call genes on an unannotated nucleotide FASTA (Pyrodigal) and return the same table the annotated path builds."""
    import hashlib

    import pyrodigal
    seqs, name, buf = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        seqs.append((name, "".join(buf)))
    total = sum(len(s) for _, s in seqs)
    finder = pyrodigal.GeneFinder(meta=total < 20000)
    if total >= 20000:
        finder.train(*[s for _, s in seqs])
    rows = []
    for ci, (contig, seq) in enumerate(seqs):
        for g in finder.find_genes(seq):
            rows.append((ci, contig, g.begin, g.end, "+" if g.strand > 0 else "-", g.translate().rstrip("*")))
    rows.sort(key=lambda r: (r[0], r[2], r[3]))
    df = pd.DataFrame({"pos": range(len(rows)), "contig": [r[1] for r in rows],
                       "contig_idx": [r[0] for r in rows], "start": [r[2] for r in rows], "end": [r[3] for r in rows],
                       "strand": [r[4] for r in rows],
                       "locus": [f"{r[1]}_{i}" for i, r in enumerate(rows)],
                       "protein_id": [f"{r[1]}_{i}" for i, r in enumerate(rows)],
                       "seq_hash": [hashlib.md5(r[5].encode()).hexdigest() for r in rows],
                       "length": [len(r[5]) for r in rows]})
    return df, {h: r[5] for h, r in zip(df.seq_hash, rows)}


def prepare_fastas(paths):
    tables, seqs, organism = {}, {}, {}
    for f in paths:
        name = Path(f).stem
        df, sq = proteins_from_fasta(f)
        tables[name], organism[name] = df, f"(unannotated FASTA: {Path(f).name})"
        seqs.update(sq)
        print(f"{name}: {len(df)} genes called on {df.contig_idx.nunique()} contigs", flush=True)
    return tables, seqs, organism


def prepare_proteins(accs):
    tables, seqs, organism = {}, {}, {}
    for acc in accs:
        known = (PROC / "proteins" / f"{acc}.parquet").exists()
        dl.OUT = pt.RAW = (ROOT / "data/raw/genomes") if known else (CACHE / "genomes")
        pt.OUT = (PROC / "proteins") if known else (CACHE / "proteins")
        pt.OUT.mkdir(parents=True, exist_ok=True)
        status = dl.download_one(acc)
        if status[1] != "ok":
            raise SystemExit(f"{acc}: download failed or no NCBI annotation ({status[1]})")
        uniq = {}
        pt.build(acc, uniq)
        seqs.update(uniq)
        tables[acc] = pd.read_parquet(pt.OUT / f"{acc}.parquet")
        rep = dl.OUT / acc / "assembly_report.json"
        organism[acc] = json.loads(rep.read_text().splitlines()[0]).get("organism", {}).get("organismName", "") \
            if rep.exists() else ""
    return tables, seqs, organism


def ko_features(tables, seqs):
    _, panel = ko_mod.build_panel()
    hmms, thr, stype = ko_mod.load_hmms(panel)
    hashes = list(seqs)
    hits = ko_mod.scan(hmms, thr, stype, hashes, [seqs[h] for h in hashes], log_every=0)
    return ko_mod.presence_matrix(hits, {a: t.seq_hash for a, t in tables.items()}, panel.KO)


def ko_features_for(tables, seqs, kos):
    """KO presence matrix restricted to an explicit KO list (the released pathway model's feature columns).

    Same KofamScan rules and thresholds as the training annotation; only the set of profiles scanned differs.
    """
    import pandas as _pd
    ko_list = _pd.read_csv(ko_mod.KOFAM / "ko_list.gz", sep="\t")
    ko_list = ko_list[ko_list.knum.isin(set(kos))].copy()
    ko_list["threshold"] = _pd.to_numeric(ko_list.threshold, errors="coerce")
    panel = ko_list[["knum", "threshold", "score_type"]].rename(columns={"knum": "KO"})
    hmms, thr, stype = ko_mod.load_hmms(panel)
    hashes = list(seqs)
    hits = ko_mod.scan(hmms, thr, stype, hashes, [seqs[h] for h in hashes], log_every=0)
    X = ko_mod.presence_matrix(hits, {a: t.seq_hash for a, t in tables.items()}, panel.KO)
    return X.reindex(columns=list(kos)).fillna(0.0)


def pathway_model_predict(tables, seqs):
    """Predicted energy-metabolism pathways from the released pathway model (results/models/pathways/)."""
    import joblib as _joblib
    rel = json.loads((RES / "models/pathways/RELEASED.json").read_text())
    mdir = RES / "models/pathways" / rel["arm"]
    meta = json.loads((mdir / "metadata.json").read_text())
    bundle = _joblib.load(mdir / "model.joblib")
    kos = sorted({c for b in bundle.values() for c in b["model"].columns if c.startswith("K")})
    X = ko_features_for(tables, seqs, kos)
    if any(c.startswith("mod_") for b in bundle.values() for c in b["model"].columns):
        pm = __import__("29_pathway_modules")
        mods = pd.read_csv(PROC / "kegg_modules_energy.tsv", sep="\t")
        M = pd.DataFrame(index=X.index, columns=["mod_" + m for m in mods.module], dtype=float)
        for acc in X.index:
            kset = set(X.columns[X.loc[acc] > 0])
            for m in mods.itertuples():
                steps = pm.parse_steps(m.definition)
                M.at[acc, "mod_" + m.module] = float(np.mean([pm.step_present(st, kset) for st in steps])) \
                    if steps else 0.0
        X = X.join(M)
    rows = []
    for t, b in bundle.items():
        p = b["model"].predict_proba(X)
        for acc, pi in zip(X.index, p):
            rows.append(dict(genome_id=acc, pathway=t, score=float(pi),
                             calibrated_score=float(apply_calibrator(b.get("calibrator"), [pi])[0]),
                             threshold=float(b["threshold"]), predicted=int(pi >= b["threshold"]),
                             arm=rel["arm"], call_reliable=t in meta.get("reliable_traits", []),
                             auprc_family_cv=meta.get("auprc_family", {}).get(t)))
    return pd.DataFrame(rows).sort_values(["genome_id", "score"], ascending=[True, False])


def embedding_features(tables, seqs, need_bags):
    import torch
    emb_mod = __import__("03_embed_esm2")
    ge = __import__("07_genome_embeddings")
    dev = "mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu")
    hashes = list(seqs)
    E = emb_mod.embed_sequences([seqs[h] for h in hashes], dev)
    index = pd.Index(hashes)
    bf, bf_mag = ge.load_bacformer(dev), ge.load_bacformer(dev, mag=True)
    frames, bags = [], {}
    for acc, prot in tables.items():
        e = E[index.get_indexer(prot.seq_hash)].astype(np.float32)
        reps, _ = ge.genome_representations(prot, e, bf, dev, bf_mag)
        frames.append(pd.DataFrame({f"{name}_{i}": [v] for name, vec in reps.items() for i, v in enumerate(vec)},
                                   index=[acc]))
        if need_bags:
            bags[acc] = e
    return pd.concat(frames), bags


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("accessions", nargs="*")
    ap.add_argument("--file")
    ap.add_argument("--fasta", action="append", default=[], help="unannotated nucleotide FASTA; genes are called")
    ap.add_argument("--arm", default="released")
    ap.add_argument("--out", default="predictions.csv")
    ap.add_argument("--embeddings", help="parquet of precomputed Bacformer-large genome vectors (rows = genome ids), "
                                         "needed by arm B2XLE_large; computed on a GPU by the kaggle/ embedding kernels")
    ap.add_argument("--pathways", default="auto", choices=["auto", "model", "modules", "both", "none"],
                    help="energy-metabolism pathway output: trained model, KEGG module completeness readout, both, "
                         "or none. 'auto' = model if one is released, else modules")
    ap.add_argument("--no-pathways", action="store_true", help=argparse.SUPPRESS)  # kept: earlier CLI
    args = ap.parse_args()
    accs = list(args.accessions)
    if args.file:
        accs += [line.strip() for line in open(args.file) if line.strip()]
    if not accs and not args.fasta:
        ap.error("give at least one accession or --fasta")
    arm = json.loads((RES / "models/RELEASED.json").read_text())["released_arm"] if args.arm == "released" \
        else args.arm
    mdir = RES / "models" / arm
    meta = json.loads((mdir / "metadata.json").read_text())
    traits = meta["traits"]

    tables, seqs, organism = prepare_proteins(accs) if accs else ({}, {}, {})
    if args.fasta:
        ft, fs, fo = prepare_fastas(args.fasta)
        tables.update(ft)
        seqs.update(fs)
        organism.update(fo)
    bundle = joblib.load(mdir / "model.joblib") if (mdir / "model.joblib").exists() else None
    probs, thrs = {}, {}
    if arm in KO_ARMS:
        X = _isolated("ko_features", tables, seqs)
        for t in traits:
            b = bundle[t]
            if b["type"] == "rule":
                groups = RULE_GROUPS[t]
                score = np.max([X[k].mean(axis=1).to_numpy() for k in groups.values()], axis=0) if groups \
                    else np.zeros(len(X))
                probs[t] = score
            else:
                probs[t] = b["model"].predict_proba(X)
            thrs[t] = b["threshold"]
        index = X.index
    elif arm == "B2XLE_large":  # B2XL + Bacformer-large: KO panel (computed here) + precomputed genome vectors
        if not args.embeddings:
            raise SystemExit("arm B2XLE_large needs --embeddings <genome_features_bacformer_large_mean.parquet>: "
                             "ESM-C + Bacformer-large vectors require a GPU and are computed by a kaggle/ kernel")
        emb = pd.read_parquet(args.embeddings)
        X = _isolated("ko_features", tables, seqs).join(emb, how="inner")
        missing = sorted(set(tables) - set(X.index))
        if missing:
            raise SystemExit(f"no Bacformer-large vector for {missing}")
        index = X.index
        for t in traits:
            probs[t], thrs[t] = np.asarray(bundle[t]["model"].predict_proba(X)), bundle[t]["threshold"]
    elif arm in FUSION_ARMS:  # B2E: energy-panel KO model + Bacformer-base probe (same features as training)
        X_emb = _isolated("embedding_features", tables, seqs, need_bags=False)[0]
        X = _isolated("ko_features", tables, seqs).join(X_emb, how="inner")
        index = X.index
        for t in traits:
            probs[t], thrs[t] = bundle[t]["model"].predict_proba(X), bundle[t]["threshold"]
    elif arm in EMB_ARMS:
        X, bags = _isolated("embedding_features", tables, seqs, need_bags=(arm == "D"))
        hits_for_bags = _isolated("ko_mod_hits_for_bags", tables, seqs) if arm == "D" else None
        index = X.index
        if arm == "D":
            import torch
            d_mod = __import__("06b_arm_D_abmil")
            ck = torch.load(mdir / "model.pt", weights_only=False)
            if ck["config"]["prefilter"]:
                hits = hits_for_bags
                bag_list = [bags[a][tables[a].seq_hash.isin(hits).to_numpy()] if tables[a].seq_hash.isin(hits).any()
                            else bags[a][:1] for a in index]
            else:
                bag_list = [bags[a] for a in index]
            P = []
            for sd in ck["state_dicts"]:
                m = d_mod.ABMIL(bag_list[0].shape[1], len(traits))
                m.load_state_dict(sd)
                P.append(d_mod.predict(m.eval(), bag_list, np.arange(len(bag_list)))[0])
            P = np.mean(P, axis=0)
            for j, t in enumerate(traits):
                probs[t], thrs[t] = P[:, j], ck["thresholds"][j]
        else:
            for t in traits:
                probs[t], thrs[t] = bundle[t]["model"].predict_proba(X), bundle[t]["threshold"]
    else:
        raise SystemExit(f"unknown arm {arm}")

    calibrators = {t: (bundle[t].get("calibrator") if bundle and t in bundle else None) for t in traits}
    labels = pd.read_csv(PROC / "labels.tsv", sep="\t").set_index("genome_id")
    train_genera = set(labels.loc[labels.labelled, "genus"])
    rows = []
    for t in traits:
        for acc, p in zip(index, probs[t]):
            in_train = bool(acc in labels.index and labels.at[acc, "labelled"])
            genus = labels.at[acc, "genus"] if acc in labels.index else (organism[acc].replace("Candidatus ", "")
                                                                         .split() or [""])[0]
            cal = calibrators.get(t)
            rows.append(dict(genome_id=acc, organism=organism[acc], trait=t, score=float(p),
                             calibrated_score=float(apply_calibrator(cal, [p])[0]),
                             threshold=float(thrs[t]), predicted=int(p >= thrs[t]), arm=arm,
                             in_training_set=in_train, genus_in_training_set=genus in train_genera,
                             training_label=(labels.at[acc, t] if in_train else np.nan)))
    out = pd.DataFrame(rows)
    out.to_csv(args.out, index=False)
    mode = "none" if args.no_pathways else args.pathways
    if mode == "auto":
        mode = "model" if (RES / "models/pathways/RELEASED.json").exists() else "modules"
    stem = str(Path(args.out).with_suffix(""))
    if mode in ("model", "both"):
        pw = _isolated("pathway_model_predict", tables, seqs)
        pw.to_csv(stem + "_pathways_predicted.csv", index=False, float_format="%.4f")
        print(f"written {stem}_pathways_predicted.csv ({pw.pathway.nunique()} pathways x "
              f"{pw.genome_id.nunique()} genomes)")
        print(pw[pw.predicted == 1].pivot_table(index="pathway", columns="genome_id", values="score",
                                                aggfunc="first").round(2).to_string())
    if mode in ("modules", "both"):
        pw = _isolated("pathway_readout", tables, seqs)
        pw.to_csv(stem + "_pathways.csv", index=False, float_format="%.3f")
        print(f"written {stem}_pathways.csv ({len(pw)} module calls)")
    print(out.pivot(index="genome_id", columns="trait", values="predicted").T.to_string())
    print(f"written {args.out}")


def pathway_readout(tables, seqs):
    """KEGG energy-metabolism module completeness per genome, from the same KO annotation (src/29_pathway_modules.py)."""
    pm = __import__("29_pathway_modules")
    X = ko_features(tables, seqs)
    mods = pd.read_csv(PROC / "kegg_modules_energy.tsv", sep="\t")
    rows = []
    for acc in X.index:
        kos = set(X.columns[X.loc[acc] > 0])
        for m in mods.itertuples():
            steps = pm.parse_steps(m.definition)
            if not steps:
                continue
            frac = float(np.mean([pm.step_present(st, kos) for st in steps]))
            if frac >= 0.5:
                rows.append(dict(genome_id=acc, module=m.module, module_name=m.name, completeness=frac,
                                 n_steps=len(steps), call="complete" if frac >= 0.8 else "partial"))
    return pd.DataFrame(rows).sort_values(["genome_id", "completeness"], ascending=[True, False])


def ko_mod_hits_for_bags(tables, seqs):
    """Hashes of proteins with a KOfam energy-panel hit (ABMIL bag prefilter)."""
    _, panel = ko_mod.build_panel()
    hmms, thr, stype = ko_mod.load_hmms(panel)
    hashes = list(seqs)
    return set(ko_mod.scan(hmms, thr, stype, hashes, [seqs[h] for h in hashes], log_every=0).seq_hash)


if __name__ == "__main__":
    main()
