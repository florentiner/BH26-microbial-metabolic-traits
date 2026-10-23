"""A (docs/preregistration_hybrid_AB.md): METABOLIC-G v4 per-HMM presence calls (its 315 function-table rows),
recomputed with pyhmmer so they can be produced for every genome:
  * custom HMMs (METABOLIC_hmm_db) with the template's cut-off and score type; KO HMMs with KOfam 2026-08-25 cut-offs
    (ko_list; '-' -> 50 full), exactly as METABOLIC-G.pl merges them;
  * 'full' -> sequence score >= cut-off; 'domain' -> best single-domain score >= cut-off;
  * motif check (regex, X = any residue) and motif-pair check (check.hmm score >= anti check.hmm score, > 0);
  * a row with several HMMs is present when any of them is present.
Proteins: RefSeq protein.faa (training, external, new species). Output: data/processed/genome_features_metabolic_hmm.parquet
(training) and data/additional/metabolic_hmm/{external,newspecies}.parquet; `validate` compares with real METABOLIC
output on the 100 primary genomes.

usage: python src/103_metabolic_features.py build | validate
"""
import gzip
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MET = ROOT / "data/additional/metabolic/METABOLIC"
TPL = MET / "METABOLIC_template_and_database"
KOFAM = ROOT / "data/raw/kofam"
OUTD = ROOT / "data/additional/metabolic_hmm"


def table():
    t = pd.read_csv(TPL / "hmm_table_template.txt", sep="\t")
    t["row"] = [f"{c}|{f}|{g}" for c, f, g in zip(t.Category, t.Function, t["Gene abbreviation"])]
    return t


def thresholds(t):
    thr = {}
    for hm, th in zip(t["Hmm file"].astype(str), t["Hmm detecting threshold"].astype(str)):
        for h, x in zip(re.split(r",\s*", hm), re.split(r",\s*", th)):
            if h and h != "nan" and "|" in x:
                v, s = x.split("|")
                thr[h.strip()] = (float(v), s.strip())
    ko = pd.read_csv(KOFAM / "ko_list.gz", sep="\t")
    for k, v, s in zip(ko.knum, ko.threshold, ko.score_type):
        h = f"{k}.hmm"
        if h in thr:
            thr[h] = (50.0, "full") if str(v) == "-" else (float(v), str(s))
    return thr


def motifs():
    m = dict(l.strip().split(":", 1) for l in open(TPL / "motif.txt") if ":" in l)
    pair = dict(l.strip().split(":", 1) for l in open(TPL / "motif.pair.txt") if ":" in l)
    rx = {k: re.compile(v.replace("X", "[ARNDCQEGHILKMFPSTWYV]")) for k, v in m.items()}
    return rx, pair


def load_hmm(path):
    import pyhmmer
    with pyhmmer.plan7.HMMFile(str(path)) as f:
        return list(f)


def hmm_path(h):
    p = MET / "METABOLIC_hmm_db" / h
    return p if p.exists() else KOFAM / "profiles" / h


def scan_genome(args):
    acc, faa, hmm_names = args
    import pyhmmer
    t = table()
    thr = thresholds(t)
    rx, pair = motifs()
    alph = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(str(faa), digital=True, alphabet=alph) as sf:
        seqs = sf.read_block()
    text = {s.name.decode() if isinstance(s.name, bytes) else s.name: s.textize().sequence for s in seqs}
    hmms = []
    for h in hmm_names:
        p = hmm_path(h)
        if p.exists():
            hmms += load_hmm(p)
    names = [h for h in hmm_names if hmm_path(h).exists()]
    present = set()
    check_cache = {}
    for name, top in zip(names, pyhmmer.hmmsearch(hmms, seqs, cpus=1)):
        v, st = thr[name]
        base = name[:-4]
        for hit in top:
            score = hit.score if st == "full" else hit.best_domain.score
            if score < v:
                continue
            pid = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            if base in rx:
                if not rx[base].search(text[pid]):
                    continue
            elif base in pair:
                def chk(b):
                    key = (b, pid)
                    if key not in check_cache:
                        ch = load_hmm(MET / "METABOLIC_hmm_db" / f"{b}.check.hmm")
                        one = pyhmmer.easel.TextSequence(name=pid.encode(), sequence=text[pid]).digitize(alph)
                        res = list(pyhmmer.hmmsearch(ch, pyhmmer.easel.DigitalSequenceBlock(alph, [one])))
                        check_cache[key] = max((h.score for h in res[0]), default=0.0)
                    return check_cache[key]
                s1, s2 = chk(base), chk(pair[base])
                if not (s1 >= s2 and s1 != 0):
                    continue
            present.add(name)
            break
    row = {}
    for r, hm in zip(t.row, t["Hmm file"].astype(str)):
        hs = [h.strip() for h in re.split(r",\s*", hm) if h.strip() and h != "nan"]
        row[r] = int(any(h in present for h in hs))
    return acc, row


def proteome(acc):
    for d in ["data/raw/genomes", "data/external/genomes", "data/predict_cache/genomes"]:
        p = ROOT / d / acc / "protein.faa"
        if p.exists():
            return p
    return None


def build_set(accs, out, jobs=8):
    t = table()
    names = sorted({h.strip() for hm in t["Hmm file"].astype(str) for h in re.split(r",\s*", hm)
                    if h.strip() and h != "nan"})
    todo = [(a, proteome(a), names) for a in accs if proteome(a) is not None]
    rows = {}
    with ProcessPoolExecutor(jobs) as ex:
        for i, (a, r) in enumerate(ex.map(scan_genome, todo, chunksize=2), 1):
            rows[a] = r
            if i % 50 == 0:
                print(out.name, i, "/", len(todo), flush=True)
    df = pd.DataFrame.from_dict(rows, orient="index").astype(np.float32)
    df.columns = [f"metabolic_{c}" for c in df.columns]
    df.index.name = "genome_id"
    df.to_parquet(out)
    print(out, df.shape)


def build():
    OUTD.mkdir(parents=True, exist_ok=True)
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    s = pd.read_csv(ROOT / "validation_new_organisms/scored178_B2.csv")
    blind = sorted(s[s.blind == True].accession.unique())  # noqa: E712
    build_set(blind, OUTD / "newspecies.parquet")
    validate()
    build_set(list(lab[lab.labelled == True].genome_id), ROOT / "data/processed/genome_features_metabolic_hmm.parquet")  # noqa: E712
    ext = pd.read_csv(ROOT / "data/external/labels.tsv", sep="\t")
    build_set(list(ext.genome_id), OUTD / "external.parquet")


def validate():
    mine = pd.read_parquet(OUTD / "newspecies.parquet")
    w = pd.read_csv(ROOT / "data/additional/metabolic/full_out/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet1.tsv",
                    sep="\t")
    t = table()
    w["row"] = [f"{c}|{f}|{g}" for c, f, g in zip(w.Category, w.Function, w["Gene abbreviation"])]
    cols = [c for c in w.columns if c.endswith("Hmm presence")]
    agree, n, per_row = 0, 0, {}
    for c in cols:
        g = c.replace(" Hmm presence", "")
        if g not in mine.index:
            continue
        for r, v in zip(w.row, w[c]):
            m = mine.at[g, f"metabolic_{r}"] if f"metabolic_{r}" in mine.columns else np.nan
            if np.isnan(m):
                continue
            ok = int(m) == int(v == "Present")
            agree += ok
            n += 1
            per_row.setdefault(r, []).append(ok)
    bad = sorted(((np.mean(v), k) for k, v in per_row.items()))[:8]
    print(f"agreement with real METABOLIC on {len(cols)} genomes: {agree}/{n} = {agree / n:.4f}")
    print("lowest-agreement rows:", [(round(a, 2), k) for a, k in bad])
    (OUTD / "validation.txt").write_text(f"{agree}/{n} = {agree / n:.4f}\n")


if __name__ == "__main__":
    {"build": build, "validate": validate}[sys.argv[1]]()
