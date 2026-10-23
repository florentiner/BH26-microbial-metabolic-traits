"""KOfam HMM annotation of every protein for a pre-specified KO panel (KofamScan logic).

Panel (fixed before any modelling): all KOs under KEGG BRITE ko00001 category
"09102 Energy metabolism" UNION the 84 marker KOs in ko_presence_586genomes.tsv.
A full 27,866-profile KofamScan over 1.8 M proteins is not tractable on the local machine.

Assignment rule = KofamScan default (reported with '*'): hit if
  score_type == 'full'   -> full-sequence bit score >= KO threshold
  score_type == 'domain' -> best-domain bit score   >= KO threshold
KOs without a predefined threshold are not assigned (KofamScan default without -E).
Bit scores do not depend on database size, so searching all unique sequences at
once is equivalent to per-genome KofamScan runs.

Output:
  data/processed/kofam_hits.parquet        seq_hash, KO, score, threshold
  data/processed/ko_panel.tsv              KO, in_energy_brite, is_marker, threshold, score_type, definition
  data/processed/ko_matrix_panel.tsv       genome x KO presence (0/1) for the panel
"""
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyhmmer

ROOT = Path(__file__).resolve().parents[1]
KOFAM = ROOT / "data/raw/kofam"
CPUS = 8


def energy_kos():
    tree = json.load(open(KOFAM / "ko00001.json"))
    out = set()

    def walk(node, in_energy):
        in_energy = in_energy or node["name"].startswith("09102")
        for c in node.get("children", []):
            walk(c, in_energy)
        if "children" not in node and in_energy:
            m = re.match(r"(K\d{5})", node["name"])
            if m:
                out.add(m.group(1))

    walk(tree, False)
    return out


def build_panel():
    tsv = pd.read_csv(ROOT / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig")
    markers = {c.split()[0] for c in tsv.columns if re.match(r"K\d{5} ", c)}
    energy = energy_kos()
    ko_list = pd.read_csv(KOFAM / "ko_list.gz", sep="\t")
    ko_list = ko_list[ko_list.knum.isin(energy | markers)].copy()
    ko_list["threshold"] = pd.to_numeric(ko_list.threshold, errors="coerce")
    panel = ko_list[["knum", "threshold", "score_type", "definition"]].rename(columns={"knum": "KO"})
    panel["in_energy_brite"] = panel.KO.isin(energy)
    panel["is_marker"] = panel.KO.isin(markers)
    missing_markers = markers - set(panel.KO)
    print(f"panel KOs={len(panel)} (energy={len(energy)}, markers={len(markers)}, markers missing from KOfam={missing_markers})")
    return tsv, panel


def load_hmms(panel):
    scorable = panel.dropna(subset=["threshold"])
    hmms = []
    for ko in scorable.KO:
        with pyhmmer.plan7.HMMFile(KOFAM / "profiles" / f"{ko}.hmm") as fh:
            hmms.append(fh.read())
    return hmms, dict(zip(scorable.KO, scorable.threshold)), dict(zip(scorable.KO, scorable.score_type))


def scan(hmms, thr, stype, hashes, sequences, cpus=CPUS, log_every=50):
    """KofamScan-rule hits for the given sequences -> DataFrame(seq_hash, KO, score, threshold)."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = pyhmmer.easel.TextSequenceBlock(
        pyhmmer.easel.TextSequence(name=h.encode(), sequence=s) for h, s in zip(hashes, sequences)
    ).digitize(alphabet)
    rows, t0 = [], time.time()
    for k, hits in enumerate(pyhmmer.hmmer.hmmsearch(hmms, seqs, cpus=cpus, E=10.0)):
        ko = hits.query.name.decode() if isinstance(hits.query.name, bytes) else hits.query.name
        ko = ko.split(".")[0]
        t = thr[ko]
        for h in hits:
            score = h.score if stype[ko] == "full" else h.best_domain.score
            if score >= t:
                name = h.name.decode() if isinstance(h.name, bytes) else h.name
                rows.append((name, ko, float(score), t))
        if log_every and k % log_every == 0:
            print(f"{k + 1}/{len(hmms)} {ko} elapsed {time.time() - t0:.0f}s hits_total={len(rows)}", flush=True)
    return pd.DataFrame(rows, columns=["seq_hash", "KO", "score", "threshold"])


def presence_matrix(hits, genome_hashes, kos):
    """genome_hashes: dict genome_id -> iterable of seq_hash. Returns genome x KO 0/1 DataFrame."""
    ko_of = hits.groupby("seq_hash").KO.apply(set).to_dict()
    mat = pd.DataFrame(0, index=pd.Index(list(genome_hashes), name="genome_id"), columns=sorted(kos), dtype=np.int8)
    for acc, hs in genome_hashes.items():
        present = set().union(*[ko_of.get(h, set()) for h in set(hs)])
        mat.loc[acc, sorted(present & set(kos))] = 1
    return mat


def main():
    tsv, panel = build_panel()
    panel.to_csv(ROOT / "data/processed/ko_panel.tsv", sep="\t", index=False)
    print("KOs without threshold (never assigned):", sorted(set(panel.KO) - set(panel.dropna(subset=["threshold"]).KO)))
    hmms, thr, stype = load_hmms(panel)
    u = pd.read_parquet(ROOT / "data/processed/unique_sequences.parquet", columns=["seq_hash", "sequence"])
    print("sequences", len(u), flush=True)
    hits = scan(hmms, thr, stype, u.seq_hash.tolist(), u.sequence.tolist())
    del u
    hits.to_parquet(ROOT / "data/processed/kofam_hits.parquet", index=False)
    genome_hashes = {acc: pd.read_parquet(ROOT / f"data/processed/proteins/{acc}.parquet", columns=["seq_hash"]).seq_hash
                     for acc in tsv.genome_id}
    mat = presence_matrix(hits, genome_hashes, panel.KO)
    mat.to_csv(ROOT / "data/processed/ko_matrix_panel.tsv", sep="\t")
    print("done", mat.shape, "KOs present in >=1 genome:", int((mat.sum() > 0).sum()))


if __name__ == "__main__":
    main()
