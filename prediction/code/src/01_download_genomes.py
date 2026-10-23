"""Download RefSeq proteome (protein.faa) and annotation (GFF3) for every genome_id.

Input : ko_presence_586genomes.tsv (genome_id column)
Output: data/raw/genomes/<acc>/{protein.faa,genomic.gff,assembly_report.json}
        data/raw/genomes/download_manifest.tsv  (genome_id, status, n_proteins, faa_sha256)
Idempotent: a genome is complete only when both files exist (written atomically via rename).
"""
import hashlib
import io
import sys
import time
import urllib.error
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/raw/genomes"
API = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
       "?include_annotation_type=PROT_FASTA&include_annotation_type=GENOME_GFF")
N_WORKERS = 4  # NCBI anonymous limit is 5 requests/s


def fetch(url, tries=6):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=600) as r:
                return r.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            wait = 10 * (i + 1)
            print(f"  retry {i + 1} for {url[:90]} after {e}; sleeping {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError(f"failed after {tries} tries: {url}")


def download_one(acc):
    d = OUT / acc
    faa, gff = d / "protein.faa", d / "genomic.gff"
    if not (faa.exists() and gff.exists()):
        d.mkdir(parents=True, exist_ok=True)
        try:
            z = zipfile.ZipFile(io.BytesIO(fetch(API.format(acc=acc))))
        except Exception as e:  # recorded in the manifest, run exits non-zero
            print(acc, "ERROR", e, flush=True)
            return acc, f"error:{e}", 0, ""
        members = {n.rsplit("/", 1)[-1]: n for n in z.namelist()}
        missing = [b for b in ("protein.faa", "genomic.gff") if b not in members]
        if missing:
            print(acc, "MISSING", missing, flush=True)
            return acc, "missing:" + ",".join(missing), 0, ""
        if "assembly_data_report.jsonl" in members:
            (d / "assembly_report.json").write_bytes(z.read(members["assembly_data_report.jsonl"]))
        for base, target in (("genomic.gff", gff), ("protein.faa", faa)):
            tmp = target.with_suffix(target.suffix + ".tmp")
            tmp.write_bytes(z.read(members[base]))
            tmp.rename(target)
        time.sleep(0.5)
    data = faa.read_bytes()
    return acc, "ok", data.count(b">"), hashlib.sha256(data).hexdigest()


def main():
    ids = pd.read_csv(ROOT / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig")["genome_id"].tolist()
    rows = []
    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        for k, row in enumerate(ex.map(download_one, ids)):
            rows.append(row)
            if k % 25 == 0:
                print(f"{k + 1}/{len(ids)} {row}", flush=True)
    m = pd.DataFrame(rows, columns=["genome_id", "status", "n_proteins", "faa_sha256"])
    m.to_csv(OUT / "download_manifest.tsv", sep="\t", index=False)
    print(m.status.value_counts())
    print(m.n_proteins.describe())
    sys.exit(0 if (m.status == "ok").all() else 1)


if __name__ == "__main__":
    main()
