"""Full-genome KofamScan KO profiles from the RDF Portal BH26 Microbes SPARQL endpoint (the data behind TogoMCP).

For every genome: all protein->KO KofamScan statements with score, threshold, E-value, significance flag.
Used for (a) baseline 2 (LLM over TogoMCP KO evidence) and (b) arm B3 (full significant-KO profile).
Output: data/processed/togo_kofam_<set>.parquet   (genome_id, protein, KO, score, threshold, evalue, significant)
        data/processed/togo_ko_matrix_<set>.tsv    (genome x KO, significant hits only)
        data/processed/togo_fetch_log_<set>.tsv
usage: python src/20_togo_ko_profiles.py [--set internal|external]
"""
import json
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ENDPOINT = "https://rdfportal.org/microbes/sparql"
QUERY = """PREFIX sio: <http://semanticscience.org/resource/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
SELECT ?prot ?ko ?score ?thr ?ev ?sig WHERE {
 ?prot <http://purl.obolibrary.org/obo/RO_0002162> <http://identifiers.org/refseq.gcf/%s> .
 ?st rdf:subject ?prot ; rdf:object ?ko ; sio:SIO_000008 ?a1, ?a2, ?a3, ?a4 .
 ?a1 rdfs:label "significant" ; sio:SIO_000300 ?sig .
 ?a2 rdfs:label "score" ; sio:SIO_000300 ?score .
 ?a3 rdfs:label "threshold" ; sio:SIO_000300 ?thr .
 ?a4 rdfs:label "E-value" ; sio:SIO_000300 ?ev .
}"""


def fetch(acc):
    data = urllib.parse.urlencode({"query": QUERY % acc}).encode()
    for i in range(6):
        try:
            req = urllib.request.Request(ENDPOINT, data=data,
                                         headers={"Accept": "application/sparql-results+json"})
            with urllib.request.urlopen(req, timeout=600) as r:
                rows = json.load(r)["results"]["bindings"]
            out = [(acc, b["prot"]["value"].rsplit("/", 1)[-1], b["ko"]["value"].rsplit("/", 1)[-1],
                    float(b["score"]["value"]), float(b["thr"]["value"]) if b["thr"]["value"] not in ("", "-") else None,
                    b["ev"]["value"], b["sig"]["value"] == "true") for b in rows]
            return acc, out, "ok" if out else "empty"
        except Exception as e:
            err = str(e)
            time.sleep(10 * (i + 1))
    return acc, [], f"error:{err}"


def main():
    which = sys.argv[sys.argv.index("--set") + 1] if "--set" in sys.argv else "internal"
    if which == "internal":
        ids = pd.read_csv(ROOT / "ko_presence_586genomes.tsv", sep="\t", encoding="utf-8-sig").genome_id.tolist()
    else:
        st = pd.read_csv(ROOT / "data/external/download_manifest.tsv", sep="\t")
        ids = st.loc[st.status == "ok", "genome_id"].tolist()
    rows, log = [], []
    with ThreadPoolExecutor(3) as ex:
        for k, (acc, out, status) in enumerate(ex.map(fetch, ids)):
            rows += out
            log.append((acc, status, len(out)))
            if k % 50 == 0:
                print(k, acc, status, len(out), flush=True)
    df = pd.DataFrame(rows, columns=["genome_id", "protein", "KO", "score", "threshold", "evalue", "significant"])
    df.to_parquet(ROOT / f"data/processed/togo_kofam_{which}.parquet", index=False)
    lg = pd.DataFrame(log, columns=["genome_id", "status", "n_statements"])
    lg.to_csv(ROOT / f"data/processed/togo_fetch_log_{which}.tsv", sep="\t", index=False)
    sig = df[df.significant].drop_duplicates(["genome_id", "KO"])
    mat = pd.crosstab(sig.genome_id, sig.KO).clip(upper=1).reindex(ids).fillna(0).astype("int8")
    mat.index.name = "genome_id"
    mat.to_csv(ROOT / f"data/processed/togo_ko_matrix_{which}.tsv", sep="\t")
    print(lg.status.value_counts().to_dict(), "matrix", mat.shape)


if __name__ == "__main__":
    main()
