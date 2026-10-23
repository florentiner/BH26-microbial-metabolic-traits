"""Concordance of marker-KO annotations: provided table vs local KOfam reimplementation vs TogoMCP (BH26 RDF).

Output: results/qc/togo_concordance.json
"""
import json

import pandas as pd

from redox_common import PROC, RES, load_features


def main():
    log = pd.read_csv(PROC / "togo_fetch_log_internal.tsv", sep="\t")
    ok = log.loc[log.status == "ok", "genome_id"]
    prov = load_features("markers_provided", ok)
    ours = load_features("markers_recomputed", ok)[prov.columns]
    t = pd.read_parquet(PROC / "togo_kofam_internal.parquet")
    t = t[t.KO.isin(prov.columns)]
    res = {"n_genomes_in_togo": int(len(ok)), "n_genomes_missing_in_togo": int((log.status != "ok").sum()),
           "n_marker_kos": int(prov.shape[1])}
    for name, df in (("togo_significant", t[t.significant]), ("togo_any_hit", t)):
        m = pd.crosstab(df.genome_id, df.KO).clip(upper=1).reindex(index=ok, columns=prov.columns).fillna(0)
        for ref_name, ref in (("provided", prov), ("local_kofam", ours)):
            res[f"{name}_vs_{ref_name}"] = {
                "agreement": float((ref.values == m.values).mean()),
                f"{ref_name}1_{name}0": int(((ref == 1) & (m == 0)).values.sum()),
                f"{ref_name}0_{name}1": int(((ref == 0) & (m == 1)).values.sum())}
    (RES / "qc").mkdir(parents=True, exist_ok=True)
    (RES / "qc/togo_concordance.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
