"""Independent check that external genomes share no near-identical proteome with the training set (audit 4, m-6).

Jaccard of exact protein-sequence sets (seq_hash) between every novel external genome and every labelled training
genome. Output: results/qc/external_novelty_jaccard.json
"""
import json

import pandas as pd

from redox_common import PROC, RES, ROOT, load_labelled


def main():
    lab, _, _ = load_labelled()
    ext = pd.read_csv(ROOT / "data/external/lineage.tsv", sep="\t")
    ext = ext[ext.novel]
    train = {a: set(pd.read_parquet(PROC / "proteins" / f"{a}.parquet", columns=["seq_hash"]).seq_hash)
             for a in lab.genome_id}
    best = []
    for a in ext.genome_id:
        s = set(pd.read_parquet(ROOT / f"data/external/proteins/{a}.parquet", columns=["seq_hash"]).seq_hash)
        j, who = 0.0, None
        for b, t in train.items():
            inter = len(s & t)
            if inter:
                v = inter / len(s | t)
                if v > j:
                    j, who = v, b
        best.append((a, who, j))
    d = pd.DataFrame(best, columns=["external_genome", "closest_training_genome", "jaccard"])
    res = {"n_external_novel": len(d), "n_training": len(train), "max_jaccard": float(d.jaccard.max()),
           "mean_jaccard": float(d.jaccard.mean()), "n_above_0.05": int((d.jaccard > 0.05).sum()),
           "n_above_0.2": int((d.jaccard > 0.2).sum()),
           "argmax": d.loc[d.jaccard.idxmax(), ["external_genome", "closest_training_genome"]].to_dict()}
    (RES / "qc").mkdir(parents=True, exist_ok=True)
    (RES / "qc/external_novelty_jaccard.json").write_text(json.dumps(res, indent=2))
    d.sort_values("jaccard", ascending=False).head(20).to_csv(RES / "qc/external_novelty_top20.csv", index=False)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
