"""How good *can* metrics be on these labels? Label-consistency ceiling analysis (research loop).

1. Near-identical genome pairs: proteome Jaccard on exact protein sequences (seq_hash sets). For pairs above
   0.5 (conspecific-level sharing), per-trait label discordance. Biologically near-identical genomes should
   share category labels; discordance marks annotation noise that no genome-based model can fit.
2. Identical-feature conflicts: genomes with identical energy-panel KO vectors (B2 space) but different labels;
   the minority labels in such groups are irreducible errors for any model on B2 features.
3. Annotation-consistency reference (NOT a valid predictor; uses no hold-out): each genome scored by the
   similarity-weighted labels of its k=5 most proteome-similar *other* labelled genomes. Its AUPRC shows how
   predictable the labels are from the closest relatives when nothing is held out, so an out-of-clade model
   is not expected to exceed it by much. (Observed: it is *below* the out-of-clade KO models, so relatedness alone
   does not bound performance; the dataset is effectively dereplicated.)
4. Within-genus label heterogeneity: genera with >=2 labelled genomes whose members disagree on a trait. This is
   either real strain-level variation or incomplete (closed-world) annotation; both cap achievable precision.
Output: results/qc/label_ceiling_pairs.csv, results/label_ceiling_summary.csv
"""
import numpy as np
import pandas as pd

from redox_common import PROC, RES, load_features, load_labelled, safe_ap

K = 5


def main():
    lab, traits, Y = load_labelled()
    sets = [set(pd.read_parquet(PROC / "proteins" / f"{a}.parquet", columns=["seq_hash"]).seq_hash)
            for a in lab.genome_id]
    n = len(sets)
    J = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            inter = len(sets[i] & sets[j])
            if inter:
                J[i, j] = J[j, i] = inter / len(sets[i] | sets[j])
    pairs = [(i, j, J[i, j]) for i in range(n) for j in range(i + 1, n) if J[i, j] >= 0.5]
    prow = []
    for i, j, jac in pairs:
        prow.append(dict(genome_a=lab.genome_id[i], genome_b=lab.genome_id[j], organism_a=lab.organism_name[i],
                         organism_b=lab.organism_name[j], jaccard=jac,
                         **{t: f"{int(Y[i, k])}/{int(Y[j, k])}" for k, t in enumerate(traits)}))
    (RES / "qc").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(prow).to_csv(RES / "qc/label_ceiling_pairs.csv", index=False, float_format="%.3f")

    X = load_features("energy_panel", lab.genome_id)
    key = X.astype(int).astype(str).agg("".join, axis=1).to_numpy()
    summary = []
    for k, t in enumerate(traits):
        y = Y[:, k]
        disc = [y[i] != y[j] for i, j, _ in pairs if y[i] + y[j] >= 1]
        conflicts = 0
        for g in pd.unique(key):
            m = key == g
            if m.sum() > 1:
                conflicts += min(y[m].sum(), (1 - y[m]).sum())
        W = J.copy()
        np.fill_diagonal(W, -1)
        score = np.zeros(n)
        for i in range(n):
            nb = np.argsort(-W[i])[:K]
            w = np.maximum(W[i, nb], 1e-9)
            score[i] = (w * y[nb]).sum() / w.sum()
        gen = pd.DataFrame({"g": lab.genus, "y": y}).groupby("g").y.agg(["size", "sum"])
        multi = gen[gen["size"] >= 2]
        mixed = multi[(multi["sum"] > 0) & (multi["sum"] < multi["size"])]
        summary.append(dict(trait=t, n_pos=int(y.sum()), n_genera_ge2_genomes=len(multi),
                            n_genera_with_positive=int((multi["sum"] > 0).sum()),
                            n_genera_mixed_labels=len(mixed),
                            frac_positive_genera_mixed=len(mixed) / max(int((multi["sum"] > 0).sum()), 1),
                            n_near_identical_pairs_with_a_positive=len(disc),
                            discordance_rate_near_identical=float(np.mean(disc)) if disc else np.nan,
                            identical_B2_vector_irreducible_errors=int(conflicts),
                            annotation_consistency_auprc_knn_no_holdout=safe_ap(y, score)))
    s = pd.DataFrame(summary)
    s.to_csv(RES / "label_ceiling_summary.csv", index=False, float_format="%.4f")
    print(f"near-identical pairs (Jaccard>=0.5): {len(pairs)}")
    print(s.to_string())
    print("macro annotation-consistency AUPRC (no hold-out, kNN proteome):",
          round(s.annotation_consistency_auprc_knn_no_holdout.mean(), 4))


if __name__ == "__main__":
    main()
