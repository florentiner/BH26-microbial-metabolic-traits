"""Is "one specific KO beats the whole group" real, or selection bias?

Test B picked, per functional group, the single KO with the best MCC - and scored it on the same genomes it was
picked on, which is optimistic by construction. Here the choice is made on training genera only and scored on
held-out genera (5 folds grouped by genus, 20 repeats with different fold assignments). If the advantage of the
single KO survives, the "definition too broad" conclusion is real; the selection stability (how often the same
KO is chosen) says whether the choice is a property of the biology or of the sample.

Output: results/jaccard/single_ko_crossval.csv
usage: python src/43_single_ko_crossval.py
"""
import importlib.util
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/jaccard"
SEED = 20260918
REPEATS = 20

spec = importlib.util.spec_from_file_location("jc", ROOT / "src/40_jaccard_correlations.py")
jc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(jc)


def pooled_mcc(pred, y):
    m, _ = jc.mcc_phi(pred, y)
    return m


def main():
    X, groups, role, ranks, info = jc.load_html(ROOT / "jaccard_similarity.html")
    lab = pd.read_csv(ROOT / "data/processed/labels.tsv", sep="\t")
    lab = lab[lab.labelled] if "labelled" in lab.columns else lab
    lab = lab.set_index("genome_id")
    ids = [g for g in X.index if g in lab.index]
    genera = ranks["genus"].reindex(ids).fillna("?").to_numpy()
    rng = np.random.default_rng(SEED)
    rows = []
    for name, kos in groups.items():
        trait = jc.GROUP_TRAIT[name]
        if trait not in lab.columns:
            continue
        y = lab.loc[ids, trait].to_numpy(dtype=int)
        Xg = X.loc[ids, kos].to_numpy()
        anyko = Xg.max(axis=1)
        mcc_single, mcc_any, chosen = [], [], Counter()
        for r in range(REPEATS):
            # re-label genera with a random permutation so GroupKFold gives a different split each repeat
            ug = np.unique(genera)
            remap = dict(zip(ug, rng.permutation(len(ug))))
            gcode = np.array([remap[g] for g in genera])
            pred_single = np.zeros(len(y), dtype=int)
            for tr, te in GroupKFold(n_splits=5).split(Xg, y, gcode):
                best_j, best_m = 0, -2.0
                for j in range(Xg.shape[1]):
                    if Xg[tr, j].sum() == 0:
                        continue
                    m = pooled_mcc(Xg[tr, j], y[tr])
                    if not np.isnan(m) and m > best_m:
                        best_j, best_m = j, m
                chosen[kos[best_j]] += 1
                pred_single[te] = Xg[te, best_j]
            mcc_single.append(pooled_mcc(pred_single, y))
            mcc_any.append(pooled_mcc(anyko, y))  # the group rule has nothing to select, so no CV needed
        ms, ma = np.array(mcc_single, dtype=float), np.array(mcc_any, dtype=float)
        top, top_n = chosen.most_common(1)[0]
        rows.append(dict(functional_group=name, phenotype=trait, n_kos=len(kos), n_pos=int(y.sum()),
                         mcc_any_ko_rule=float(np.nanmean(ma)),
                         mcc_single_ko_out_of_genus=float(np.nanmean(ms)),
                         mcc_single_ko_ci_low=float(np.nanpercentile(ms, 2.5)),
                         mcc_single_ko_ci_high=float(np.nanpercentile(ms, 97.5)),
                         gain_out_of_genus=float(np.nanmean(ms) - np.nanmean(ma)),
                         most_chosen_ko=top, selection_stability=top_n / sum(chosen.values()),
                         kos_ever_chosen=", ".join(f"{k}:{v}" for k, v in chosen.most_common(4))))
    out = pd.DataFrame(rows).sort_values("gain_out_of_genus", ascending=False)
    out["verdict"] = np.where(out.mcc_single_ko_ci_low > out.mcc_any_ko_rule, "single KO better (robust)",
                              np.where(out.mcc_single_ko_ci_high < out.mcc_any_ko_rule, "group rule better",
                                       "no robust difference"))
    out.to_csv(OUT / "single_ko_crossval.csv", index=False, float_format="%.4f")
    print(out[["functional_group", "n_pos", "mcc_any_ko_rule", "mcc_single_ko_out_of_genus",
               "mcc_single_ko_ci_low", "mcc_single_ko_ci_high", "gain_out_of_genus", "most_chosen_ko",
               "selection_stability", "verdict"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
