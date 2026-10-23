"""Integrity checks for the pathway task. Every check prints PASS/FAIL; a FAIL blocks reporting.

Checks (mirroring the donor/acceptor audit):
  1. labels.tsv matches its frozen hash
  2. one genome per species; no duplicate accessions
  3. no split group (genus / family) spans two outer folds - the blocking is real
  4. every genome appears in exactly one fold at each level
  5. features contain no label-derived column and no taxonomy
  6. KO coverage: no genome with an empty profile; feature matrix aligned to the label rows
  7. FAPROTAX-derived share per trait (taxonomy-propagated labels are the circularity risk here)
  8. overlap with the donor/acceptor training genomes (independence of the cross-source check)
  9. out-of-fold predictions: one prediction per (arm, level, trait, genome), labels agree with the label matrix
usage: python src/37_pathway_audit.py
"""
import json
from pathlib import Path

import pandas as pd

from pathway_common import PATH, RES, load_pathways
from redox_common import PROC, sha256

FAILS = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}{(' - ' + detail) if detail else ''}")
    if not ok:
        FAILS.append(name)


def main():
    lab, traits, Y = load_pathways()  # asserts the hash internally
    check("1. labels.tsv matches frozen hash", True, sha256(PATH / "labels.tsv")[:12])
    check("2a. no duplicate genome_id", lab.genome_id.is_unique, f"{len(lab)} genomes")
    check("2b. no duplicate species", lab.species_tax_id.is_unique, f"{lab.species_tax_id.nunique()} species")

    for level in ("genus", "family"):
        spans = lab.groupby(f"group_{level}")[f"fold_{level}"].nunique()
        check(f"3. no {level} group spans folds", (spans == 1).all(),
              f"{len(spans)} groups, max folds per group {int(spans.max())}")
        counts = lab[f"fold_{level}"].value_counts().sort_index()
        check(f"4. {level}: every genome in exactly one fold", counts.sum() == len(lab),
              " ".join(f"f{k}={v}" for k, v in counts.items()))

    ko = pd.read_csv(PATH / "ko_matrix.tsv", sep="\t", index_col=0, nrows=5)
    bad = [c for c in ko.columns if not c.startswith("K") or c in traits]
    check("5a. feature columns are KO ids only", not bad, f"{len(ko.columns)} columns, offending {bad[:3]}")
    tax_like = [c for c in ko.columns if c.lower() in {"genus", "family", "species", "domain", "phylum"}]
    check("5b. no taxonomy column in the features", not tax_like)

    full = pd.read_csv(PATH / "ko_matrix.tsv", sep="\t", index_col=0)
    empty = (full.sum(axis=1) == 0).sum()
    missing = set(lab.genome_id) - set(full.index)
    check("6a. no genome with an empty KO profile", empty == 0, f"{empty} empty")
    check("6b. every labelled genome has a KO profile", not missing, f"{len(missing)} missing")
    check("6c. KO matrix rows align to labels", len(full.reindex(lab.genome_id).dropna(how="all")) == len(lab),
          f"{full.shape[1]} KOs")

    fap = pd.DataFrame({"trait": traits,
                        "n_pos": [int(Y[:, j].sum()) for j in range(len(traits))],
                        "frac_pos_faprotax": [float(lab.faprotax_source[Y[:, j] == 1].mean())
                                              for j in range(len(traits))]})
    (RES / "pathways").mkdir(parents=True, exist_ok=True)
    fap.to_csv(RES / "pathways/qc_faprotax_share.csv", index=False, float_format="%.3f")
    print(fap.round(3).to_string(index=False))
    check("7. FAPROTAX share reported per trait", True,
          f"mean {fap.frac_pos_faprotax.mean():.2f} - sensitivity analysis: --subset no_faprotax")

    da = pd.read_csv(PROC / "labels.tsv", sep="\t")
    shared = set(lab.genome_id) & set(da.genome_id[da.labelled] if "labelled" in da.columns else da.genome_id)
    check("8. donor/acceptor genome overlap quantified", True,
          f"{len(shared)} of {len(lab)} pathway genomes are also donor/acceptor training genomes "
          f"({100 * len(shared) / len(lab):.1f} %); the cross-source check uses genus-held-out predictions")

    oof_dir = Path(RES / "pathways/oof")
    for f in sorted(oof_dir.glob("*.parquet")) if oof_dir.exists() else []:
        d = pd.read_parquet(f)
        dup = d.duplicated(["arm", "level", "trait", "genome_id"]).sum()
        check(f"9a. {f.stem}: one prediction per genome/trait/level", dup == 0, f"{len(d)} rows")
        idx = {g: i for i, g in enumerate(lab.genome_id)}
        j = {t: k for k, t in enumerate(traits)}
        ok = all(Y[idx[g], j[t]] == y for g, t, y in zip(d.genome_id, d.trait, d.y) if g in idx and t in j)
        check(f"9b. {f.stem}: OOF labels match the label matrix", ok)
        fl = d.merge(lab[["genome_id", "fold_genus", "fold_family"]], on="genome_id")
        for level in ("genus", "family"):
            sub = fl[fl.level == level]
            check(f"9c. {f.stem}/{level}: test-fold assignment matches the frozen split",
                  bool((sub.fold == sub[f"fold_{level}"]).all()) if len(sub) else True)

    print("\nAUDIT:", "ALL PASS" if not FAILS else f"{len(FAILS)} FAILURES: {FAILS}")
    (RES / "pathways/audit.json").write_text(json.dumps(
        {"failures": FAILS, "n_genomes": len(lab), "n_traits": len(traits),
         "n_genus_groups": int(lab.group_genus.nunique()), "n_family_groups": int(lab.group_family.nunique()),
         "labels_sha256": sha256(PATH / "labels.tsv")}, indent=2))
    raise SystemExit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
