"""Energy-metabolism pathway calls from the proteome: KEGG module completeness (project deliverable 3).

The supplied labels are organism-level donor/acceptor *categories*, so no pathway-level model can be trained from
them. Pathway presence is instead a deterministic readout of the same KO annotation: for every KEGG module in
CLASS "Energy metabolism" (plus carbon-fixation and methane modules), the module's DEFINITION is parsed and the
fraction of satisfied steps is computed from the genome's KOs.

Definition grammar (KEGG): space = consecutive steps (AND), comma = alternatives (OR), '+' = required complex
subunit (AND), '-' = optional component (ignored), parentheses group. A step counts as present if its expression
evaluates true with the genome's KO set; module completeness = satisfied steps / total steps.

Outputs:
  data/processed/kegg_modules_energy.tsv         module, name, class, n_steps, definition
  results/pathway_module_completeness.tsv        genome x module completeness (0-1), all 586 genomes
  results/pathway_calls.csv                      long form: complete (>=0.8) / partial (0.5-0.8) / absent modules
  results/qc/pathway_vs_labels.csv               agreement of selected modules with the curated trait labels
usage: python src/29_pathway_modules.py [--min-complete 0.8]
"""
import json
import re
import sys
import time
import urllib.request

import numpy as np
import pandas as pd

from redox_common import PROC, RES, ROOT, load_labelled
from sklearn.metrics import average_precision_score

KEGG = "https://rest.kegg.jp"
CACHE = ROOT / "data/raw/kegg_modules"
KEEP_CLASS = ("Energy metabolism", "Carbon fixation", "Methane metabolism", "Nitrogen metabolism", "Sulfur metabolism")
# modules whose presence should track a curated trait, used as a sanity check (not a prediction)
MODULE_TRAIT = {
    "M00567": "acceptor_CO2_methanogen_acetogen",   # methanogenesis, CO2 => CH4
    "M00357": "acceptor_CO2_methanogen_acetogen",   # methanogenesis, acetate => CH4
    "M00377": "acceptor_CO2_methanogen_acetogen",   # Wood-Ljungdahl (acetogens)
    "M00596": "acceptor_sulfate",                   # dissimilatory sulfate reduction
    "M00529": "acceptor_nitrate",                   # denitrification
    "M00530": "acceptor_nitrate",                   # dissimilatory nitrate reduction to ammonia
    "M00528": "donor_ammonia_oxidation",            # nitrification
    "M00804": "donor_nitrite_oxidation",            # comammox (ammonia => nitrate)
    "M00595": "donor_sulfur_oxidation",             # thiosulfate oxidation by SOX complex
}


def fetch(url, path):
    if path.exists():
        return path.read_text()
    for i in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                t = r.read().decode()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(t)
            time.sleep(0.34)  # KEGG REST courtesy limit
            return t
        except Exception:
            time.sleep(5 * (i + 1))
    raise RuntimeError(url)


def parse_steps(definition):
    """Split a KEGG module DEFINITION into top-level steps (space-separated at depth 0)."""
    steps, depth, cur = [], 0, []
    for ch in definition:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == " " and depth == 0:
            if cur:
                steps.append("".join(cur))
                cur = []
        else:
            cur.append(ch)
    if cur:
        steps.append("".join(cur))
    return [s for s in steps if s.strip("-")]


def step_present(expr, kos):
    """Evaluate one step: ',' = OR, '+' = AND (complex), '-' = optional (dropped), parentheses group."""
    expr = expr.strip()
    while expr.startswith("(") and expr.endswith(")"):
        depth = 0
        balanced = True
        for i, ch in enumerate(expr):
            depth += (ch == "(") - (ch == ")")
            if depth == 0 and i < len(expr) - 1:
                balanced = False
                break
        if not balanced:
            break
        expr = expr[1:-1].strip()

    def split_top(s, sep):
        out, depth, cur = [], 0, []
        for ch in s:
            depth += (ch == "(") - (ch == ")")
            if ch == sep and depth == 0:
                out.append("".join(cur))
                cur = []
            else:
                cur.append(ch)
        out.append("".join(cur))
        return out

    alts = split_top(expr, ",")
    if len(alts) > 1:
        return any(step_present(a, kos) for a in alts)
    parts = split_top(expr, "+")
    if len(parts) > 1:
        required = [parts[0]] + [p for p in parts[1:] if not p.startswith("-")]
        return all(step_present(p.lstrip("-"), kos) for p in required)
    opt = split_top(expr, "-")
    if len(opt) > 1:
        return step_present(opt[0], kos)
    m = re.findall(r"K\d{5}", expr)
    return bool(m) and any(k in kos for k in m)


def main():
    min_complete = float(sys.argv[sys.argv.index("--min-complete") + 1]) if "--min-complete" in sys.argv else 0.8
    CACHE.mkdir(parents=True, exist_ok=True)
    listing = fetch(f"{KEGG}/list/module", CACHE / "module_list.txt")
    ids = [ln.split("\t")[0] for ln in listing.strip().split("\n")]
    rows = []
    for k, mid in enumerate(ids):
        txt = fetch(f"{KEGG}/get/{mid}", CACHE / f"{mid}.txt")
        name = re.search(r"^NAME\s+(.+)$", txt, re.M)
        cls = re.search(r"^CLASS\s+(.+)$", txt, re.M)
        dfn = re.search(r"^DEFINITION\s+(.+(?:\n\s{12}.+)*)$", txt, re.M)
        if not (name and cls and dfn):
            continue
        cls_s = cls.group(1)
        if not any(c in cls_s for c in KEEP_CLASS):
            continue
        definition = " ".join(dfn.group(1).split())
        rows.append(dict(module=mid, name=name.group(1).strip(), kegg_class=cls_s.strip(),
                         definition=definition, n_steps=len(parse_steps(definition))))
        if k % 100 == 0:
            print(f"{k}/{len(ids)} modules scanned, kept {len(rows)}", flush=True)
    mods = pd.DataFrame(rows)
    mods.to_csv(PROC / "kegg_modules_energy.tsv", sep="\t", index=False)
    print(f"energy-related modules: {len(mods)}")

    togo = PROC / "togo_ko_matrix_internal.tsv"  # full KO profile where available, else the 852-KO panel
    panel = pd.read_csv(PROC / "ko_matrix_panel.tsv", sep="\t", index_col=0)
    full = pd.read_csv(togo, sep="\t", index_col=0) if togo.exists() else None
    ko_sets = {}
    for acc in panel.index:
        s = set(panel.columns[panel.loc[acc] > 0])
        if full is not None and acc in full.index:
            s |= set(full.columns[full.loc[acc] > 0])
        ko_sets[acc] = s

    comp = pd.DataFrame(index=panel.index, columns=mods.module, dtype=float)
    for _, m in mods.iterrows():
        steps = parse_steps(m.definition)
        for acc, kos in ko_sets.items():
            comp.at[acc, m.module] = np.mean([step_present(st, kos) for st in steps]) if steps else np.nan
    comp.index.name = "genome_id"
    comp.to_csv(RES / "pathway_module_completeness.tsv", sep="\t", float_format="%.3f")

    name_of = dict(zip(mods.module, mods.name))
    long = comp.stack().reset_index()
    long.columns = ["genome_id", "module", "completeness"]
    long["module_name"] = long.module.map(name_of)
    long["call"] = np.where(long.completeness >= min_complete, "complete",
                            np.where(long.completeness >= 0.5, "partial", "absent"))
    long[long.call != "absent"].sort_values(["genome_id", "completeness"], ascending=[True, False]).to_csv(
        RES / "pathway_calls.csv", index=False, float_format="%.3f")

    lab, traits, Y = load_labelled()
    chk = []
    for mid, trait in MODULE_TRAIT.items():
        if mid not in comp.columns or trait not in traits:
            continue
        y = Y[:, traits.index(trait)]
        c = comp.loc[lab.genome_id, mid].to_numpy(dtype=float)
        chk.append(dict(module=mid, module_name=name_of[mid], trait=trait, n_pos=int(y.sum()),
                        mean_completeness_positives=float(np.nanmean(c[y == 1])),
                        mean_completeness_negatives=float(np.nanmean(c[y == 0])),
                        auprc_completeness_vs_label=float(average_precision_score(y, np.nan_to_num(c))),
                        frac_positives_complete=float(np.nanmean(c[y == 1] >= min_complete)),
                        frac_negatives_complete=float(np.nanmean(c[y == 0] >= min_complete))))
    chk = pd.DataFrame(chk)
    (RES / "qc").mkdir(parents=True, exist_ok=True)
    chk.to_csv(RES / "qc/pathway_vs_labels.csv", index=False, float_format="%.3f")
    print(chk.round(3).to_string(index=False))
    print(f"\nmodules complete (>= {min_complete}) per genome: median "
          f"{int(np.median((comp >= min_complete).sum(axis=1)))}")
    json.dump({"n_modules": int(len(mods)), "min_complete": min_complete,
               "median_complete_modules_per_genome": int(np.median((comp >= min_complete).sum(axis=1)))},
              open(RES / "qc/pathway_modules_summary.json", "w"), indent=2)


if __name__ == "__main__":
    main()
