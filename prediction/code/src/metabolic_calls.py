"""Category calls from the pyhmmer re-implementation of METABOLIC-G (src/103_metabolic_features.py, per-HMM-row
presence, columns 'metabolic_<Category>|<Function>|<Gene>'), with the pre-registered mapping of src/95_metabolic_score.py.
A function is present when any of its rows is present; sulfite reduction also counts dsrA/dsrB and
sulfur oxidation does not, as METABOLIC's worksheet 2 does (verified: 100% agreement with real METABOLIC output, src/118_metabolic_benchmarks.py validate)."""
import importlib

import pandas as pd

m95 = importlib.import_module("95_metabolic_score")


def calls_from_rows(R):
    parts = [c[len("metabolic_"):].split("|") for c in R.columns]
    cat = pd.Series([p[0] for p in parts], index=R.columns)
    fun = pd.Series([p[1].strip() if len(p) > 1 else "" for p in parts], index=R.columns)
    out = {}
    for k, pairs in m95.MAP.items():
        funs = {f for _, f in pairs}
        cols = [c for c in R.columns if fun[c] in funs]
        if "Sulfur oxidation" in funs:  # ...and does not count dsrA/dsrB towards sulfur oxidation
            cols = [c for c in cols if not (fun[c] == "Sulfur oxidation" and c.split("|")[-1] in ("dsrA", "dsrB"))]
        if "Sulfite reduction" in funs:  # METABOLIC's worksheet 2 also counts dsrA/dsrB towards sulfite reduction
            cols += [c for c in R.columns if cat[c] == "Sulfur cycling" and c.split("|")[-1] in ("dsrA", "dsrB")]
        assert cols, k
        out[k] = (R[cols] > 0).any(axis=1).astype(int)
    ox = [c for c in R.columns if cat[c].startswith("Oxygen metabolism") and any(x in fun[c] for x in m95.OXY_KEYS)]
    out["acceptor_oxygen"] = (R[ox] > 0).any(axis=1).astype(int)
    return pd.DataFrame(out)
