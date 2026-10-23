"""Per-trait 95% intervals for Supplementary Table S9 (pathway model trained only on the phenotype collection).

Re-runs the cross-validation part of src/124_pathway_hackathon_only.py (same genomes, folds, configs and seeds) and
bootstraps the per-trait AUPRC difference (collection-only minus final EnergyKO-M) over genera, as for the macro value.
Output: results/additional/pathway_hackathon_only/per_trait_ci.csv
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pathway_common import ko_matrix, load_pathways, module_completeness  # noqa: E402
from redox_common import load_labelled, nested_cv, safe_ap  # noqa: E402

arms = importlib.import_module("32_pathway_arms")
SEED, B = 20261018, 2000
lab, traits, Y = load_pathways()
da, _, _ = load_labelled()
keep = lab.genome_id.isin(set(da.genome_id)).to_numpy()
labS, YS = lab[keep].reset_index(drop=True), Y[keep]
use = [t for j, t in enumerate(traits) if YS[:, j].sum() >= 25]
J = [traits.index(t) for t in use]
Xdf = pd.concat([ko_matrix(labS.genome_id), module_completeness(labS.genome_id)], axis=1)
Xdf = Xdf.loc[:, Xdf.std(axis=0) > 0]
X = Xdf.to_numpy(dtype=np.float32)
outs = Parallel(n_jobs=8)(delayed(nested_cv)("PB3_hack", X, labS, [t], YS[:, [j]], arms.configs(), arms.fit_predict, "genus")
                          for t, j in zip(use, J))
oofS = pd.concat(outs).groupby(["trait", "genome_id"]).prob.mean()
full = pd.read_parquet(ROOT / "results/pathways/oof/PB3.parquet")
full = full[(full.level == "genus") & full.genome_id.isin(labS.genome_id)].groupby(["trait", "genome_id"]).prob.mean()
ids, g = labS.genome_id.to_numpy(), labS.group_genus.to_numpy()
ug = np.unique(g)
rng = np.random.default_rng(SEED)
boots = [np.concatenate([np.where(g == u)[0] for u in rng.choice(ug, len(ug))]) for _ in range(B)]
rows = []
for t, j in zip(use, J):
    a, b = oofS.loc[t].reindex(ids).to_numpy(), full.loc[t].reindex(ids).to_numpy()
    y = YS[:, j]
    d = [safe_ap(y[i], a[i]) - safe_ap(y[i], b[i]) for i in boots if 0 < y[i].sum() < len(i)]
    rows.append(dict(trait=t, collection_only=safe_ap(y, a), final=safe_ap(y, b), diff=safe_ap(y, a) - safe_ap(y, b),
                     ci_low=np.nanpercentile(d, 2.5), ci_high=np.nanpercentile(d, 97.5)))
out = pd.DataFrame(rows)
out.to_csv(ROOT / "results/additional/pathway_hackathon_only/per_trait_ci.csv", index=False, float_format="%.4f")
print(out.round(3).to_string(index=False))
