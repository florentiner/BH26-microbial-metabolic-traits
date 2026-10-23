"""Analysis P (docs/preregistration_P_pathway_hackathon_organisms.md): EnergyKO-M trained only on the hackathon's
organisms (the 174 donor/acceptor training genomes that have pathway labels), compared with the released EnergyKO-M on
the same genomes and traits (genus- and family-blocked CV) and on the new-species growth tests.
Output: results/additional/pathway_hackathon_only/results.json"""
import importlib
import json
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
tf = importlib.import_module("35_pathway_train_final")
nov = importlib.import_module("84_novel_pathway_comparison")
OUT = ROOT / "results/additional/pathway_hackathon_only"
OUT.mkdir(parents=True, exist_ok=True)
SEED, B = 20261018, 2000

lab, traits, Y = load_pathways()
da, _, _ = load_labelled()
keep = lab.genome_id.isin(set(da.genome_id)).to_numpy()
labS, YS = lab[keep].reset_index(drop=True), Y[keep]
use = [t for j, t in enumerate(traits) if YS[:, j].sum() >= 25]
J = [traits.index(t) for t in use]
print(f"hackathon organisms with pathway labels: {len(labS)}; traits with >=25 positives: {use}", flush=True)
Xdf = pd.concat([ko_matrix(labS.genome_id), module_completeness(labS.genome_id)], axis=1)
Xdf = Xdf.loc[:, Xdf.std(axis=0) > 0]
X = Xdf.to_numpy(dtype=np.float32)

# ---------------- CV on the 174 genomes (paper's frozen folds restricted to them)
outs = Parallel(n_jobs=7)(delayed(nested_cv)("PB3_hack", X, labS, [t], YS[:, [j]], arms.configs(), arms.fit_predict, lvl)
                          for lvl in ("genus", "family") for t, j in zip(use, J))
oofS = pd.concat(outs).groupby(["level", "trait", "genome_id"]).prob.mean()
full = pd.read_parquet(ROOT / "results/pathways/oof/PB3.parquet")
full = full[full.genome_id.isin(labS.genome_id)].groupby(["level", "trait", "genome_id"]).prob.mean()
res = {"n_genomes": int(len(labS)), "traits": use, "positives": {t: int(YS[:, j].sum()) for t, j in zip(use, J)}}
for lvl in ("genus", "family"):
    g = labS[f"group_{lvl}"].to_numpy()
    ids = labS.genome_id.to_numpy()

    def macro(idx, src):
        v = [safe_ap(YS[idx, j], src.loc[(lvl, t)].reindex(ids[idx]).to_numpy()) for t, j in zip(use, J) if YS[idx, j].sum() > 0]
        return float(np.nanmean(v))
    allidx = np.arange(len(ids))
    a, b = macro(allidx, oofS), macro(allidx, full)
    rng = np.random.default_rng(SEED)
    ug = np.unique(g)
    diffs = []
    for _ in range(B):
        idx = np.concatenate([np.where(g == u)[0] for u in rng.choice(ug, len(ug))])
        diffs.append(macro(idx, oofS) - macro(idx, full))
    per = {t: {"hackathon_only": safe_ap(YS[:, j], oofS.loc[(lvl, t)].reindex(ids).to_numpy()),
               "full": safe_ap(YS[:, j], full.loc[(lvl, t)].reindex(ids).to_numpy())} for t, j in zip(use, J)}
    res[f"cv_{lvl}"] = {"macro_auprc_hackathon_only": a, "macro_auprc_full_model_same_genomes": b, "diff": a - b,
                        "diff_ci": [float(x) for x in np.percentile(diffs, [2.5, 97.5])], "per_trait": per}
# full model on all 3,378 genomes for the same traits (context)
fa = pd.read_parquet(ROOT / "results/pathways/oof/PB3.parquet")
fa = fa[fa.level == "genus"].groupby(["trait", "genome_id"]).prob.mean()
res["cv_genus_full_model_all_3378_same_traits"] = float(np.mean(
    [safe_ap(Y[:, j], fa.loc[t].reindex(lab.genome_id).to_numpy()) for t, j in zip(use, J)]))

# ---------------- new species: final model on the 174 genomes
bundle = {}
for t, j in zip(use, J):
    m, thr, cfg, inner, cal = tf.select_and_fit(Xdf, YS[:, j], labS.group_genus.to_numpy())
    bundle[t] = (m, thr)
accs = open(ROOT / "validation_new_organisms/balanced_pathway_100.txt").read().split()
recs = json.loads((ROOT / "validation_new_organisms/validation_set.json").read_text())
Xv = nov.features(accs)
d = nov.calls(recs, set(accs))
d = d[d.trait.isin(use)].reset_index(drop=True)
import joblib  # noqa: E402
rel = joblib.load(ROOT / "results/models/pathways/PB3/model.joblib")
for name, get in (("hack", lambda t: bundle[t][0]), ("full", lambda t: rel[t]["model"])):
    sc = {}
    for t in use:
        m = get(t)
        sc[t] = pd.Series(np.asarray(m.predict_proba(Xv.reindex(columns=m.columns).fillna(0.0))), index=Xv.index)
    d[f"score_{name}"] = [sc[t][a] for a, t in zip(d.accession, d.trait)]
a = nov.ranking(d, "score_hack")[0]
b = nov.ranking(d, "score_full")[0]
gp = {p: x for p, x in d.groupby("paper")}
keys = np.array(list(gp))
rng = np.random.default_rng(SEED)
bs = []
for _ in range(B):
    s = pd.concat([gp[k] for k in rng.choice(keys, len(keys))])
    try:
        bs.append(nov.ranking(s, "score_hack")[0] - nov.ranking(s, "score_full")[0])
    except Exception:
        pass
res["new_species"] = {"tests": int(len(d)), "traits_with_both_classes": nov.ranking(d, "score_full")[2],
                      "macro_auroc_hackathon_only": a, "macro_auroc_full": b, "diff": a - b,
                      "diff_ci": [float(x) for x in np.nanpercentile(bs, [2.5, 97.5])]}
(OUT / "results.json").write_text(json.dumps(res, indent=2, default=float))
print(json.dumps(res, indent=1, default=lambda x: round(float(x), 3)))
