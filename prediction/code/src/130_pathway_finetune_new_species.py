"""Score the final fine-tuned Bacformer-large pathway model (arm FTP_BFL2) on the pathway new-species test.

Input: job_FTP_BFL2_final_-1.json from the Kaggle kernel your-kaggle-user/redox-pw-finetune-final
(kaggle/pathway_finetune_final/), downloaded to data/kaggle/pw_finetune_final/. Decision thresholds per trait are
the best-MCC thresholds on the genus-level cross-validation predictions of the same arm
(data/kaggle/pw_finetune/job_FTP_BFL2_genus_*.json), as for the other final models. Tests, metrics and bootstrap are
those of src/84_novel_pathway_comparison.py.

usage: python src/130_pathway_finetune_new_species.py
Output: validation_new_organisms/novel100_pathway_finetune.csv (+ adds score_FT/call_FT to
        novel100_pathway_extra_scores.csv)
"""
import glob
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pathway_common import load_pathways  # noqa: E402
from redox_common import best_mcc_threshold  # noqa: E402

V = ROOT / "validation_new_organisms"
cmp84 = importlib.import_module("84_novel_pathway_comparison")


def main():
    lab, traits, Y = load_pathways()
    oof = np.full(Y.shape, np.nan)
    for f in glob.glob(str(ROOT / "data/kaggle/pw_finetune/job_FTP_BFL2_genus_*.json")):
        j = json.load(open(f))
        oof[j["test_idx"]] = np.array(j["test_pred"])
    assert not np.isnan(oof).any(), "genus-level CV predictions incomplete"
    thr = {t: best_mcc_threshold(Y[:, k], oof[:, k]) for k, t in enumerate(traits)}
    fin = json.load(open(ROOT / "data/kaggle/pw_finetune_final/job_FTP_BFL2_final_-1.json"))
    S = pd.DataFrame(np.array(fin["ext_scores"]), index=fin["ext_ids"], columns=traits)
    d = pd.read_csv(V / "novel100_pathway_extra_scores.csv")
    missing = sorted(set(d.accession) - set(S.index))
    print(f"final model scored {len(S)} genomes; tests without a score: {d.accession.isin(missing).sum()} "
          f"({len(missing)} genomes {missing})")
    d["score_FT"] = [S.at[a, t] if a in S.index else np.nan for a, t in zip(d.accession, d.trait)]
    d["call_FT"] = [int(s >= thr[t]) if pd.notna(s) else np.nan for s, t in zip(d.score_FT, d.trait)]
    d.to_csv(V / "novel100_pathway_extra_scores.csv", index=False)
    rows = []
    for sub, dd in [("all tests", d), ("tests with an FT score", d[d.score_FT.notna()])]:
        for a in ["PB3", "FT", "ET", "tax", "rule"]:
            x = dd[dd[f"score_{a}"].notna()]
            au, ap, nt = cmp84.ranking(x, f"score_{a}")
            lo, hi = cmp84.boot(x, lambda z: cmp84.ranking(z, f"score_{a}")[0])
            rows.append(dict(subset=sub, model=a, tests=len(x), auroc=au, auroc_lo=lo, auroc_hi=hi, auprc=ap,
                             traits=nt, mcc=cmp84.mcc(x.expected, x[f"call_{a}"].astype(int))))
    t = pd.DataFrame(rows)
    t.to_csv(V / "novel100_pathway_finetune.csv", index=False, float_format="%.4f")
    print(t.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
