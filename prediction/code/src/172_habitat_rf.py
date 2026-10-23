"""X4 habitat case study repeated with the final model RedoxKO (random forest, src/predict.py --arm RF).

Same 300 MAGs (data/additional/habitat/selection.tsv), same pre-registered expectations X4-E1..E5, same tests and
BH correction as src/90_habitat_case.py (registered with B2, now RedoxKO-LX); only the model changes. Calls are
RedoxKO's released calls; categories without a released call use RedoxKO's thresholded score, as the original did.
usage: python src/172_habitat_rf.py download|predict|analyze
Output: data/additional/habitat/pred_rf/, results/additional/rf_habitat/x4_predicted_fraction_by_habitat.csv,
        x4_tests.csv
"""
import importlib
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
h = importlib.import_module("90_habitat_case")
PRED = h.WORK / "pred_rf"
h.OUT = ROOT / "results/additional/rf_habitat"


def predict(batch=20):
    PRED.mkdir(parents=True, exist_ok=True)
    st = pd.read_csv(h.WORK / "download_status.tsv", sep="\t")
    accs = sorted(st[st.status == "ok"].accession)
    done = {p.stem[2:] for p in PRED.glob("b_*.csv")}
    todo = [a for a in accs if a not in done]
    for i in range(0, len(todo), batch):
        chunk = todo[i:i + batch]
        out = PRED / f"batch_{i:05d}_{chunk[0]}.csv"
        cmd = [sys.executable, str(ROOT / "src/predict.py"), "--arm", "RF", "--pathways", "none", "--out", str(out)]
        for a in chunk:
            cmd += ["--fasta", str(h.FNA / f"{a}.fna")]
        t0 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print("FAILED", chunk[0], r.stderr[-800:], flush=True)
            continue
        for a in chunk:
            (PRED / f"b_{a}.csv").write_text("done\n")
        print(f"{i + len(chunk)}/{len(todo)} in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    h.PRED = PRED
    h.OUT.mkdir(parents=True, exist_ok=True)
    {"download": h.download, "predict": predict, "analyze": h.analyze}[sys.argv[1]]()
