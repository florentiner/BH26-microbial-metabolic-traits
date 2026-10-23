"""One-command prediction with confidence labels.

Runs src/predict.py unchanged (its SHA-256 is frozen for the registered blind tests), with all arguments passed through,
and then adds the confidence label of src/confidence_labels.py to every row of the donor and acceptor table.
Labels are added only for RedoxKO (--arm RF, which is also the released default); for other arms the table is left
as predict.py wrote it.

usage: python src/predict_with_labels.py GCF_000008665.1 [GCF_...] [--fasta genome.fna] [--out predictions.csv] ...
"""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import confidence_labels  # noqa: E402


def option(args, name, default):
    return args[args.index(name) + 1] if name in args and args.index(name) + 1 < len(args) else default


def main():
    args = sys.argv[1:]
    subprocess.run([sys.executable, str(ROOT / "src/predict.py"), *args], check=True)
    arm = option(args, "--arm", "released")
    if arm == "released":
        arm = json.loads((ROOT / "results/models/RELEASED.json").read_text())["released_arm"]
    out = option(args, "--out", "predictions.csv")
    if arm != "RF":
        print(f"confidence labels are calibrated for RedoxKO (RF) only, not added for arm {arm}")
        return
    pred = confidence_labels.add_labels(pd.read_csv(out), confidence_labels.load_thresholds())
    pred.to_csv(out, index=False)
    print(f"added confidence labels to {out}")
    print(pred.pivot_table(index="trait", columns="confidence_label", values="genome_id", aggfunc="count", fill_value=0)
          .to_string())


if __name__ == "__main__":
    main()
