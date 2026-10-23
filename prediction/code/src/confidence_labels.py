"""Add confidence labels to the predictions of the final model RedoxKO.

Reads the table written by src/predict.py (--arm RF, the default released model) and adds, for every genome and
category, a confidence label from class-conditional (Mondrian) split-conformal prediction: "confident yes",
"confident no", "uncertain" (both answers possible) or "no label" (neither answer reaches the threshold). Labels are
given only in the six categories that get yes/no calls under the release rule (reliable_traits in metadata.json); the
other categories get "score only", because a confident label there would contradict the release rule.
The per-category thresholds were calibrated on the scores that the 358 training genomes received in family-level
cross-validation, with a target coverage of 90% (src/149_conformal_rf.py, Supplementary Table S10), and are stored
with the model in results/models/RF/conformal_thresholds.json. src/predict.py itself is not changed.

usage: python src/confidence_labels.py predictions.csv [--out labeled.csv]
       python src/confidence_labels.py --check    (reproduces the coverage of Supplementary Table S10)
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
THRESHOLDS = ROOT / "results/models/RF/conformal_thresholds.json"


def load_thresholds():
    return {t: {int(c): v for c, v in q.items()} for t, q in json.loads(THRESHOLDS.read_text())["thresholds"].items()}


def label(trait, score, q):
    """Mondrian conformal set from the raw RedoxKO score: 1 is included when 1 - score <= q1, 0 when score <= q0."""
    if trait not in q:
        return ""
    yes = (1 - score) <= q[trait][1]
    no = score <= q[trait][0]
    if yes and no:
        return "uncertain"
    return "confident yes" if yes else "confident no" if no else "no label"


def add_labels(pred, q):
    rel = set(json.loads((ROOT / "results/models/RF/metadata.json").read_text())["reliable_traits"])
    pred = pred.copy()
    if "arm" in pred and set(pred.arm.unique()) - {"RF"}:
        raise SystemExit("confidence labels are calibrated for RedoxKO only (predict.py --arm RF)")
    pred["confidence_label"] = [label(t, s, q) if t in rel else "score only" for t, s in zip(pred.trait, pred.score)]
    return pred


def check():
    """Recompute the labels on the evaluation sets of Supplementary Table S10 (six call categories) and compare."""
    q = load_thresholds()
    ref = pd.read_csv(ROOT / "results/additional/rf_x3_conformal_released.csv").query("scope == 'six'").set_index("set")
    ext = pd.read_csv(ROOT / "results/external_validation_missing_rows_predictions.csv").query("arm == 'Random forest'")
    ext = ext.merge(pd.read_csv(ROOT / "results/external_validation_predictions.csv").query("arm == 'B2'")[
        ["genome_id", "trait", "y"]], on=["genome_id", "trait"]).dropna(subset=["score", "y"])
    new = pd.read_csv(ROOT / "results/additional/rf_paper/x1_calls_rf.csv").rename(
        columns={"category": "trait", "expected": "y", "score_RF": "score", "accession": "genome_id"})
    for name, d in (("External genomes", ext), ("New species", new)):
        d = add_labels(d[["genome_id", "trait", "y", "score"]], q)
        d = d[d.confidence_label != "score only"]
        covered = np.where(d.y == 1, d.confidence_label.isin(["confident yes", "uncertain"]),
                           d.confidence_label.isin(["confident no", "uncertain"]))
        conf = d[d.confidence_label.isin(["confident yes", "confident no"])]
        acc = ((conf.confidence_label == "confident yes").astype(int) == conf.y).mean()
        print(f"{name}: coverage {covered.mean():.4f} (Table S10 {ref.at[name, 'coverage']:.4f}), "
              f"confident {len(conf) / len(d):.4f} ({ref.at[name, 'confident']:.4f}), "
              f"correct if confident {acc:.4f} ({ref.at[name, 'correct_if_conf']:.4f})")
        assert abs(covered.mean() - ref.at[name, "coverage"]) < 1e-9, name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("predictions", nargs="?")
    ap.add_argument("--out")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.check:
        return check()
    pred = add_labels(pd.read_csv(a.predictions), load_thresholds())
    out = a.out or str(Path(a.predictions).with_suffix("")) + "_labeled.csv"
    pred.to_csv(out, index=False)
    print(f"written {out}")
    print(pred.pivot_table(index="trait", columns="confidence_label", values="genome_id", aggfunc="count", fill_value=0))


if __name__ == "__main__":
    main()
