# Pre-registration — A+B: hybrid of the KO model with METABOLIC's curated HMMs

Written 2026-10-15 before any hybrid feature was computed for a training genome or any hybrid model fitted. Hash in
`docs/preregistration_additional_experiments.sha256`.

## Why, and the fairness constraint
On the primary 100 new species METABOLIC beat B2 overall, mainly on oxygen respiration and hydrogen oxidation
(`docs/RESULTS_additional.md`, X1). Because that observation came from the primary 100, claims about the hybrid rest
only on data not used to motivate it: **primary test = the 43 held-out blind new species (75 growth tests)**;
secondary = the external 562 (Madin labels). The primary 100 is reported only. Training data: the same 358 genomes,
the same closed-world labels, the same frozen genus/family folds. No test genome is used for any choice below.

## A. Feature hybrid (arm B2M)
* Features: the 852-KO panel (as B2) **plus** METABOLIC v4's per-HMM presence calls (its function-table HMMs:
  custom HMMs with METABOLIC's own score cut-offs and score types, KOfam HMMs with the KOfam 2026-08-25 cut-offs,
  and its motif / motif-pair checks), recomputed with pyhmmer for every genome (`src/103_metabolic_features.py`).
* Fidelity check before use: the recomputed presence calls are compared with real METABOLIC-G output for the 100
  primary genomes (worksheet 1). If agreement over all HMM × genome cells is < 0.97, the features are not used and
  this is reported.
* Model: identical procedure to B2 (`06_run_arms.py` generic path, logistic regression or XGBoost chosen by inner
  grouped CV, MCC-optimal threshold on inner out-of-fold predictions); final model refitted on all 358 genomes into
  `results/models_hybrid/B2M/` (released models untouched).

## B. Decision hybrid (B2M+M)
For each category, use B2M's call where B2 currently releases a call (oxygen, sulfate, sulfur reduction, sulfur
oxidation, methanogen/acetogen — the released set fixed before this study); for every other category use
METABOLIC-G's function call with the Addendum A mapping (real METABOLIC run on the 43 held-out genomes).

## Evaluation
* Internal: B2M vs B2 macro AUPRC, genus and family level, paired cluster bootstrap (as for every arm).
* Held-out 43: accuracy, MCC, sensitivity, specificity for B2, B2M, B2M+M and METABOLIC on the same 75 calls; paired
  cluster bootstrap over papers.
* External 562: B2M vs B2, macro AUROC / AUPRC / MCC, paired bootstrap over genera.

## Predictions
* H-P1: B2M internal genus-level macro AUPRC ≥ B2 (direction; CI reported).
* H-P2: on the held-out 43, B2M+M MCC ≥ METABOLIC MCC (direction; CI reported). With 75 calls only large
  differences can be resolved; a tie is the expected outcome.
* H-P3: on the held-out 43, B2M+M MCC > B2 MCC.
* H-P4: on the external 562, B2M macro AUROC is not lower than B2's by more than 0.02.
