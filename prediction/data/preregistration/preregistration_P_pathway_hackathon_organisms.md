# Pre-registration: P, the pathway model trained only on the hackathon's organisms

Registered 18 October 2026 before the analysis was run. SHA-256 in `docs/preregistration_P_pathway_hackathon_organisms.sha256`.
Descriptive sensitivity analysis; released models are not changed.

**Question.** How much do pathway metrics change if EnergyKO-M is trained only on the organisms of the hackathon's
donor/acceptor collection (the 174 of 358 that have pathway labels) instead of all 3,378 genomes?

**Known before registration.** Only label counts: in the 174 genomes, five traits have at least 25 positives
(sulfur oxidation 40, nitrate reduction 39, methanogenesis 32, hydrogen oxidation 31, sulfur reduction 30).

**Design.** EnergyKO-M procedure unchanged (full KO profile + module completeness, logistic regression or XGBoost per
trait chosen by inner grouped CV, `src/32_pathway_arms.py` / `src/35_pathway_train_final.py`), fitted on the 174
genomes, for the traits with at least 25 positives among them.
* **Cross-validation:** the paper's frozen genus- and family-blocked pathway folds restricted to the 174 genomes.
  Compared with the released model's out-of-fold predictions on the **same 174 genomes and traits**
  (`results/pathways/oof/PB3.parquet`). Macro AUPRC; paired genus-cluster bootstrap, 2,000 resamples, seed 20261018.
* **New species:** the growth tests of the 99 pathway species for these traits; macro AUROC of the 174-genome model vs
  the released EnergyKO-M; paper bootstrap.
* **Expected:** similar or slightly lower point estimates with much wider intervals (about 60%).
