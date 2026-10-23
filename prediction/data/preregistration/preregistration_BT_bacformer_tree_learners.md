# Pre-registration: BT, tree, tabular and AutoML learners on the Bacformer-large genome vector

Registered 19 October 2026 before any of these models was fitted. SHA-256 in
`docs/preregistration_BT_bacformer_tree_learners.sha256`. Descriptive sensitivity analysis; released models are not
changed and nothing in the paper's model choice depends on the outcome.

**Question.** The paper classifies the Bacformer-large genome vector with logistic regression (arm `E_large`) and
reports that a multilayer perceptron, a trained head and fine-tuning did no better. Do the learners that were run on
the KO panel (random forest, extra-trees, XGBoost, TabPFN v2, FLAML model search) do better on the Bacformer vector?

**Known before registration.** Only the existing Bacformer-large results (Supplementary Table S5): logistic regression
0.577 [0.534, 0.623] genus, 0.531 family; MLP 0.569 / 0.519; trained head 0.556 / 0.532; fine-tuned 0.549 / 0.535.
None of the five learners below has been run on this vector.

**Design.** Donor/acceptor task, 358 labeled genomes, 14 categories, features = the 960-dimensional Bacformer-large
mean vector (`data/processed/genome_features_bacformer_large_mean.parquet`), constant columns dropped.
Everything else identical to the KO-panel runs of the same learners:
* frozen genus- and family-blocked outer folds, inner grouped CV for configuration and MCC threshold
  (`redox_common.nested_cv`);
* RF and ET: `47_trees_xai.py` configs (500 trees, `min_samples_leaf` 1 or 5, balanced subsample, 3 seeds);
* XGBoost: the RedoxKO grid (`max_depth` 2 or 4, 300 trees, learning rate 0.05, 3 seeds, `scale_pos_weight`);
* TabPFN v2: `65_tabpfn.py` (4 estimators, MPS, pretraining limits ignored because 960 > 500 features);
* FLAML: `63_automl.py` (30 s budget per search, inner and outer alike).
Arms are named `RF_BF`, `ET_BF`, `XGB_BF`, `TPFN_BF`, `AM_BF`; outputs in `results/trees_xai/da_bf/oof/`.

**Endpoints.** Primary: macro AUPRC (14 categories) in genus-level CV, paired difference learner − `E_large` with a
genus-cluster bootstrap (`redox_common.paired_macro_delta`, 1,000 resamples, seed 20260917). Secondary: the same at
family level (family-cluster bootstrap). Five comparisons per level; reported with and without Bonferroni adjustment of
the interval (97.5% → 99.5% two-sided). No learner can replace a released model.

**Expected.** No learner clearly better than logistic regression (all paired intervals include zero or lie below it);
trees probably slightly worse on a dense embedding (probability about 70%). A learner "clearly better" means a paired
genus-level interval entirely above zero.

**Reporting.** All five rows go into Supplementary Table S5 whatever the outcome, and the Methods sentence on the
Bacformer learners is updated to name them.
