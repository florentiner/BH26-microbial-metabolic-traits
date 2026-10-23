# Pre-registration — electron donor / acceptor prediction from genomes

Written 2026-09-17, after labels and splits were frozen (`data/processed/splits.sha256` =
`ae0cdef16851fbf003fc0a4b6c1dbb30e9a5907745b587470dacfd8407ff1146`) and **before any model
was fitted or any prediction was scored**. Adapted from `redox-model-search-protocol.html`.
The SHA-256 of this file is recorded in `results/run_manifest.json` by the evaluation script.

## Task
Input: RefSeq genome accession → proteome (NCBI Datasets) → features.
Output: 14 binary traits (7 electron-donor categories, 7 electron-acceptor categories),
formulated as **multi-label classification** (independent sigmoid per trait; not seq2seq —
the label vocabulary is a fixed set of 14 categories and N = 358 labelled genomes is far too small
to learn a generative decoder).

## Data
* 358 labelled genomes (database entries); 228 "not in database" genomes are unlabelled and used
  only as the application set.
* 0 = negative under a closed-world assumption (limitation).
* All 14 traits have ≥10 positives and ≥10 negatives (nitrite oxidation: 12 positives — interpret with care).

## Evaluation design
* Outer CV: 5 frozen folds, grouped by **genus** (primary) and by **family** (stricter).
* Hyper-parameters and decision thresholds: selected on an inner grouped CV (same grouping level)
  inside the outer training folds only. Test folds are never used for selection.
* Metrics per trait × arm × split level: AUPRC (primary), AUPRC lift over prevalence, AUROC,
  MCC and F1 at the inner-CV threshold, Brier score. 95% CIs by **cluster bootstrap over
  groups** (1000 resamples) of pooled out-of-fold predictions.
* Paired comparisons between arms: cluster-bootstrap distribution of ΔAUPRC; macro-average over traits.
* Headline number: `delta_vs_null` = AUPRC(arm) − AUPRC(phylogenetic nearest-neighbour null).
* Seeds: 5 for neural arms; 3 for stochastic tree arms.

## Nulls
* **Prevalence null**: constant training prevalence (AUPRC ≈ prevalence).
* **Taxonomic nearest-neighbour null**: prediction = mean label of training genomes sharing the
  lowest taxonomic rank with the test genome (genus → family → order → class → phylum → domain →
  global prevalence). No phylogenetic tree is available, so taxonomy stands in for phylogeny.

## Arms (fixed budget, declared now)
| Arm | Features | Model | Budget |
|---|---|---|---|
| A  | 84 marker KOs from recomputed KOfam annotation | expert rules: category positive if any of its marker complexes reaches module fraction ≥ θ; θ ∈ {0.5, 1.0} selected on inner CV | 2 configs |
| B1 | 84 marker KOs (as provided in `ko_presence_586genomes.tsv`) | L2 logistic regression (C grid 4) and gradient boosting (depth {2,4}, 3 seeds) | 6 configs |
| B2 | 852-KO energy-metabolism panel (KOfam, recomputed from proteome) | same as B1 | 6 configs |
| C  | ESM-2 t12 35M protein embeddings, mean- and max-pooled per genome | L2 logistic probe on standardised features (C grid 4) | 8 configs |
| D  | same protein embeddings, set of proteins | gated-attention MIL (ABMIL), multi-trait, masked BCE | ≤12 configs, 5 seeds (only if time permits) |
| E  | Bacformer-base (masked-complete-genomes), frozen, genome-coordinate-ordered proteins, mean-pooled | L2 logistic probe (C grid 4) | 4 configs |

| B1r | 84 marker KOs **recomputed** from the proteome (subset of the B2 matrix) | same as B1 | 6 configs |

B1 uses the table provided with the labels; B1r/B2/A use only features derivable from the genome
sequence, which is what the released predictor must consume. B1 vs B1r quantifies annotation
reproducibility. Rule groups for arm A follow the marker-gene column blocks of the provided table
(mapping in `src/redox_common.py::RULE_GROUPS`); traits with no marker genes (iron oxidation, iron
reduction) receive a constant rule score. Group score = fraction of the group's KOs present; trait
score = max over its groups; θ ∈ {any KO (>0), 0.5, 1.0}.

Deviations from the protocol declared up front:
* Bacformer-**large** / ESM-C 300M were benchmarked at 15 proteins/s on the local Apple M3 Pro
  (≈34 h for 1.8 M proteins); Bacformer-**base** with its own encoder ESM-2 35M is used instead.
  C and E therefore share the encoder (encoder-matching rule is respected).
* Arm F (fine-tuning) is run only if E beats C on held-out genus (stopping rule); arm G (retrieval
  seed sets) and arm H (LLM over KOs through MCP) are out of scope for this model-training work.
* No phylogenetic D statistic (no tree); taxonomic-NN null reported instead.
* Full KofamScan (27,866 KOs) is not tractable locally; the energy-metabolism KO panel is used.

## Model selection for the released model
The released model for each trait is the arm with the highest genus-level macro AUPRC across
traits (single arm for all traits, to avoid per-trait winner's-curse selection); its
configuration is re-selected by grouped CV on all 358 labelled genomes and refitted on all of them.
Its reported performance is the outer-CV estimate — not a refit score.

## Predictions (recorded before running)
* P1. Arm B (KO + model) ties or beats every embedding arm (C, E) on macro AUPRC at genus level.
* P2. Arm B2 (energy panel) ≥ B1 (84 markers) — more KOs help traits whose markers are ambiguous
  (nitrite oxidation nxrAB = narGH KOs; iron oxidation/reduction have no KO markers).
* P3. Iron oxidation, iron reduction, "other" donors and "other" acceptors are the worst traits for
  every arm (no dedicated KO markers).
* P4. Arm E ≈ arm C (synteny adds little for category-level traits at N = 358).
* P5. On held-out family, the taxonomic-NN null is beaten by learned arms on the marker-rich traits
  (ammonia oxidation, sulfate reduction, methanogenesis/acetogenesis, oxygen) but not on iron or
  "other" traits.
* P6. Performance drops from genus- to family-level hold-out for every learned arm.
