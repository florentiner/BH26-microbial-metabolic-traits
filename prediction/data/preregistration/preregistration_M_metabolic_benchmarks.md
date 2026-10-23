# Pre-registration: M, where RedoxKO and METABOLIC each win

Registered 17 October 2026, before any METABOLIC call was compared with labels or with RedoxKO on the external
genomes or on MAGs. SHA-256 in `docs/preregistration_M_metabolic_benchmarks.sha256`. Released models are not changed.

## Question
On the 255 new-species tests METABOLIC gave better yes/no calls (MCC), whereas RedoxKO ranked better and was more
specific. Does this pattern hold on a larger, independent set, and do the two tools differ in robustness to genome
incompleteness and in speed? All results are reported as a scorecard of wins and losses, whatever their direction.

## What was looked at before registration
* All new-species comparisons in the paper (primary 100, withheld 43), including METABOLIC's sensitivity (0.80),
  specificity (0.72), precision (0.81) and macro AUROC (0.64) on the primary set, against RedoxKO's 0.52, 0.84, 0.83
  and 0.73.
* METABOLIC calls from the pyhmmer re-implementation (`src/103_metabolic_features.py`) mapped to categories by
  `src/metabolic_calls.py`. The mapping rule (any row of a mapped function; sulfite reduction also counts *dsrA*/*dsrB*
  and sulfur oxidation does not, as METABOLIC's worksheet 2 does) was fixed by comparison with real METABOLIC output
  only: 100% agreement on 143 genomes x 12 categories. No label was used for this.
* RedoxKO's external scores (`results/external_validation_predictions.csv`, arm B2), already published in the paper.
* Not looked at: METABOLIC calls against external labels, and METABOLIC calls on any MAG or paired isolate.

## Tools and data
* **RedoxKO:** released model; raw scores; calls at its own thresholds for every category (as in the tool
  comparison). Scores on MAGs and isolates: `results/additional/x2_pair_calls.csv` (already computed).
* **METABOLIC:** category calls (0/1) as above; as a ranking score, the 0/1 call.
* **External set:** the 562 genomes of unseen genera, 9 categories (oxygen on 466 genomes with labels), labels from
  the Madin compilation.
* **MAG pairs:** the MAG-isolate pairs of `results/additional/x2_pair_calls.csv`. MAG and isolate nucleotide FASTA are
  downloaded again from NCBI (as in `src/88_real_mags.py`), proteins are called with the same Pyrodigal path
  (`src/predict.py proteins_from_fasta`) for both, and METABOLIC rows are scanned with `src/103` `scan_genome`.

## Endpoints
Bootstrap: 2,000 resamples, seed 20261017; external clustered by genus, MAG pairs by species, new species by paper.
Macro = mean over categories with both classes.
* **M1 (primary, external):** macro AUROC, RedoxKO minus METABOLIC. Success: lower 95% bound above 0.
* **M1b:** macro MCC (calls), RedoxKO minus METABOLIC, external; all 9 categories, and separately RedoxKO's 5 released
  categories (of which 4 are external categories). Reported two-sided.
* **M1c:** pooled sensitivity, specificity and precision of calls for both tools, external.
* **M2 (scorecard, descriptive):** per category, on the external set and on the 100 primary new species: AUROC, MCC,
  sensitivity and specificity for both tools, with bootstrap intervals of the differences. A category is a "win" for
  a tool when the 95% interval of the difference excludes 0 in its favor. BH-adjusted q values are also reported.
* **M3 (MAG robustness):** agreement of category calls between MAG and isolate, by completeness bin, for both tools on
  the 12 METABOLIC-mapped categories. Primary: 50-70% bin, agreement RedoxKO minus METABOLIC; success: lower 95% bound
  above 0. Secondary: fraction of isolate-positive calls lost in the MAG, per tool and bin.
* **M4 (speed, descriptive):** wall time per genome. METABOLIC-G v4: from the logged run on the 100 primary genomes
  (1 h 39 min 59 s, 6 threads, `data/additional/metabolic/full_run.log`). RedoxKO: the panel scan and full prediction
  times reported in the paper, measured on the same workstation.
* **Expected results, stated in advance:** M1 likely positive (about 60-70%), since a 0/1 call has limited AUROC;
  M1b likely favors METABOLIC or shows no difference; M3 uncertain (about 50-60% that RedoxKO is more robust); M4
  favors RedoxKO.

## Reporting
All endpoints are reported in `docs/RESULTS_Bclean_A_D.md` (section M) and the supplement. The paper gains one
"where each tool wins" table if the results are informative in either direction.
