# Addendum B — DRAM category mapping (fixed before DRAM distill was run)

Written 2026-10-14 ~23:30 MSK after reading DRAM's distill forms (github.com/WrightonLabCSU/DRAM, master, data/*.tsv)
and before any genome was distilled. Input: KOfam annotations of the 100 primary new species (KOfam release
2026-08-25, per-profile thresholds as in DRAM's KOfam mode) for every KO referenced by DRAM's forms, distilled with
DRAM's own `distill` code. A category is predicted positive when any listed DRAM function is present in DRAM's
product output; KEGG-module (ETC) entries count as present at coverage ≥ 0.5 (DRAM reports coverage, not a call).

| Our category | DRAM product function |
|---|---|
| donor_sulfur_oxidation | Thiosulfate oxidation by SOX complex, thiosulfate => sulfate |
| donor_ammonia_oxidation | Bacterial/Archaeal ammonia oxidation |
| donor_nitrite_oxidation | nitrite=> nitrate |
| donor_hydrogen_oxidation | not covered by DRAM's distillate (scored negative, reported as not covered) |
| donor_CO_oxidation | not covered |
| donor_iron_oxidation | not covered |
| acceptor_oxygen | ETC Complex IV: M00154, M00155, M00416, M00417 (low affinity) or M00156 (cbb3); bd oxidase M00153 excluded |
| acceptor_nitrate | nitrate => nitrite |
| acceptor_sulfate | dissimilatory sulfate reduction (and oxidation) sulfate => sulfide |
| acceptor_sulfur | thiosulfate => sulfite; tetrathionate => thiosulfate |
| acceptor_CO2_methanogen_acetogen | Methanogenesis "Key functional gene" (mcrABG); acetogenesis not covered |
| acceptor_iron | not covered |

Predictions as in Addendum A (A-P1, A-P2) apply to DRAM as well.
