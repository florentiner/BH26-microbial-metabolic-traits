#!/bin/sh
# Links the released models and data into code/ in the layout the scripts expect (results/models, data/processed,
# data/pathways). Run once from the repository root:  sh code/setup_layout.sh
# Then put the KOfam profiles and ko_list.gz (release 25 August 2026) under code/data/raw/kofam/.
set -e
cd "$(dirname "$0")"
mkdir -p results/models/pathways data/processed data/pathways data/raw/kofam
ln -sfn ../../../models/RedoxKO results/models/RF
ln -sfn ../../../models/RedoxKO-LX results/models/B2
ln -sfn ../../../models/RedoxKO-ens-BF results/models/B2XLE_large
ln -sfn ../../../models/RELEASED_donor_acceptor.json results/models/RELEASED.json
ln -sfn ../../../../models/EnergyKO-M results/models/pathways/PB3
ln -sfn ../../../../models/EnergyKO-BF results/models/pathways/PE
ln -sfn ../../../../models/RELEASED_pathways.json results/models/pathways/RELEASED.json
D=../../../data
for f in labels.tsv labels_curated_D.tsv ko_matrix_panel.tsv gtdb_taxonomy.tsv genome_qc.tsv; do
  ln -sfn $D/donor_acceptor/$f data/processed/$f; done
ln -sfn $D/donor_acceptor/ko_panel_852.tsv data/processed/ko_panel.tsv
ln -sfn $D/donor_acceptor/bacformer_large_vectors.parquet data/processed/genome_features_bacformer_large_mean.parquet
ln -sfn $D/pathways/kegg_energy_modules.tsv data/processed/kegg_modules_energy.tsv
for f in ../data/donor_acceptor/splits/*; do ln -sfn $D/donor_acceptor/splits/$(basename $f) data/processed/$(basename $f); done
ln -sfn $D/pathways/labels_and_folds.tsv data/pathways/labels.tsv
for f in labels.sha256 species_accessions.tsv module_completeness.tsv; do ln -sfn $D/pathways/$f data/pathways/$f; done
ln -sfn $D/pathways/ko_matrix_full.tsv data/pathways/ko_matrix.tsv
ln -sfn ../data/donor_acceptor/ko_presence_586genomes.tsv ko_presence_586genomes.tsv
echo "layout ready; add KOfam profiles to code/data/raw/kofam/"
