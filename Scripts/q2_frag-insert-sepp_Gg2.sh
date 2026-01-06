#!/usr/bin/env bash

# CONFIGS
QIIME_ENV="qiime2-metagenome-2024.10"
V1V3_FASTA="../Data/16S/Fasta/179426_V1V3_ASVs.fasta"
V4_FASTA="../Data/16S/Fasta/174951_V4_ASVs.fasta"
OUT_DIR="../Analyses/16S/SEPP"
GG2_SEPP_REF="../Taxonomy/Greengenes2/2022.10.backbone.sepp-reference.qza"
THREADS=16

# ACTIVATE CONDA + QIIME 2
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${QIIME_ENV}"

echo "Activated environment: ${QIIME_ENV}"
which qiime
qiime --version

mkdir -p "${OUT_DIR}"

# IMPORT FASTA → QZA
echo "Importing V1–V3 ASVs..."
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "${V1V3_FASTA}" \
  --output-path "${OUT_DIR}/179426_V1V3_ASVs.qza"

echo "Importing V4 ASVs..."
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path "${V4_FASTA}" \
  --output-path "${OUT_DIR}/174951_V4_ASVs.qza"

# SEPP FRAGMENT INSERTION
echo "Running SEPP for V1–V3..."
qiime fragment-insertion sepp \
  --i-representative-sequences "${OUT_DIR}/179426_V1V3_ASVs.qza" \
  --i-reference-database "${GG2_SEPP_REF}" \
  --o-tree "${OUT_DIR}/179426_V1V3_sepp_tree.qza" \
  --o-placements "${OUT_DIR}/179426_V1V3_sepp_placements.qza" \
  --p-threads "${THREADS}"

echo "Running SEPP for V4..."
qiime fragment-insertion sepp \
  --i-representative-sequences "${OUT_DIR}/174951_V4_ASVs.qza" \
  --i-reference-database "${GG2_SEPP_REF}" \
  --o-tree "${OUT_DIR}/174951_V4_sepp_tree.qza" \
  --o-placements "${OUT_DIR}/174951_V4_sepp_placements.qza" \
  --p-threads "${THREADS}"

# GREENGENES2 TAXONOMY CLASSIFICATION (GENUS-LEVEL)
GG2_CLASSIFIER="../Taxonomy/Greengenes2/2022.10.backbone.full-length.nb.sklearn-1.4.2.qza"

echo "Running Greengenes2 taxonomy classification for V1–V3..."
qiime feature-classifier classify-sklearn \
  --i-classifier "${GG2_CLASSIFIER}" \
  --i-reads "${OUT_DIR}/179426_V1V3_ASVs.qza" \
  --o-classification "${OUT_DIR}/179426_V1V3_GG2_taxonomy.qza"

echo "Running Greengenes2 taxonomy classification for V4..."
qiime feature-classifier classify-sklearn \
  --i-classifier "${GG2_CLASSIFIER}" \
  --i-reads "${OUT_DIR}/174951_V4_ASVs.qza" \
  --o-classification "${OUT_DIR}/174951_V4_GG2_taxonomy.qza"


echo "SEPP + Greengenes2 workflow completed successfully."
