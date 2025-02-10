#!/bin/bash -l
#SBATCH --job-name=mmvec_V1V3_V4_merged
#SBATCH --output=../logs/mmvec_%j.out
#SBATCH --error=../logs/mmvec_%j.err
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=short

# source activate mmvec_qiime2-2020.6_barnacle2
# to edit mmvec source code: /home/yac027/mambaforge3/envs/mmvec_qiime2-2020.6_barnacle2/lib/python3.6/site-packages/mmvec/

echo '*** STEP 1: Converting biom files to qza ***'
qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path ../tables/16S_MERGED-Avg-tbl_metaB-matched_relative_Healthy.biom \
  --output-path ../tables/16S_MERGED-Avg-tbl_metaB-matched_relative_Healthy.qza
qiime tools import \
  --type 'FeatureTable[Frequency]' \
  --input-path ../tables/metaB_MERGED-Avg-tbl_matched_relative_Healthy.biom \
  --output-path ../tables/metaB_MERGED-Avg-tbl_matched_relative_Healthy.qza


echo '*** STEP 2: Checking qza files ***'
qiime tools peek ../tables/16S_MERGED-Avg-tbl_metaB-matched_relative_Healthy.qza | tee ../logs/16S_MERGED-Avg-tbl_metaB-matched_peek_relative_Healthy.log
qiime tools peek ../tables/metaB_MERGED-Avg-tbl_matched_relative_Healthy.qza | tee ../logs/metabolomics-tbl_MERGED-matched_peek_relative_Healthy.log


echo '*** STEP 3: Running mmvec paired omics ***'
qiime mmvec paired-omics \
  --i-microbes ../tables/16S_MERGED-Avg-tbl_metaB-matched_relative_Healthy.qza \
  --i-metabolites ../tables/metaB_MERGED-Avg-tbl_matched_relative_Healthy.qza \
  --p-learning-rate 1e-3 \
  --p-epochs 500 \
  --output-dir ../outputs/Healthy/


echo '*** STEP 4: Checking mmvec outputs ***'
qiime tools peek ../outputs/Healthy/conditionals.qza | tee ../logs/MERGED-metaB_conditionals_peek_relative_Healthy.log
qiime tools peek ../outputs/Healthy/conditional_biplot.qza | tee ../logs/MERGED-metaB_conditional_biplot_peek_relative_Healthy.log


echo '*** STEP 5: Running qiime metadata tabulate ***'
qiime metadata tabulate \
        --m-input-file ../outputs/Healthy/conditionals.qza \
        --o-visualization ../outputs/Healthy/conditionals-viz.qzv 2>&1 | tee ../logs/MERGED-metaB_conditionals-viz_relative_Healthy.log

echo '*** STEP 6: Generate model convergence summaries ***'
qiime mmvec summarize-single \
        --i-model-stats ../outputs/Healthy/model_stats.qza \
        --o-visualization ../outputs/Healthy/model-summary.qzv


echo '*** STEP 7: Generating ranks file ***'
qiime mmvec heatmap \
  --i-ranks ../outputs/Healthy/conditionals.qza \
  --o-visualization ../outputs/Healthy/heatmap-ranks.qzv \
  --p-color-palette "RdBu_r" \
  --p-x-labels True \
  --p-y-labels True \
  --p-level 50 | tee ../logs/MERGED_heatmap-ranks_relative_Healthy.log


echo '*** STEP 8: Running emperor biplot ***'
qiime emperor biplot \
    --i-biplot ../outputs/Healthy/conditional_biplot.qza \
    --m-sample-metadata-file ../metadata/metabolite_info.csv \
    --o-visualization ../outputs/Healthy/emperor-biplot.qzv \
    --p-number-of-features 10 2>&1 | tee ../logs/MERGED_emperor-biplot_relative_Healthy.log


echo '*** STEP 9: Export MMVEC conditionals output to tsv ***'
qiime tools export \
  --input-path ../outputs/Healthy/conditionals.qza \
  --output-path ../outputs/Healthy/exported_conditionals

echo '*** Finished ***'
rm -rf log*
