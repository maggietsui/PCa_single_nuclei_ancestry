#!/bin/bash
#$ -cwd
#$ -l gpu_mem=32000M
#$ -l h_rt=10:00:00
#$ -o /wynton/scratch/mtsui/AA_integration/cellbender/logs
#$ -e /wynton/scratch/mtsui/AA_integration/cellbender/logs
#$ -t 1-17
set -e

module load Sali anaconda
source activate cellbender

module load cuda

export CUDA_VISIBLE_DEVICES=$SGE_GPU

# inputs are path to raw_feature_bc_matrix.h5
readarray -t INPUTFILES < ~/AA_Jam_single_nuc_analysis/scripts/cellbender_inputs_multiome_realign.txt
filename="${INPUTFILES[$SGE_TASK_ID - 1]}"

# Get sample ID from path
sample=$(echo $filename  | cut -d "/" -f6)

if [ $sample == '3___11T' ]; then
    sample="PCa_3_11T"
fi

filtered_barcodes=/wynton/scratch/hsong1/CZI_JAM_Realign/"$sample"_filtered_feature_bc_matrix/barcodes.tsv.gz
# get cellranger filtered cell count
expected_cells=$(zcat "$filtered_barcodes" | wc -l)

# Need to create folder for sample so checkpoint doesn't get
# overwritten
cd /wynton/scratch/mtsui/AA_integration/cellbender/with_expected_cellcount

mkdir "$sample" 

# Go into sample folder to save checkpoint here
cd "$sample"

cellbender remove-background \
	--cuda \
	--input $filename \
    --learning-rate 0.0000125 \
    --expected-cells $expected_cells \
	--output /wynton/scratch/mtsui/AA_integration/cellbender/with_expected_cellcount/"$sample"/"$sample"_cellbender.h5


# To load result into anndata
# from cellbender.remove_background.downstream import anndata_from_h5
# import anndata
# adata = anndata_from_h5("/wynton/group/fhuang/cellbender/PCa11N_PCa12N_PCa14N_FF/raw_feature_bc_matrix_cellbender.h5")