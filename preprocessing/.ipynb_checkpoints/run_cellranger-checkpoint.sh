#!/bin/bash
#$ -cwd
#$ -pe smp 8
#$ -l mem_free=6G
#$ -l h_rt=200:00:00
#$ -o /wynton/scratch/mtsui/AA_integration/logs
#$ -e /wynton/scratch/mtsui/AA_integration/logs
#$ -t 1-1

set -e

module load CBI
module load cellranger/7.0.0

readarray -t INPUTFILES < ~/AA_Jam_single_nuc_analysis/scripts/cellranger_fastqs.txt
fastq_dir="${INPUTFILES[$SGE_TASK_ID - 1]}"

# get sample name
sample=$(echo "$fastq_dir" | awk -F "/" '{print $NF}')
cellranger count \
    --nosecondary \
    --disable-ui \
    --id "$sample" \
    --fastqs "$fastq_dir" \
    --sample "$sample" \
    --transcriptome /wynton/group/fhuang/Cellranger/hg38 \
    --jobmode /wynton/group/fhuang/Cellranger/sge.template \
    --localcores 8

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"