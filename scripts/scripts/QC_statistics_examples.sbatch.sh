#!/bin/bash

#SBATCH --job-name=QC_statistics_examples
#SBATCH --output=/gpfs/projects/bsc83/Projects/scRNAseq/aripol1/sc-eqtlgen-consortium-pipeline/ongoing/wg1-qc_filtering/scripts/logs/outputs/QC_statistics_examples.%A_%a.out
#SBATCH --error=/gpfs/projects/bsc83/Projects/scRNAseq/aripol1/sc-eqtlgen-consortium-pipeline/ongoing/wg1-qc_filtering/scripts/logs/errors/QC_statistics_examples.%A_%a.err
#SBATCH --workdir=/gpfs/projects/bsc83/Projects/scRNAseq/aripol1/sc-eqtlgen-consortium-pipeline/ongoing/wg1-qc_filtering
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --tasks-per-node=1
#SBATCH --qos=bsc_ls
#SBATCH --time=02:00:00
#SBATCH --array=1-8

echo '----------------------'
date
echo '----------------------'
echo ' '

Rscript /gpfs/projects/bsc83/Projects/scRNAseq/aripol1/sc-eqtlgen-consortium-pipeline/ongoing/wg1-qc_filtering/QC_statistics.R \
`sed -n ${SLURM_ARRAY_TASK_ID}p $ARRAY_FILE | cut -f1` \
`sed -n ${SLURM_ARRAY_TASK_ID}p $ARRAY_FILE | cut -f2` \
`sed -n ${SLURM_ARRAY_TASK_ID}p $ARRAY_FILE | cut -f3` \
`sed -n ${SLURM_ARRAY_TASK_ID}p $ARRAY_FILE | cut -f4` \
`sed -n ${SLURM_ARRAY_TASK_ID}p $ARRAY_FILE | cut -f5` \
`sed -n ${SLURM_ARRAY_TASK_ID}p $ARRAY_FILE | cut -f6` \
`sed -n ${SLURM_ARRAY_TASK_ID}p $ARRAY_FILE | cut -f7` \
`sed -n ${SLURM_ARRAY_TASK_ID}p $ARRAY_FILE | cut -f8` \
`sed -n ${SLURM_ARRAY_TASK_ID}p $ARRAY_FILE | cut -f9` \
`sed -n ${SLURM_ARRAY_TASK_ID}p $ARRAY_FILE | cut -f10` \
`sed -n ${SLURM_ARRAY_TASK_ID}p $ARRAY_FILE | cut -f11` \
`sed -n ${SLURM_ARRAY_TASK_ID}p $ARRAY_FILE | cut -f12` \


echo ' '
echo '----------------------'
date
echo '----------------------'
echo ' '


