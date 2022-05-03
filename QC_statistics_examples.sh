## MN4 SLURM SBATCH Array submission script
## With /gpfs/projects/bsc83/utils/wrappers/Marenostrum.sbatch_wrapper.v1.sh -h you display help

/gpfs/projects/bsc83/utils/wrappers/Marenostrum.sbatch_wrapper.v1.sh \
	-j QC_statistics_examples \
	-o /gpfs/projects/bsc83/Projects/scRNAseq/aripol1/sc-eqtlgen-consortium-pipeline/ongoing/wg1-qc_filtering/scripts \
	-r `pwd` \
	-n 1 \
	-u 48 \
	-q bsc_ls \
	-w "02:00:00" \
	-a True \
	-f /gpfs/projects/bsc83/Projects/scRNAseq/aripol1/sc-eqtlgen-consortium-pipeline/ongoing/wg1-qc_filtering/QC_statistics_examples.tab \
	-c 'Rscript /gpfs/projects/bsc83/Projects/scRNAseq/aripol1/sc-eqtlgen-consortium-pipeline/ongoing/wg1-qc_filtering/QC_statistics.R' \
