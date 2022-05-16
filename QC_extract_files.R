#!/usr/bin/env Rscript

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--dataset"), action="store", default=NULL, type='character',
              help="Name of the dataset."),
  make_option(c("--level"), action="store", default=NULL, type='character',
              help="Level at which the QC filtering is peformed. By default, it will be performed at the whole dataset. You can chose to peform it at the level of a specific metadata variable (such as, lane/pool)."),
  make_option(c("--qc_mad"), action="store", default="nCount_RNA_lower_1_5.percent.mt_upper_1_5", type='character',
              help="Name of the dataset."),
  make_option(c("--in_dir"), action="store", default="QC_statistics", type='character',
              help="Input directory: QC_statistics.R output directory. By default, it is named QC_statistics.")
)
opt = parse_args(OptionParser(option_list=option_list))

#################################### Set Variables ####################################
# Dataset
dataset <- opt$dataset
print(paste0('Running QC_statistics.R on: ', dataset))

# Input directory
cwd <- getwd()
setwd(cwd)
in.dir <- paste0(opt$in_dir, '/', dataset, '/', opt$qc_mad, '/')

# Loading functions
functions_fn <- 'scripts/QC_functions.R'
print(paste0('Loading functions from: ',functions_fn))
source(functions_fn)

# Output directory
out.dir <- paste0(opt$in_dir, '.files/')
print(paste0('Creating output directory in: ',out.dir))
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
out_ss.dir <- paste0(out.dir, 'summary_statistics/')
if(!dir.exists(out_ss.dir)){dir.create(out_ss.dir, recursive = T)}
out_hp.dir <- paste0(out.dir, 'heatmap_plots/')
if(!dir.exists(out_hp.dir)){dir.create(out_hp.dir, recursive = T)}

#################################### Extract files ####################################
# Summary statistics
ss_1 <- paste0(in.dir, 'by_dataset/dataset/tag.rds')
ss_2 <- paste0(in.dir, 'by_dataset/', opt$level, '/tag.rds')
ss_3 <- paste0(in.dir, 'by_metadata/dataset/tag.rds')
ss_4 <- paste0(in.dir, 'by_metadata/', opt$level, '/tag.rds')
statistics.list <- list(ss_1, ss_2, ss_3, ss_4)

# Heatmap plots
hp_1 <- paste0(in.dir, 'by_metadata/dataset/cell/')
hp_1.list <- as.list(list.files(hp_1, pattern='*_Outlier_pct.label_n.cluster_md.pdf', full.names = T))
hp_2 <- paste0(in.dir, 'by_metadata/dataset/donor/', opt$level, '_Outlier_pct.label_n.pdf')
hp_3 <- paste0(in.dir, 'by_metadata/', opt$level, '/donor/', opt$level, '_Outlier_pct.label_n.pdf')
heatmap_plots.list <- c(hp_1.list, list(hp_2, hp_3))

# Extract files
files_list <- list(summary_statistics = statistics.list,
                   heatmap_plots = heatmap_plots.list)
outdirs_list <- list(summary_statistics = out_ss.dir,
                     heatmap_plots = out_hp.dir)

extract_files.res <- sapply(names(files_list), function(i) extract_files(type = i, 
                                                                         outdirs = outdirs_list,
                                                                         files = files_list), simplify = F)