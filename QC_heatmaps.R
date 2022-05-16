#!/usr/bin/env Rscript

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--datasets"), action="store", default=NULL, type='character',
              help="Input filename"),
  make_option(c("--level"), action="store", default="dataset", type='character',
              help="QC filtering peformed at the level of the whole dataset or at specific variable."),
  make_option(c("--qc_mad"), action="store", default="qc_mad.tab", type='character',
              help="QC metrics file"),
  make_option(c("--in_dir"), action="store", default='QC_statistics', type='character',
              help="Input directory: directory with summary stats for each dataset."),
  make_option(c("--out_dir"), action="store", default="QC_heatmaps", type='character',
              help="Output directory.")
)
opt = parse_args(OptionParser(option_list=option_list))

#################################### Set Variables ####################################
# Input directory
cwd <- getwd()
setwd(cwd)
main.dir <- paste0(cwd, '/')
in.dir <- paste0(opt$in_dir, '/')

# opt$datasets <- 'consortium_datasets.tab'
  
# Loading functions
functions_fn <- 'scripts/QC_functions.R'
print(paste0('Loading functions from: ',functions_fn))
source(functions_fn)

# Datasets
datasets_fn <- paste0(main.dir, opt$datasets)
datasets <- read.table(datasets_fn, header = F)$V1
datasets_label <- paste0(datasets,collapse = '.')

# Read QC-bound-MAD file
print(paste0('Reading QC metrics and MAD combinations file in: ',opt$qc_mad))
qc_mad.fn <- paste0(opt$qc_mad)
qc.mad_df <- read.table(qc_mad.fn, header = T)
qc.mad.temp_list <- split(qc.mad_df, qc.mad_df$QC_metric)
qc.mad_list <- lapply(qc.mad.temp_list, function(i) qc_mad_format(i))
qc.label <- paste(paste(qc.mad_df$QC_metric, qc.mad_df$bound, qc.mad_df$MAD_min, qc.mad_df$MAD_max, sep='_'),collapse='.')

# MAD combination order
mad_label.list <- lapply(qc.mad_list, function(x) paste0('MAD_',x$mad))
if(length(mad_label.list)==2){
  mad.order <- unlist(lapply(mad_label.list[[1]], function(i) paste0(i,':',mad_label.list[[2]])))
}else{
  mad.order_df <- expand.grid(mad_label.list)
  mad.order_df <- mad.order_df[order(mad.order_df$Var2),]
  mad.order_df$label <- paste(mad.order_df$Var1, mad.order_df$Var2, mad.order_df$Var3, sep = ':')
  mad.order <- unique(mad.order_df$label)
}

# Output directory
out.dir <- paste0(main.dir, opt$out_dir, '/', datasets_label, '/', qc.label, '/', opt$level, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

# Print report
print(paste0('Considering ', nrow(qc.mad_df), ' QC metrics: ', paste(qc.mad_df$QC_metric, collapse = ', ')))
info <- lapply(names(qc.mad_list), function(qc){
  print(qc)
  print(paste0('MAD min: ', min(qc.mad_list[[qc]][["mad"]])))
  print(paste0('MAD max: ', max(qc.mad_list[[qc]][["mad"]])))
  cat('\n')
  return(NULL)
})

#################################### Heatmaps across datasets ####################################
# Read input files in a list
print(paste0('Reading summarized data from: ', in.dir))
system.time(summarized.list <- sapply(datasets, function(i) read_in(dataset = i), simplify = F))
system.time(summarized.df <- do.call("rbind",summarized.list))

# Split dataframe by tag
system.time(tag.list <- split(summarized.df, summarized.df$tag))
  
# Pheatmap
## parameters
values <- c('n','prop','pct')
tags <- names(tag.list)
logical_vec <- c(TRUE,FALSE)
display_numbers_vec <- c('none', 'raw', 'n')
width.var <- ifelse(length(datasets)<=5, 5.5, length(datasets))

## apply function
system.time(pheatmap.res <- lapply(tags, function(t)
  lapply(values, function(v)
    lapply(display_numbers_vec, function(d)
      lapply(logical_vec, function(c) pheatmap_across_datasets(tag = t, 
                                                               value = v,
                                                               display_numbers = d,
                                                               cluster_mads = c,
                                                               tag_list = tag.list,
                                                               out_dir = out.dir,
                                                               width_var = width.var))))))