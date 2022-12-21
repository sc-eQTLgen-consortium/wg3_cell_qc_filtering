#!/usr/bin/env Rscript

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--dataset"), action="store", default=NULL, type='character',
              help="Name of the dataset."),
  make_option(c("--level"), action="store", default="dataset", type='character',
              help="Level at which the QC filtering is peformed. By default, it will be performed at the whole dataset. You can chose to peform it at the level of a specific metadata variable (such as, lane/pool)."),
  make_option(c("--md_vars"), action="store", default=NULL, type='character',
              help="Metadata variables file (optional). By default, the QC statistics will be calculated at the whole dataset. You can choose to calcualte them by metadata variable."),
  make_option(c("--qc_mad"), action="store", default="qc_mad.tab", type='character',
              help="QC metrics file. By default, the QC metrics-thresholds will be: lower nCount_RNA and upper percent.mt, from MAD 1 to 5."),
  make_option(c("--azimuth_pairing"), action="store", default="azimuth_l1_l2.csv", type='character',
              help="Azimuth l1-l2 pairing file."),
  make_option(c("--downsampling"), action="store", default=NULL, type='character',
              help="Downsample file with metadata variable and number of cells (optional)."),
  make_option(c("--in_dir"), action="store", default=NULL, type='character',
              help="Input directory: WG2 output directory."),
  make_option(c("--out_dir"), action="store", default="QC_statistics", type='character',
              help="Output directory.")
)
opt = parse_args(OptionParser(option_list=option_list))

#################################### Set Variables ####################################
 # Input directory
cwd <- getwd()
setwd(cwd)
in.dir <- paste0(opt$in_dir, '/')

# Loading functions
functions_fn <- '/tools/wg1-qc_filtering/scripts/QC_functions.R'
print(paste0('Loading functions from: ',functions_fn))
source(functions_fn)

# Dataset
dataset <- opt$dataset
print(paste0('Running QC_statistics.R on: ', dataset))

# Read azimuth_l1_l2 pairing file
print(paste0('Reading azimuth_l1_l2 pairing file in: ',opt$azimuth_pairing))
pairs.fn <- opt$azimuth_pairing
pairs.df <- read.csv2(pairs.fn)
colnames(pairs.df) <- c('l1','l2')

# Read QC-bound-MAD file
print(paste0('Reading QC metrics and MAD combinations file in: ',opt$qc_mad))
qc_mad.fn <- paste0(opt$qc_mad)
qc.mad_df <- read.table(qc_mad.fn, header = T)
qc.mad.temp_list <- split(qc.mad_df, qc.mad_df$QC_metric)
qc.mad_list <- lapply(qc.mad.temp_list, function(i) qc_mad_format(i))

# Setting MAD combination order
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
out.label <- paste(paste(qc.mad_df$QC_metric, qc.mad_df$bound, qc.mad_df$MAD_min, qc.mad_df$MAD_max, sep='_'),collapse='.')
out.dir <- paste0(opt$out_dir, '/', dataset, '/', out.label,'/')
if(!is.null(opt$md_vars)){
  out.dir <- paste0(out.dir, 'by_metadata/')
}else{
  out.dir <- paste0(out.dir, 'by_dataset/')
}
out.dir <- paste0(out.dir, opt$level, '/')
if(!is.null(opt$downsampling)){
  downsampling.fn <- opt$downsampling
  downsampling_df <- read.table(downsampling.fn, header = T)
  n <- downsampling_df$n
  md_ds <- downsampling_df$md_var
  print(paste0('Downsampling ', n, ' cells from ', md_ds))
  label_ds <- paste(md_ds, n, sep='_')
  out.dir <- paste0(out.dir, label_ds, '/')
}
print(paste0('Creating output directory in: ',out.dir))
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

#################################### Summarize ####################################
# Read input data (metadata or seurat object in case metadata is not there)
md.fn <- paste0(in.dir, dataset, '/step4_reduce/metadata.reduced_data.RDS')
so.fn <- paste0(in.dir, dataset, '/step4_reduce/reduced_data.RDS')
if(file.exists(md.fn)){
  print(paste0('Reading metadata in: ', md.fn))
  system.time(so_metadata <- readRDS(md.fn))
}else if(file.exists(so.fn)){
  print(paste0('Reading seurat object in: ', so.fn))
  system.time(so <- readRDS(so.fn))
  print('Extracting metadata from the seurat object...')
  system.time(so_metadata <- so@meta.data)
}else{
  err_in <- paste0('Neither metadata object (metadata.reduced_data.RDS) nor seurat object (reduced_data.RDS) are found in the input directory: ',in.dir)
  stop(err_in) 
}

# Add Azimuth l1 classification to the metadata
print('Adding Azimuth l1 classification to the metadata...')
system.time(so_metadata <- add_azimuth_l1(so_md = so_metadata, pairing = pairs.df))

# Downsampling (optional)
if(!is.null(opt$downsampling)){
  # split metadata by metdata variable to downsample
  system.time(so_md.list <- split(so_metadata, so_metadata[[md_ds]]))
  
  # apply function to each element of the metadata list
  system.time(so_md_ds.list <- lapply(so_md.list, function(i) downsample_by_metadata(so_md = i,
                                                                                     n_cells = n,
                                                                                     md_var = md_ds)))
  # rbind
  system.time(so_metadata <- do.call("rbind", c(so_md_ds.list, make.row.names = FALSE)))
  rownames(so_metadata) <- so_metadata$Barcode
}

# Calculate MADs and assigning Outlier/NotOutlier tags per cell barcode
print('Calculating MADs and assigning Outlier/NotOutlier tags per cell barcode...')
system.time(qc_tag.df <- qc_tag(so_md = so_metadata, 
                                qc_mad_list = qc.mad_list,
                                filter_level = opt$level))
# Save extended output
qc_tag_fn <- paste0(out.dir,'qc_tag.rds')
print(paste0('Saving extended output in: ', qc_tag_fn))
system.time(saveRDS(qc_tag.df, qc_tag_fn))

# Summarize by dataset or by metadata
print('Summarizing the data...')
if(!is.null(opt$md_vars)){
  # Defining metadata variables
  md_vars.fn <- opt$md_vars
  print(paste0('By metadata variables in: ',md_vars.fn))
  md_df <- read.table(md_vars.fn, header = T)
  md_vec <- unique(md_df$md_var)
  md_list <- sapply(unique(md_df$type), function(i) unique(md_df[md_df$type==i,]$md_var), simplify = F)
  md_n.vec.order <- sapply(md_vec, function(i) table(so_metadata[[i]]), simplify = F)
  md_vec.order <- sapply(md_vec, function(i) names(sort(table(so_metadata[[i]]),decreasing = T)), simplify = F)
  md_vec.order.fn <- paste0(out.dir, 'md_order.rds')
  print(paste0('Saving metadata levels ordered by nCells in: ', md_vec.order.fn))
  system.time(saveRDS(md_vec.order, md_vec.order.fn))
  
  # Summarize
  system.time(count_by_dataset.res <- sapply(names(md_list), function(i) count_by_dataset(qc_tag = qc_tag.df,
                                                                                          md_type = i), simplify = F))
  system.time(tag.out <- lapply(count_by_dataset.res, function(md_type)
    lapply(md_type, function(md_var) split(md_var, md_var$tag))))
  
  # Pheatmap
  print('Plotting heatmaps per each dataset and metadata type...')

  ## parameters
  width.list <- lapply(md_vec.order, function(x){
    x.len <- length(x)
    if(x.len<5){
      width <- 5
    }else{
      width <- round(x.len*0.75,2)
    }
    return(width)
  })
  md_types <- names(md_list)

  ## main pheatmap function
  system.time(pheatmap_main.res <- sapply(md_types, function(i) pheatmap_main(md_type = i,
                                                                              tag_list = tag.out, 
                                                                              out_dir = out.dir,
                                                                              width_list = width.list), simplify = F))
}else{
  print('By dataset...')
  
  # Summarize
  system.time(tag.out <- count_by_dataset(qc_tag = qc_tag.df))
}

# Save summarized output
tag_fn <- paste0(out.dir,'tag.rds')
print(paste0('Saving summarized output in: ', tag_fn))
system.time(saveRDS(tag.out, tag_fn))
