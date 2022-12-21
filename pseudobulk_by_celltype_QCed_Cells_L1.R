#!/usr/bin/env Rscript

# setting working directory (cluster or local)
path_cluster <- '/gpfs/projects/bsc83/'
path_em <- '/Users/Aida/Desktop/bsc/'
path_opensuse <- '/home/aripol1/Desktop/bsc/'

if(file.exists(path_cluster)){
  setwd(paste(path_cluster))
  # .libPaths(c(.libPaths(),"/gpfs/apps/MN4/R/3.6.1-Rcpp_1.0.2/INTEL/lib64/R/library"))
  # .libPaths(c(.libPaths(),"/gpfs/apps/MN4/R/3.6.1/INTEL/lib64/R/library"))
}else if(file.exists(path_em)){
  setwd(paste(path_em))
}else if(file.exists(path_opensuse)){
  setwd(paste(path_opensuse))
}

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--dataset"), action="store", default="wg2_ng2018", type='character',
              help="Dataset."),
  make_option(c("--md_vars"), action="store", default="00.pseudobulk_by_celltype.wg2_ng2018.md_vars.tab", type='character',
              help="Donor metadata variables."),
  make_option(c("--cell_level"), action="store", default=NA, type='character',
              help="Azimuth l1 or l2."),
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="cell types in each of the cell_level"),
  make_option(c("--filter_type"), action="store", default=NA, type='character',
              help="cell-predicted.celltype.l1 // cell-predicted.celltype.l2"),
  make_option(c("--aggregate_fun"), action="store", default='sum', type='character',
              help="sum or mean --> if mean, we need to normalize before."),
  make_option(c("--assay"), action="store", default="RNA", type='character',
              help="Seurat object assay"),
  make_option(c("--mad_comb"), action="store", default="nCount_RNA_lower_1_5.percent.mt_upper_1_5", type='character',
              help="MAD combination (nCount_RNA and percent.mt)"),
  make_option(c("--so_indir"), action="store", default='wg1-preprocessing/00.so_split_by_celltype/', type='character',
              help="Input main directory"),
  make_option(c("--filt_indir"), action="store", default='wg1-qc_filtering-testing_datasets-MAD_Corr/MAD_Corr/', type='character',
              help="Barcodes kept main directory"),
  make_option(c("--out_dir"), action="store", default='wg1-preprocessing/00.pseudobulk_by_celltype/', type='character',
              help="Output main directory"))
opt = parse_args(OptionParser(option_list=option_list))

# Packages
shhh(library(data.table))
shhh(library(Seurat))
shhh(library(lme4))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(RColorBrewer))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(ggplot2))
shhh(library(ggrepel))
shhh(library(Hmisc))
shhh(library(Matrix))
#shhh(library(Matrix.utils))
shhh(library(edgeR))
shhh(library(limma))
shhh(library(textTinyR))
shhh(library(pbapply))

#################### Set Variables and load Data #################### 
# Main directory
main.dir <- 'Projects/scRNAseq/aripol1/sc-eqtlgen-consortium-pipeline/ongoing/'

# Dataset
# opt$dataset <- 'wg2_onek1k_subset'
# opt$dataset <- 'wg2_ng2018'
# opt$dataset <- 'wg2_wijst2020.V2'
# opt$dataset <- 'wg2_wijst2020.V3'

# Donor metadata variables
# opt$md_vars <- '00.pseudobulk_by_celltype.wg2_ng2018.md_vars.tab'
md_vars.fn <- paste0(main.dir, 'wg1-preprocessing/scripts/', opt$md_vars)
md_vars.df <- read.delim(md_vars.fn)
donor_id <- md_vars.df[md_vars.df$type=='donor',]$md_var
md_vars <- md_vars.df[md_vars.df$type=='other',]$md_var

# Azimuth cell level
# opt$cell_level <- 'predicted.celltype.l1'
# opt$cell_level <- 'predicted.celltype.l2'

# Azimuth cell type
# opt$cell_type <- 'CD4T' #if cell_level=predicted.celltype.l1
# opt$cell_type <- 'NK' #if cell_level=predicted.celltype.l1
# opt$cell_type <- 'T_other' #if cell_level=predicted.celltype.l1
# opt$cell_type <- 'Doublet' #if cell_level=predicted.celltype.l1
# opt$cell_type <- 'CD8_TCM' #if cell_level=predicted.celltype.l2
# opt$cell_type <- 'NK_CD56bright' #if cell_level=predicted.celltype.l2
# opt$cell_type <- 'Eryth' #if cell_level=predicted.celltype.l1 / predicted.celltype.l2

# Assay
# opt$assay <- 'RNA'
# opt$assay <- 'predicted_ADT'

# Aggregate function
# opt$aggregate_fun <- 'sum'
# opt$aggregate_fun <- 'mean'

# Output directory
# opt$filter_type <- 'cell-predicted.celltype.l1'
# opt$filter_type <- 'cell-predicted.celltype.l2'
filt_type <- str_split_fixed(opt$filter_type,'-',2)[,1]
filt_var <- str_split_fixed(opt$filter_type,'-',2)[,2]
out.tdir <- paste0(main.dir, opt$out_dir, '/', opt$dataset, '/', opt$mad_comb, '/', opt$cell_level, '/', opt$cell_type, '/', opt$assay, '/', filt_type, '/', filt_var, '/')
out.dir <- paste0(main.dir, opt$out_dir, '/', opt$dataset, '/', opt$mad_comb, '/', opt$cell_level, '/', opt$cell_type, '/', opt$assay, '/', filt_type, '/', filt_var, '/', opt$aggregate_fun, '/')
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}
print(paste0('Main results directory: ',out.dir))

# Report
print(paste0('Dataset: ', opt$dataset))
print(paste0('Cell level: ', opt$cell_level))
print(paste0('Cell type: ', opt$cell_type))
print(paste0('Filter type: ', opt$filter_type))
print(paste0('MAD combination: ', opt$mad_comb))
print(paste0('Aggegate function: ', opt$aggregate_fun))
cat('\n\n')

#################### Read data ####################
# Seurat object
so.indir <- paste0(main.dir, opt$so_indir, '/', opt$dataset, '/', opt$cell_level, '/')
pbmc_fn <- paste0(so.indir, opt$cell_type, '.rds')
print(paste0('Reading pbmc seurat object file in: ',pbmc_fn))
system.time(pbmc <- readRDS(pbmc_fn))
DefaultAssay(pbmc) <- opt$assay

# Check: compare Azimuth l1 vs. l2
table(pbmc@meta.data$predicted.celltype.l1, pbmc@meta.data$predicted.celltype.l2)
round(prop.table(table(pbmc@meta.data$predicted.celltype.l1, pbmc@meta.data$predicted.celltype.l2),margin = 1),2)

# # Downsampling (testing)
# set.seed(123)
# cells_to_subset <- sample(colnames(pbmc), 500)
# pbmc <- subset(pbmc, cells = cells_to_subset)
# pbmc@meta.data <- droplevels(pbmc@meta.data)

# Kept barcodes
filt.indir <- paste0(main.dir, opt$filt_indir, '/', opt$mad_comb, '/', opt$dataset, '/', filt_type, '/')
barcodes_fn <- paste0(filt.indir, filt_var, '.NotOutlier.mad_corr.rds')
print(paste0('Reading barcodes (kept) in: ',barcodes_fn))
system.time(barcodes_all_list <- readRDS(barcodes_fn))
names(barcodes_all_list[[1]])
barcodes_list <- lapply(barcodes_all_list, function(x){
  cell_type <- opt$cell_type
  if(!opt$cell_type%in%c('NK_CD56bright','T_other')){
    cell_type <- gsub('_', ' ', opt$cell_type)
  }
  bc_ct <- x[[cell_type]]
  return(bc_ct)
})
print('Summary of barcodes (kept):')
summary(unlist(lapply(barcodes_list, length)))
cat('\n\n')

## filtered
n_barcodes <- unlist(lapply(barcodes_list, length))
n_max <- max(n_barcodes)
bc_max <- which(n_barcodes==n_max)[1]
bc_notmax <- n_barcodes[which(n_barcodes!=n_max)]
bc_filtered <- c(names(bc_notmax), names(bc_max))
barcodes.filtered_list <- barcodes_list[names(barcodes_list)%in%bc_filtered]

## save
bc_fn <- paste0(out.tdir, 'barcodes_list.rds')
if(!file.exists(bc_fn)){
  print(paste0('Saving barcodes kept in: ', bc_fn))
  saveRDS(barcodes_list, bc_fn)
}
bc.f_fn <- paste0(out.tdir, 'barcodes.filtered_list.rds')
if(!file.exists(bc.f_fn)){
  print(paste0('Saving filtered barcodes kept in: ', bc.f_fn))
  saveRDS(barcodes.filtered_list, bc.f_fn)
}

############### Subset seurat object ################
# Function
# filt <- names(barcodes.filtered_list)[1]
# filt <- names(barcodes.filtered_list)[82]
# so = pbmc
# bc_list = barcodes.filtered_list
# out_dir = out.tdir
so_by_filt <- function(filt, so = pbmc, bc_list = barcodes.filtered_list, out_dir = out.tdir){
  barcodes <- bc_list[[filt]]
  if(!is.null(barcodes)){
    so_subset <- so[,barcodes]
    kept_cells <-  ncol(so_subset)
  }else{
    print('NO kept cells...')
    so_subset <- NULL
    kept_cells <- 0
  }
  total_cells <- ncol(so)
  prop_kept <- round(kept_cells/total_cells*100,3)
  
  # Report
  print(filt)
  print(paste0('Total nCells: ', total_cells))
  print(paste0('Kept nCells: ', kept_cells))
  print(paste0('Kept propCells: ', prop_kept, '%'))
  
  # # Save filtered seurat objects
  # if(!is.null(so_subset)){
  #   so_fn <- paste0(out_dir, filt, '.so.rds')
  #   print(paste0('Saving filtered seurat objects in: ',so_fn))
  #   system.time(saveRDS(so_subset, so_fn))
  # }
  cat('\n') 
  return(so_subset)
}

# Apply Function
system.time(so_by_filt.list <- sapply(names(barcodes.filtered_list), function(i) so_by_filt(i), simplify = FALSE))

# Filter out filters without any cells
so_by_filt.list <- Filter(Negate(is.null), so_by_filt.list)

# Save (if not exists)
so_by_filt.fn <- paste0(out.tdir, 'so_by_filt.rds')
if(!file.exists(so_by_filt.fn)){
  print(paste0('Saving seurat object list in: ', so_by_filt.fn))
  system.time(saveRDS(so_by_filt.list, so_by_filt.fn))
}

#################### Convert sc matrix to pseudo-bulk matrix (by donor) #################### 
# Function
# filt <- names(so_by_filt.list)[1]
# filt <- names(so_by_filt.list)[length(so_by_filt.list)]
# aggregate_fun = opt$aggregate_fun
# aggregates = donor_id
# so_list = so_by_filt.list
# so_assay = opt$assay
# vars_model = md_vars
# out_dir = out.dir
pseudobulk.func <- function(filt, aggregate_fun = opt$aggregate_fun, aggregates = donor_id, so_list = so_by_filt.list, so_assay = opt$assay, vars_model = md_vars, out_dir = out.dir){
  # Seurat object
  so <- so_list[[filt]]
  
  # Create the aggregate_countMatrix (gene ~ assignment)
  ## get the metadata
  metadata <- so@meta.data
  
  ## create aggregated metadata
  # aggregate_metadata <- unique(groups) #only one variable to aggregate
  # aggregate_metadata <- unique(metadata[, aggregates]) #more than one variable to aggregate
  # aggregate_metadata <- unique(metadata[, vars_model]) #we can be removing some individual (aggregates) that have the same vars_model values
  aggregate_metadata <- droplevels(unique(metadata[, c(aggregates,vars_model)]))
  
  ## set the rownames to be the same as the countsMatrix aggregate colnames, should go correctly if data is valid
  # rownames(aggregate_metadata) <- colnames(aggregate_countMatrix) # not working with large and sparse aggregate_countMatrix (omitting column names)
  rownames(aggregate_metadata) <- as.character(unique(aggregate_metadata[[aggregates]]))
  aggregate_metadata <- aggregate_metadata[, vars_model]
  
  ## get donor ids
  IDs <- metadata[,aggregates]
  unique_ID_list <- as.list(unique(IDs))
  
  # ## create the groups to aggregate on, here it's on 'assignment' --> with Matrix.utils::aggregate.Matrix()
  # groups <- data.frame(assignment=metadata[, aggregates])
  # rownames(groups) <- rownames(metadata)
  # #example for more than one variable in the metadata (e.g., sample and timepoint)
  # aggregates <- c('assignment', 'timepoint')
  # groups <- metadata[, aggregates]
  
  ## grab the countmatrix
  countMatrix <- GetAssayData(so, assay = so_assay, slot = "counts")
  # countMatrix <- so@assays[[so_assay]]@counts
  
  # Info
  n_cells <- ncol(countMatrix)
  n_genes <- nrow(countMatrix)
  print(paste0('########### Filter: ', filt, ' ###########'))
  print(paste0('Calculating pseudo-bulk expression matrix using: ', aggregate_fun))
  print(paste0('   ### n_cells: ', n_cells))
  print(paste0('   ### n_genes (sc-data): ', n_genes))
  cat('\n')
  
  if(aggregate_fun=='mean'){
    print('Normalizing the sc-data to perform mean pseudobulk...')
    cat('\n')
    
    # PF
    sampleSumInfo = colSums(countMatrix)
    meanSampleSum = mean(sampleSumInfo)
    sampleScale = sampleSumInfo/meanSampleSum
    
    # Scale to match meanRowSums
    system.time(countMatrix <- sweep(countMatrix, 2, sampleScale, FUN="/")) #divide each column by sampleScale
    
    # Log
    countMatrix = log(countMatrix+1)
    
    # PF
    sampleSumInfo = colSums(countMatrix)
    meanSampleSum = mean(sampleSumInfo)
    sampleScale = sampleSumInfo/meanSampleSum
    
    # Scale to match meanRowSums
    system.time(countMatrix <- sweep(countMatrix, 2, sampleScale, FUN="/")) #divide each column by sampleScale
    
    print('Aggregating count matrix using lapply + textTinyR::sparse_Means() ...')
    countMatrix <- as(countMatrix, "dgCMatrix")
    system.time(aggregate_countMatrix <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){sparse_Means(countMatrix[,IDs == x, drop = FALSE], rowMeans = TRUE)})))
    cellcount <- pblapply(unique_ID_list, FUN = function(x){ncol(countMatrix[,IDs == x, drop = FALSE])})
    colnames(aggregate_countMatrix) <- names(cellcount) <- unique_ID_list
    rownames(aggregate_countMatrix) <- rownames(countMatrix)
  }else{
    print('Aggregating count matrix using lapply + textTinyR::sparse_Sums() ...')
    system.time(aggregate_countMatrix <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){sparse_Sums(countMatrix[,IDs == x, drop = FALSE], rowSums = TRUE)})))
    cellcount <- pblapply(unique_ID_list, FUN = function(x){ncol(countMatrix[,IDs == x, drop = FALSE])})
    colnames(aggregate_countMatrix) <- names(cellcount) <- unique_ID_list
    rownames(aggregate_countMatrix) <- rownames(countMatrix)
  }
  
  pseudobulk_list <- NA
  geneExpr <- NA
  aggregate_countMatrix_expr.deseq2 <- NA
  
  cat('\n')
  n_sample.pseudo <- ncol(aggregate_countMatrix)
  n_genes.pseudo <- nrow(aggregate_countMatrix)
  print(paste0('   ### n_donors (pseudo-data): ', n_sample.pseudo))
  print(paste0('   ### n_genes (pseudo-data): ', n_genes.pseudo))
  cat('\n')
  
  # Save count matrices (mean = normalized, sum = raw)
  out_fn <- paste0(out_dir, filt, '.expressionCounts.rds')
  print(paste0('Saving expression counts in: ',out_fn))
  system.time(saveRDS(countMatrix, out_fn))
  cat('\n')
  
  if(aggregate_fun=='sum'){
    print('Normalizing like bulk (limma/DESeq2) using pseudobulk-sum...')
    # (If sum)
    # Filter genes by minimum expression (number of counts) -as bulk-
    isexpr = rowSums(cpm(aggregate_countMatrix)>0.1) >= 5
    isexpr_deseq2 = rowSums(aggregate_countMatrix>0) > 10
    aggregate_countMatrix_expr <- aggregate_countMatrix[isexpr,]
    aggregate_countMatrix_expr.deseq2 <- aggregate_countMatrix[isexpr_deseq2,]
    
    # Info
    n_genes_expr <- nrow(aggregate_countMatrix_expr)
    n_genes_expr.deseq2 <- nrow(aggregate_countMatrix_expr.deseq2)
    expressed_genes.prop_round <- round((n_genes_expr/n_genes),2)
    expressed_genes.deseq2.prop_round <- round((n_genes_expr.deseq2/n_genes),2)
    print(paste0('   ### n_cells: ', n_cells))
    print(paste0('   ### Limma voom --> n_genes (expressed): ', n_genes_expr, ' (proportion of ', expressed_genes.prop_round,')'))
    print(paste0('   ### DESeq2 --> n_genes (expressed): ', n_genes_expr.deseq2, ' (proportion of ', expressed_genes.deseq2.prop_round,')'))
    if(n_genes_expr>2){
      # Normalization (standard usage of limma/voom)
      geneExpr = DGEList(aggregate_countMatrix_expr)
      geneExpr = calcNormFactors(geneExpr)
      
      # Gene expression matrices (for the correlation analysis) --> use options 2,3 and 4; but save all of them
      ## 1. raw pseudo-bulk counts (sum counts from cells from the same cell type)
      counts <- geneExpr$counts
      
      ## 2. log2 pseudo-bulk counts (sum counts from cells from the same cell type)
      counts_log <- log2(counts+0.5)
      # counts_log <- log2(counts+1)
      
      ## 3. log2CPM normalized using normalization factors -voom- (previously calculated with calcNormFactors() from edgeR)
      # from source code of voom() --> https://rdrr.io/bioc/limma/src/R/voom.R --> log2-counts-per-million
      # y <- t(log2(t(counts+0.5)/(lib.size+1)*1e6))
      # y <- normalizeBetweenArrays(y,method=normalize.method) #from https://rdrr.io/bioc/limma/src/R/norm.R --> in normalizeBetweenArrays() it seems if object=matrix, quantile normalization
      # out$E <- y #out represents 'voom.out' object in this script
      voom.out <- voom(geneExpr) #Transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observation-level weights. The data are then ready for linear modelling.
      counts_log_norm <- voom.out$E
      
      ## 4. log2 of weighted counts by number of cells per donor
      md.ss <- metadata[c(aggregates,'Barcode')]
      md.ss %>%
        group_by_at(aggregates) %>%
        count() %>%
        as.data.frame() -> cells_by_donor
      w <- cells_by_donor$n # weight
      counts_w.t <- t(counts)/w
      counts_weighted <- t(counts_w.t)
      counts_log_weighted <- log2(counts_weighted+0.5)
      # counts_log_weighted <- log2(counts_weighted+1)
      
      pseudobulk_list <- list(counts = counts,
                              counts_log = counts_log,
                              counts_log_norm = counts_log_norm,
                              counts_weighted = counts_weighted,
                              counts_log_weighted = counts_log_weighted)
    }else{
      aggregate_countMatrix_expr.deseq2 <- NA
      print(paste0('Pseudobulk have not been performed for filter: ', filt, '. voom() needs at least 2 genes to fit a mean-variance trend. We have ', n_genes_expr, ' genes expressed.'))
      # err_in <- paste0('Pseudobulk have not been performed for filter: ', filt, '. voom() needs at least 2 genes to fit a mean-variance trend. We have ', n_genes_expr, ' genes expressed.')
      # stop(err_in)
    }
    cat('\n')
  }
 
  ## output
  out <- list(pseudobulk_LimmaVoom = pseudobulk_list,
              geneExpr_LimmaVoom = geneExpr,
              geneExpr_DESeq2 = aggregate_countMatrix_expr.deseq2, 
              aggregate_countMatrix = aggregate_countMatrix,
              metadata = aggregate_metadata)
  
  # Save pseudobulk matrices
  out_fn <- paste0(out_dir, filt, '.pseudobulk_matrices.rds')
  print(paste0('Saving pseudobulk matrices in: ',out_fn))
  system.time(saveRDS(out, out_fn))
  cat('\n')
  
  return(out)
}

# Apply Function
system.time(pseudobulk.list <- sapply(names(so_by_filt.list), function(i) pseudobulk.func(i), simplify = FALSE))
# summary(rowMeans(pseudobulk.list[[1]]$aggregate_countMatrix))

# Save the whole object
pseudobulk.list_fn <- paste0(out.dir, 'pseudobulk_matrices.rds')
print(paste0('Saving all pseudobulk matrices in: ',pseudobulk.list_fn))
system.time(saveRDS(pseudobulk.list, pseudobulk.list_fn))