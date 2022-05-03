# sc-eQTLGen Quality Control (QC) Threshold Selection Committee
We provide an [add-on script](QC_statistics.R) for the QC to generate the following statistics/plots when setting different MAD (Median Absolute Deviation) combinations:

* **summary statistics** of lost/kept cells for the whole dataset and/or by metadata variable. The metadata variables that should be mandatory for every group are: cell type classification (predicted.celltype.l1/predicted.celltype.l2/scpred_prediction, predicted from WG2) and the batch (sequencing lane/pool).

* **heatmap plots** for the metadata variables showing the percentage and number of lost/kept cells.

*Of note*: To run this script you should have successfully run WG1 and WG2 sc-eQTLGen consortium pipelines.

## Contact
If you have any questions or issues, feel free to open an issue or directly email Aida Ripoll-Cladellas (aida.ripoll@bsc.es)


## Required Software
* **R** >=4.0.0 version: You need to install the packages loaded in the [add-on script](QC_statistics.R) and in the [external script](scripts/QC_functions.R).

## Required Input
This section explains the input data and it’s structure to run the [add-on script](QC_statistics.R).

### Test Data
We have provided a **test dataset** *(wg2_onek1k_subset)* that contains one pool of a 10x run from the [**OneK1K** dataset](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02293-3). Notice that it is a significantly down-sized and sub-sampled version of the whole dataset. In this test dataset, the total number of cells is 1,207 from 13 donors.

Here is the structure of the [input directory for the test dataset](/wg2-cell_type_classification/wg2_onek1k_subset/). This input directory (*/wg2-cell_type_classification/wg2_onek1k_subset/*) should have the same structure as the WG2 pipeline output directory. We will need only the files in the [step4_reduce](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/) directory:

**wg2-cell_type_classification**    
└── wg2_onek1k_subset  
    ├── cell_classification.sif  
    ├── map_hierscpred.R  
    ├── schier_workaroung.sh  
    ├── step1_split  
    │   └── OneK1K-test_dataset.RDS  
    ├── step2_azimuth  
    │   ├── OneK1K-test_dataset.RDS  
    │   ├── OneK1K-test_dataset_ref_spca.png  
    │   └── OneK1K-test_dataset_ref_umap.png  
    ├── step3_hierscpred  
    │   └── OneK1K-test_dataset.RDS  
    ├── **step4_reduce**   
    │   ├── **metadata.reduced_data.RDS**    
    │   └── **reduced_data.RDS**    
    └── step5_compare  
        ├── comparison_contingency_table.tsv  
        ├── comparison_heatmap_counts.pdf  
        └── comparison_heatmap_prop.pdf  
        
The main input for the [add-on script](QC_statistics.R) is the metadata slot ([metadata.reduced_data.RDS](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/metadata.reduced_data.RDS)) of the seurat object provided by WG2 pipeline ([reduced_data.RDS](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/reduced_data.RDS)). The WG2 pipeline is peforming the cell type classification of the non-QC filtered singlets predicted by WG1 pipeline.

* **Recommended:** We recommend you to use the WG2 seurat object's metadata slot ([metadata.reduced_data.RDS](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/metadata.reduced_data.RDS)). It will speed up the running time and memory of the script. 

* **Alternative:** You can also the WG2 seurat object ([reduced_data.RDS](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/reduced_data.RDS)). However, it will slow down the running time and memory of the script as we will need to read the full seurat object which can be very large depending on the number of cells (e.g., ~77K cells, 8.9G). 

### Required Data
**wg1-qc_filtering**  
|-- azimuth_l1_l2.csv    
|-- downsampling.tab    
|-- metadata_variables.tab    
|-- qc_mad.tab    
|-- wg2-cell_type_classification    

#### QC-MAD combinations ([qc_mad.tab](qc_mad.tab))
A tsv file that has in the:
* 1st column: QC metrics. By default, number of UMIs (*nCount_RNA*) and % of mitochondrial genes (*percent.mt*).
* 2nd column: Upper or lower threshold. By default, lower for *nCount_RNA* and upper for *percent.mt*.
* 3rd and 4rd columns: minimum and maximum MADs. By default, *minimum*=1 and *maximum*=5.

*Of note*:
* Tab separated
* It is assumed that the QC metrics are calculated in the seurat object as a result from WG1 pipeline, and thus, they are columns of the metadata slot of the seurat object.
* This file must have this header. 
* The QC-MAD combinations file provided for the test dataset is the [qc_mad.tab](/qc_mad.tab) file:

| QC_metric  | bound | MAD_min  | MAD_max |
| ------------- | ------------- | ------------- | ------------- |
| nCount_RNA  | lower  | 1  | 5 |
| percent.mt  | upper  | 1  | 5 |

#### Azimuth l1-l2 pairing file ([azimuth_l1_l2.csv](/azimuth_l1_l2.csv))
A csv file that has in the:
* 1st column: Azimuth's level 1 cell type classification (L1).
* 2nd column: Azimuth's level 2 cell type classification (L2).

*Of note*:
* Semmicolon separated.
* It is assumed that the Azimuth's level 2 classification is predicted from WG2 pipeline from WG1 pipeline, whereas Azimuth's level 1 has been manually  defined to make a broader cell type classification.  
* This file must have this header. 
* The Azimuth l1-l2 pairing file provided for the test dataset is the [azimuth_l1_l2.csv](/azimuth_l1_l2.csv) file:

| L1 | L2 |  
| ------------- | ------------- |  
| CD4T  | Treg  |  
| CD4T | CD4 Naive |  
| CD4T| CD4 TCM |  
| CD4T| CD4 TEM |   
| CD4T | CD4 CTL |   
| CD4T | CD4 Proliferating | 
| CD8T | CCD8 Naive  |  
| CD8T | CD8 TCM  |  
| CD8T | CD8 TEM  |  
| CD8T | CD8 Proliferating  |  
| T_other | MAIT  |  
| T_other | dnT  |  
| T_other | gdT  |  
| T_other | ILC  |  
| NK | NK  |  
| NK | NK Proliferating  |  
| NK | NK_CD56bright  |  
| Mono | CD14 Mono  |  
| Mono | CD16 Mono  |  
| DC | cDC1  |  
| DC | cDC2  |  
| DC | pDC  |  
| DC | ASDC  |  
| B | B naive  |  
| B | B intermediate  |  
| B | B memory  |  
| B | Plasmablast  |  
| HSPC | HSPC  |  
| Platelet | Platelet  |  
| Eryth | Eryth  |  


### Optional Data
#### Metadata variables ([metadata_variables.tab](/metadata_variables.tab))
A tsv file that has in the:
* 1st column: Metadata variable name. 
* 2nd column: Metadata variable type. 
* 3rd and 4rd columns: minimum and maximum MADs. By default, *minimum*=1 and *maximum*=5.

*Of note*:
* Tab separated.
* It is assumed that the metadata variable names are columns of the metadata file or metadata slot of the seurat object.
* This file must have this header.
* By default, the QC statistics will be summarized at the whole dataset. You can choose to summarize them by metadata variable.
* The metadata variables file provided for the test dataset is the [metadata_variables.tab](/metadata_variables.tab) file: 

| md_var  | type |  
| ------------- | ------------- |  
| Pool  | donor  |  
| Assignment  | donor  |  
| predicted.celltype.l2  | cell  |  
| scpred_prediction  | cell  |  
| predicted.celltype.l1  | cell  |  


#### Downsampling file ([downsampling.tab](/downsampling.tab))
A tsv file that has in the:
* 1st column: Metadata variable name. 
* 2nd column: Number of cells to use for downsampling every level of the specified metadata variable.

*Of note*:
* Tab separated.
* It is assumed that the metadata variable name is a column of the metadata file or metadata slot of the seurat object.
* This file must have this header.
* By default, the QC statistics will be calculated using the whole dataset. You can choose to downsample the whole dataset to a specific number of cells *(n)* for each level of a specific metadata variable *(md_var)*.
* The downsampling file provided for the test dataset is the [downsampling.tab](/downsampling.tab) file: 

| md_var  | n |  
| ------------- | ------------- |  
| predicted.celltype.l1  | 100  |  

## Running the [add-on script](QC_statistics.R)
*Of note*: The functions called in the [add-on script](QC_statistics.R) are defined in an [external script](scripts/QC_functions.R).

**1.** The first step is to clone this repository and change your current working directory.    
```
git clone https://github.com/aidarripoll/wg1-qc_filtering.git  
cd wg1-qc_filtering
```

**2.** Set common environmental variables:  
```
dataset=wg2_onek1k_subset  
input_directory=wg2-cell_type_classification
output_directory=QC_statistics_examples
```

**3.** Running the add-on script with different parameters:  

3.1. Summarize the QC statistics at the dataset level after:

* Calculating the QC statistics at the dataset level:
```
Rscript QC_statistics.R --dataset $dataset --in_dir $input_directory --out_dir $output_directory
```

* Calculating the QC statistics at the batch metadata variable (i.e., Pool) level:
```
batch_variable=Pool  
Rscript QC_statistics.R --dataset $dataset --level $batch_variable --in_dir $input_directory --out_dir $output_directory 
```

The output for each of the parameters settings is the summary statistics (*tag.rds*):
QC_statistics_examples  
└── wg2_onek1k_subset  
    └── nCount_RNA_lower_1_5.percent.mt_upper_1_5  
        └── **by_dataset**   
            ├── dataset  
            │   └── **tag.rds**  
            └── Pool  
                └── **tag.rds**  

3.2. Summarize the QC statistics at the metadata level after:
```
metadata_vars=metadata_variables.tab
```

* Calculating the QC statistics at the dataset level:
```
Rscript QC_statistics.R --dataset $dataset --md_vars $metadata_vars --in_dir $input_directory --out_dir $output_directory
```

* Calculating the QC statistics at the batch metadata variable (i.e., Pool) level:
```
batch_variable=Pool  
Rscript QC_statistics.R --dataset $dataset --level $batch_variable --md_vars $metadata_vars --in_dir $input_directory --out_dir $output_directory 
```

The outputs for each of the parameters settings are the summary statistics (*tag.rds*) and the heatmap plots for each of the metadata variable, which are organized by the metadata variable type (e.g., cell or donor):
QC_statistics_examples  
└── wg2_onek1k_subset  
    └── nCount_RNA_lower_1_5.percent.mt_upper_1_5  
        ├── by_dataset  
        │   └── Pool  
        │       └── tag.rds  
        └── **by_metadata**  
            ├── dataset  
            │   ├── **cell**  
            │   ├── **donor**  
            │   ├── md_order.rds  
            │   └── **tag.rds**  
            └── Pool  
                ├── **cell**  
                ├── **donor**  
                ├── md_order.rds  
                └── **tag.rds**  


3.3. Downsampling (optional): By default, the QC statistics will be calculated using the whole dataset. You can choose to downsample the whole dataset to a specific number of cells *(n)* for each level of a specific metadata variable *(md_var)* by adding the `--downsampling` parameter in the previous commands in the 3.1 and 3.2 sections. 
```
downsampling_file=downsampling.tab
Rscript QC_statistics.R --dataset $dataset --downsampling $downsampling_file --in_dir $input_directory --out_dir $output_directory  
Rscript QC_statistics.R --dataset $dataset --level $batch_variable --downsampling $downsampling_file --in_dir $input_directory --out_dir $output_directory  
Rscript QC_statistics.R --dataset $dataset --md_vars $metadata_vars --downsampling $downsampling_file --in_dir $input_directory --out_dir $output_directory  
Rscript QC_statistics.R --dataset $dataset --level $batch_variable --md_vars $metadata_vars --downsampling $downsampling_file --in_dir $input_directory --out_dir $output_directory
```

The outputs will be the same as in 3.1 and 3.2, for example the 3.1 outputs will be in *predicted.celltype.l1_100* directory.
by_dataset  
├── dataset  
│   └── **predicted.celltype.l1_100**  
│       └── tag.rds  
└── Pool
    ├── **predicted.celltype.l1_100**  
    │   └── tag.rds  
    └── tag.rds  


## Discussion in the QC Threshold Selection Committee
To decide the final QC threshold selection criteria, we will need that you follow the next steps:  
1. Run all the commands in the **sections 3.1 and 3.2**.
2. Send the outputs to us by:  

* **Recommended:** We recommend that you tar and gunzip all the outputs.
```
tar -cvf ${output_directory}.tar.gz $output_directory
```

* **Alternative:** You can also extract the main outputs by quickly running an [extra script](QC_extract_files.R). These outputs will be the ones we will focus on for further discussions within the QC threshold selection committee.
```
Rscript QC_extract_files.R --dataset $dataset --level $batch_variable --in_dir $output_directory
```

Then, you can tar and gunzip these files.

```
output_directory_files=${output_directory}.files
tar -cvf ${output_directory_files}.tar.gz $output_directory_files
```

## Example outputs

We have provided the two output directories for the test data *(wg2_onek1k_subset)* in a tar.gz format:

* [QC_statistics_examples.tar.gz](QC_statistics_examples.tar.gz): Outputs from running the commands in **3.1, 3.2 and 3.3** as part of **Running the add-on script** section.

* [QC_statistics_examples.files.tar.gz](QC_statistics_examples.files.tar.gz): Outputs from running the commands in **3.1 and 3.2** as part of **Running the add-on script** section, and also running the **alternative** option in the **Discussion in the QC Threshold Selection Committee** section.


You can decompress them by:
```
tar -xvf QC_statistics_examples.tar.gz
tar -xvf QC_statistics_examples.files.tar.gz
```

## Running time and memory requirements
* **[add-on script]**(QC_statistics.R):
To speed up the running time and memory requirements of the **[add-on script]**(QC_statistics.R), we recommend to submit each of the commands in **3.1 and 3.2** of the **Running the add-on script** section as an independent job on your HPC infrastructure (i.e., run each job as an element of a job array). The running time and memory requirements will depend on:  
1. The **size** of your dataset. Notice that the test dataset is a significantly down-sized and sub-sampled version of the whole dataset (# of cells=1,207 and # of donors=13). 
2. Whether you already have the **metadata slot** ([metadata.reduced_data.RDS](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/metadata.reduced_data.RDS)) of the seurat object provided by WG2 pipeline, or you only have the **whole seurat object** [reduced_data.RDS](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/reduced_data.RDS) provided by WG2 pipeline. If possible, you should use the WG2 seurat object's metadata slot ([metadata.reduced_data.RDS](/wg2-cell_type_classification/wg2_onek1k_subset/step4_reduce/metadata.reduced_data.RDS)).

*Of note*: We ran this [add-on script](QC_statistics.R) on a larger dataset in [Oelen et al, 2020, V2 dataset](https://www.biorxiv.org/content/10.1101/2021.06.04.447088v1) (# of cells=480,503 and # of donors=88) taking as the main input the metadata slot of the seurat object provided by WG2 pipeline using the following SLURM parameters: `--cpus-per-task=48` and `--nodes=1`. 

