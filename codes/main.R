#==========================================#
#                                          #
#        IMC PDAC study (PRI vs MET)       #
#                                          #
#==========================================#
## Check and install missing packages
required_packages <- unique(c(
  "readxl", "harmony", "pals", "flowCore", "ggridges", "dplyr", "ggplot2", 
  "tibble", "reshape2", "cowplot", "RColorBrewer", "FNN", "stringr", "fields", 
  "tidyr", "ComplexHeatmap", "pheatmap", "circlize", "ggpubr", "tidyverse", 
  "ggsignif", "broom", "TCGAbiolinks", "SummarizedExperiment", "DESeq2", 
  "survival", "survminer","ggh4x", "BiocManager","spatstat","reshape","Hmisc","grid",
  "dittoSeq","ConsensusClusterPlus","tiff","raster",
  "FlowSOM","terra"
))


# Install missing packages
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages) > 0) {
  BiocManager::install(missing_packages)
}

# Load all required packages
lapply(required_packages, require, character.only = TRUE)

# Function to dynamically install and load required packages
load_required_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
      }, 
      BiocManager::install(pkg, dependencies = TRUE)
      )
      library(pkg, character.only = TRUE)
    }
  }
}


## load functions ====
source('codes/functions/returnfcs.R')
source('./codes/functions/do_umap.R')
source('./codes/functions/plotUmap.R')
source('./codes/functions/plot_clustering_wrapper_newscale_med.R')
source('./codes/functions/uniform_quantile_scale.R')
source('./codes/functions/plot_scaled_heatmap.R')
source('./codes/functions/visualize_ProbabilityMask.R')
source('./codes/functions/density_byCore_V2.R')
source('./codes/functions/do_CI_quantification.R')
source('./codes/functions/clusterfcs_dim.R')
source('./codes/functions/do_CN_analysis.R')
source('./codes/functions/plot_CN_cellmask2.R')
source('./codes/functions/plot_zoomed_CN.R')
source('./codes/functions/plot_multiple_probability_masks.R')
source('./codes/functions/plot_marker_level.R')


# run pre-processing if running for the first time
Specimen_designation<- read_excel('./Config/Specimen designation.xlsx')
colnames(Specimen_designation)[1]<-"sample_id"

# source('./codes/preprocessing.R')
# source('./codes/subclustering.R')

# Having backup_output.rds, start from here ====================================

# create output folder if running for the first time 
current_dir <- getwd()
output_folder <- file.path(current_dir, "output")

# Check if the folder exists
if (!dir.exists(output_folder)) {
  # Create the folder if it does not exist
  dir.create(output_folder)
  cat("output folder is created.\n")
} else {
  cat("output folder already exists.\n")
}


## load data ====
df_output <- readRDS('./data/df_output.rds')
new_expr<- readRDS('./data/new_expr.rds')

Specimen_designation$sample_id<- factor(Specimen_designation$sample_id, levels=unique(df_output$sample_id))

celltype_markers<- c("Collagen", "CD8", "CD45RA", "KI67" ,  "CD3", "CD57", 
                     "FOXP3" , "CD4" , "CD74" ,"CD86", "CD206","VISTA", 
                     "SMAVIM", "CD163" , "CK"  ,"CD15" , "CD68",  "HLADR", 
                     "Granzyme", "DCSIGN")

# customize sitelevels 
sitelevels<- c("Pancreas", "Liver")
samplevels<- unique(as.factor(df_output$sample_id))

clusterlevels=c("Immune_Mix","CD8T","CD4T","Treg", "NK", 
                "M_I","M_II","M_III","M_IV","M_V","M_VI",
                "Neutrophil","Str_I","Str_II","Str_III","Str_IV","Str_V",
                "Str_VI","Str_VII","Tumor", "UA")

# Figure 1 - analysis workflow
# Run Figure 2 & Supplementary Figure 1
source('./codes/Figure2.R')

# Run Figure 3 & Supplementary Figure 2
source('./codes/Figure3.R')

# Run Figure 4 & Supplementary Figure 3
source('./codes/Figure4.R')

# Run Figure 5 & Supplementary Figure 4
source('./codes/Figure5.R')
