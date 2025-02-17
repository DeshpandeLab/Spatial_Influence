#==========================================#
#                                          #
#        IMC PDAC study (PRI vs MET)       #
#                                          #
#==========================================#

## load library ====
library(readxl)
library(harmony)
library(pals)
library(flowCore)
library(ggridges)
library(flowCore)
library(readxl)
library(dplyr)
library(ggplot2)
library(tibble)
library(reshape2)
library(cowplot)
library(RColorBrewer)
library(FNN)
library(stringr)
library(fields)
library(tidyr)
library(ComplexHeatmap)
library(pheatmap)
library(circlize)
library(ggpubr)
library(tidyverse)
library(ggsignif)
library(broom)

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

## load data ====
output <- readRDS('./backup/backup_output.rds')
new_expr<- readRDS('./backup/new_expr.rds')
Specimen_designation$sample_id<- factor(Specimen_designation$sample_id, levels=unique(output$sample_ids))

celltype_markers<- c("Collagen", "CD8", "CD45RA", "KI67" ,  "CD3", "CD57", 
                     "FOXP3" , "CD4" , "CD74" ,"CD86", "CD206","VISTA", 
                     "SMAVIM", "CD163" , "CK"  ,"CD15" , "CD68",  "HLADR", 
                     "Granzyme", "DCSIGN")

# customize sitelevels 
sitelevels<- c("Pancreas", "Liver")
samplevels<- unique(as.factor(output$sample_ids))


clusterlevels=c("Immune_Mix","CD8T","CD4T","Treg", "NK", 
                "M_I","M_II","M_III","M_IV","M_V","M_VI",
                "Neutrophil","Str_I","Str_II","Str_III","Str_IV","Str_V",
                "Str_VI","Str_VII","Tumor", "UA")


# Run Figure 1 & 1S
source('./codes/Figure1.R')

# Run Figure 2 & 2S
source('./codes/Figure2.R')

# Run Figure 3
source('./codes/Figure3.R')

# Run Figure 4
source('./codes/Figure4.R')


