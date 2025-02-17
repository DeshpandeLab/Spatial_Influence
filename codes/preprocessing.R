#===============================================================================
#         data pre-processing
#===============================================================================

### load data =====
work <- getwd()
output<- returnfcs(metaDataFile = paste0(work, '/Config/metadata.xlsx'),
                   panelDataFile = paste0(work, '/Config/panel.xlsx'),
                   dataDirectory = paste0(work, '/Data'))



index <- match(output$meta_data$sample_id, Specimen_designation$sample_id)
output$Patient <- Specimen_designation$Patient[index]
output$Site <- Specimen_designation$Site[index]


# levels
samplevels<-c(output$meta_data$sample_id)
sitelevels<-c("Pancreas","Liver")

### batch effect correction =====

# run harmony 
expr <- fsApply(output$fcs, exprs)  #create expr matrix 

md <- data.frame(sample_id = output$sample_ids)
md <- left_join(md, Specimen_designation[,c(1:3)], by="sample_id")
md <- md%>%mutate_all(as.factor)

# sanity check for IMC data - check for cores with marker expr = 0
expr_mean<- sapply(1:64, function(id) apply(expr[md$sample_id==id, ],2,mean))


### step1. logNormalize and run harmony wrt patient&site=====
expr_norm<- log2(expr +1)
set.seed(1234)
new_harmony_embeddings <- HarmonyMatrix(
  data_mat  = expr_norm[, output$cluster_by],
  meta_data = md,
  vars_use  = c("Patient", "Site")
)

### step2. shifting values to start from 0 (-minimum expr)===== 
min_marker <- apply(new_harmony_embeddings, 2, min)
new_expr <- sweep(new_harmony_embeddings, 2, min_marker, FUN = "-")

saveRDS(new_expr, './backup/new_expr.rds')

### step3. sanity check with histogram and umap===== 
pdf('./output/histogram.pdf', width = 10, height = 5)
columns <- output$cluster_by

plot_list <- list()
# Loop over each column and create the histograms
for (col_name in columns) {
  p0 <- ggplot(data.frame(Value = expr[, col_name]), aes(x = Value)) +
    geom_histogram(fill = "red", color = "black", bins = 100) +
    ggtitle(paste("expr -", col_name))
  
  p1 <- ggplot(data.frame(Value = new_expr[, col_name]), aes(x = Value)) +
    geom_histogram(fill = "blue", color = "black", bins = 100) +
    ggtitle(paste("hm_byPatient -", col_name))
  
  combined_plot <- plot_grid(p0, p1)
  print(combined_plot)
}
dev.off()


colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(md$sample_id)))
names(colorassigned)<-c(1:64)

umapRes_raw <- do_umap(expr = expr,
                       subtype_markers = output$cluster_by,
                       sample_ids = md$sample_id,
                       cell_clustering = md$sample_id, 
                       clusterMergeFile="", 
                       seed=1234, 
                       ncells=200,
                       sample_subset=NULL)

Specimen_designation$sample_id<- factor(Specimen_designation$sample_id, levels=samplevels)
umapRes_raw<- left_join(umapRes_raw, Specimen_designation[1:3], by="sample_id")

pdf('./output/Umaps_raw_Clusterby_rmbengin.pdf',width=8,height=8)
plotUmap(umapRes = umapRes_raw,
         color_clusters = colorassigned[names(colorassigned)!="NA"],
         subtype_markers =  output$cluster_by)
dev.off()


umapRes_harmony<- do_umap(expr = new_expr,
                          subtype_markers = output$cluster_by,
                          sample_ids = md$sample_id,
                          cell_clustering = md$sample_id, 
                          clusterMergeFile="", 
                          seed=1234, 
                          ncells=200,
                          sample_subset=NULL)

umapRes_harmony<- left_join(umapRes_harmony, Specimen_designation[1:3], by="sample_id")

pdf('./output/Umaps_logNormHarmony_Clusterby_rmbengin.pdf',width=8,height=8)
plotUmap(umapRes = umapRes_harmony,
         color_clusters = colorassigned[names(colorassigned)!="NA"],
         subtype_markers =  output$cluster_by)
dev.off()


### step4. cluster with 35x35 nodes and euclidean dist.===== 
## use only the canonical celltype markers for clustering
celltype_markers<- c("Collagen", "CD8", "CD45RA", "KI67" ,  "CD3", "CD57", 
                     "FOXP3" , "CD4" , "CD74" ,"CD86", "CD206","VISTA", 
                     "SMAVIM", "CD163" , "CK"  ,"CD15" , "CD68",  "HLADR", 
                     "Granzyme", "DCSIGN")

cell_clust<- clusterfcs_dim(fcs=new_expr, 
                            cluster_by = celltype_markers, 
                            xdims=35, ydims=35, # default = 35
                            numclusters=50, 
                            scaleoption = F)

head(cell_clust$cell_clustering)
table(cell_clust$cell_clustering)

### step5. plot heatmap=====
# plot only the cell type markers, and order with clustering 
plot_clustering_heatmap_wrapper_newscale_med(expr=new_expr,
                                             color_clusters = kovesi.rainbow_bgyrm_35_85_c69(50),
                                             cell_clustering = cell_clust$cell_clustering, 
                                             cluster_by=celltype_markers,
                                             fileName = './output/Clusteringheatmap_unannotated.pdf'); dev.off()


# new clustering using new_expr
output[(length(output)+1):(length(output)+3)] <- cell_clust
names(output)[(length(output)-2):(length(output))] <- c('code_clustering','cell_clustering','metaclusters')

### step6. cluster visualization=====
expr1 <- fsApply(output$fcs1, exprs)  #create expr matrix 
expr1 <-expr1[,c("X_position","Y_position")]

expr0<-data.frame(expr1,
                  cluster=cell_clust$cell_clustering,
                  sample_id=md$sample_id)

clusterdf_plot <- as.data.frame(expr0)


clus <- c(1:length(unique(expr0$cluster)))

pdf(file = "./output/clusterviz_unannotated.pdf", width = 30, height = 30)

for (i in 1:length(clus)){
  # Visualize merged clusters
  clusterdf_merged_plot <- clusterdf_plot[clusterdf_plot$cluster==clus[i],]
  
  cluster_plot <- ggplot(data=clusterdf_merged_plot, aes(x=X_position, y = 1000-Y_position)) +
    geom_point(size = 0.05) +
    facet_wrap(~ sample_id, scales = "fixed") +
    labs(
      x = "X",
      y = "Y",
      title = clus[i])
  
  cluster_plot_styled <- cluster_plot + theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 24, colour = "black"),
      axis.ticks = element_line(colour = "grey70", size = 0.2),
      panel.grid.major = element_line(colour = "grey70", size = 0.2),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "black", color = "black", linewidth = 1),
      strip.text = element_text(face = "bold", size = 10,color = "white")
    )  
  print(cluster_plot_styled)
}
dev.off()

## Manual annotation 
## annotated cluster file
clusterMergeFile = ("./Config/merge_new_annotated_f.xlsx") 
cluster_merging <- read_excel(clusterMergeFile) #annotated

## customize clusterlevels 
clusterlevels=c("CD8T",
                "CD4T",
                "Myeloid",
                "Neutrophil",
                "NK", 
                "Immune_Mix", 
                "Stroma",
                "Tumor", 
                "UA")

mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
output$cell_clustering1m <- cluster_merging$new_cluster[mm1]

saveRDS(output,'./backup/backup_output.rds')

