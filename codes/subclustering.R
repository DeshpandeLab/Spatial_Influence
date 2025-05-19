#--------------------------------#
#          subclustering         #
#--------------------------------#

new_expr <- readRDS('./backup/new_expr.rds') # batch effect corrected expr
output <- readRDS('./backup/backup_output.rds')
output$cell_clustering2m<- output$cell_clustering1m

colorassigned <- hcl.colors(10, palette = "dynamic")

#### 1. CD4T =====
sub_celltype <- "CD4T"
subcluster_markers<- c("CD45RA", "CD45RO" , "CD3", 
                     "FOXP3" , "CD4")
mm1<- which(output$cell_clustering1m==sub_celltype)

# extract x,y coord for visualization
coord<- fsApply(output$fcs1, exprs)[mm1, c("X_position", "Y_position")]
sub_df<- new_expr[mm1, ]

cell_clust<- clusterfcs_dim(fcs=sub_df, 
                            cluster_by = subcluster_markers, 
                            xdims=10, ydims=10, # default = 35, for smaller subset -> use 10 (35 gave error)
                            numclusters=10, 
                            scaleoption = F)
unique(cell_clust$cell_clustering)
table(cell_clust$cell_clustering)
sum(table(cell_clust$cell_clustering))

sub_df<- as.data.frame(sub_df)
sub_df2<- data.frame(sub_df, 
                     coord, 
                     sample_id=output$sample_id[mm1])
sub_df2$cell_clustering<- cell_clust$cell_clustering

plot_scaled_heatmap(expr = sub_df2,
                    markers = subcluster_markers, 
                    colorassigned = colorassigned, 
                    cell_clustering1m = as.factor(sub_df2$cell_clustering),
                    heatmapcluster = T, 
                    filename =paste0("./output/", sub_celltype, "_subclustering_heatmap_unannotated.pdf"))

# manual annotation
sub_df2 <- sub_df2 %>%
  mutate(cell_clustering1m = case_when(
    cell_clustering %in% c(1,7,4,2) ~ "CD4T",
    cell_clustering %in% c(10,8,9,6,5,3) ~ "Treg"
  ))

# assign subtype annotations for CD4T
output$cell_clustering2m[mm1]<- sub_df2$cell_clustering1m

plot_scaled_heatmap(expr = sub_df2,
                    markers = subcluster_markers, 
                    colorassigned = colorassigned, 
                    cell_clustering1m = as.factor(sub_df2$cell_clustering1m), 
                    heatmapcluster = T, 
                    filename =paste0("./output/", sub_celltype, "_subclustering_heatmap_annotated.pdf"))


## sub-cluster visualization for all sample ids
pdf(file = paste0("./output/",sub_celltype,"_subclustering_plot.pdf"), width = 30, height = 30)
cluster_plot <- ggplot(data=sub_df2, 
                       aes(x=X_position, 
                           y = 1000-Y_position, 
                           color=as.factor(cell_clustering1m))) +
  geom_point(size = 0.05) +
  facet_wrap(~ sample_id, scales = "fixed") +
  labs(
    x = "X",
    y = "Y")

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
dev.off()


#### 2. Other Myeloids =====
sub_celltype <- "Myeloid" 
subcluster_markers<- c("CD163", "CD206", "HLADR",
                     "CD68", "CD86", "DCSIGN", "ARG1") 

mm3<- which(output$cell_clustering1m==sub_celltype)

# extract x,y coord for visualization
coord<- fsApply(output$fcs1, exprs)[mm3, c("X_position", "Y_position")]
sub_df<- new_expr[mm3, ]


source('./codes/functions/clusterfcs_dim.R')
cell_clust<- clusterfcs_dim(fcs=sub_df, 
                            cluster_by = subcluster_markers, 
                            xdims=10, ydims=10, 
                            numclusters=10, 
                            scaleoption = F)

unique(cell_clust$cell_clustering)
table(cell_clust$cell_clustering)
sum(table(cell_clust$cell_clustering))

sub_df<- as.data.frame(sub_df)
sub_df2<- data.frame(sub_df, 
                     coord, 
                     sample_id=output$sample_id[mm3])
sub_df2$cell_clustering<- cell_clust$cell_clustering
table(sub_df2$cell_clustering, sub_df2$sample_id)

# manual annotation
sub_df2 <- sub_df2 %>%
  mutate(cell_clustering1m = case_when(
    cell_clustering %in% c(10,8) ~ "M_II",
    cell_clustering %in% c(9) ~ "M_VI",
    cell_clustering %in% c(6) ~ "M_V",
    cell_clustering %in% c(5, 7) ~ "M_III",
    cell_clustering %in% c(4) ~ "M_IV", 
    cell_clustering %in% c(1,2,3) ~ "M_I"
  ))


plot_scaled_heatmap(expr = sub_df2,
                    markers = subcluster_markers, 
                    colorassigned = colorassigned, 
                    cell_clustering1m = as.factor(sub_df2$cell_clustering),
                    heatmapcluster = T, 
                    filename =paste0("./output/", sub_celltype, "_subclustering_heatmap_unannotated.pdf"))

# assign subtype annotations for macrophages
output$cell_clustering2m[mm3]<- sub_df2$cell_clustering1m


# include all other cells to compare
plot_scaled_heatmap(expr = new_expr,
                    markers = subcluster_markers, 
                    colorassigned = colorassigned, 
                    cell_clustering1m = as.factor(output$cell_clustering2m), 
                    heatmapcluster = T, 
                    filename =paste0("./output/", sub_celltype, "_subclustering_heatmap_all.pdf"))


sub_df2$sample_id<-as.factor(sub_df2$sample_id)
sub_df2$cell_clustering1m<- as.factor(sub_df2$cell_clustering1m)
df<- melt(sub_df2[,c(subcluster_markers, "sample_id", "cell_clustering1m")])

pdf(paste0("./output/", sub_celltype, "_markerExpr_ridgeplot.pdf"),
    height=10, width=8)
ggplot(df[df$variable%in%subcluster_markers, ], aes(x = value, y = cell_clustering1m, fill = cell_clustering1m)) +
  geom_density_ridges()+
  facet_wrap(~variable, scales = "free", ncol = 2)
dev.off()

# subset only the macrophages
plot_scaled_heatmap(expr = sub_df2,
                    markers = subcluster_markers, 
                    colorassigned = colorassigned, 
                    cell_clustering1m = as.factor(sub_df2$cell_clustering1m), 
                    heatmapcluster = T, 
                    filename =paste0("./output/", sub_celltype, "_subclustering_heatmap_annotated.pdf"))




#### 3. Stroma =====
sub_celltype <- "Stroma" 
subcluster_markers<- c("Collagen", "Podoplanin", "SMAVIM", "HLADR", "CD74") 

mm4<- which(output$cell_clustering1m==sub_celltype)

# extract x,y coord for visualization
coord<- fsApply(output$fcs1, exprs)[mm4, c("X_position", "Y_position")]
sub_df<- new_expr[mm4, ]

source('./codes/functions/clusterfcs_dim.R')
cell_clust<- clusterfcs_dim(fcs=sub_df, 
                            cluster_by = subcluster_markers, 
                            xdims=10, ydims=10, 
                            numclusters=10, 
                            scaleoption = F)
unique(cell_clust$cell_clustering)
table(cell_clust$cell_clustering)
sum(table(cell_clust$cell_clustering))

sub_df<- as.data.frame(sub_df)
sub_df2<- data.frame(sub_df, 
                     coord, 
                     sample_id=output$sample_id[mm4])
sub_df2$cell_clustering<- cell_clust$cell_clustering

# manual annotation
sub_df2 <- sub_df2 %>%
  mutate(cell_clustering1m = case_when(
    cell_clustering %in% c(9) ~ "Str_V",
    cell_clustering %in% c(1,5) ~ "Str_II",
    cell_clustering %in% c(3) ~ "Str_III",
    cell_clustering %in% c(8) ~ "Str_VI",
    cell_clustering %in% c(10) ~ "Str_VII",
    cell_clustering %in% c(2,4) ~ "Str_I",
    cell_clustering %in% c(6,7) ~ "Str_IV"
  ))

plot_scaled_heatmap(expr = sub_df2,
                    markers = subcluster_markers, 
                    colorassigned = colorassigned, 
                    cell_clustering1m = as.factor(sub_df2$cell_clustering),
                    heatmapcluster = T, 
                    filename =paste0("./output/", sub_celltype, "_subclustering_heatmap_unannotated.pdf"))


# assign annotation for stroma subtypes 
output$cell_clustering2m[mm4]<- sub_df2$cell_clustering1m

# include all other cells to compare
plot_scaled_heatmap(expr = new_expr,
                    markers = subcluster_markers, 
                    colorassigned = colorassigned, 
                    cell_clustering1m = factor(output$cell_clustering2m), 
                    heatmapcluster = T, 
                    filename =paste0("./output/", sub_celltype, "_subclustering_heatmap_all.pdf"))


sub_df2$sample_id<-as.factor(sub_df2$sample_id)
sub_df2$cell_clustering1m<- as.factor(sub_df2$cell_clustering1m)

df<- melt(sub_df2[,c(subcluster_markers, "sample_id", "cell_clustering1m")])

pdf(paste0("./output/", sub_celltype, "_markerExpr_ridgeplot.pdf"),
    height=10, width=8)
ggplot(df[df$variable%in%subcluster_markers, ], aes(x = value, y = cell_clustering1m, fill = cell_clustering1m)) +
  geom_density_ridges()+
  facet_wrap(~variable, scales = "free", ncol = 2)
dev.off()

plot_scaled_heatmap(expr = sub_df2,
                    markers = subcluster_markers, 
                    colorassigned = colorassigned, 
                    cell_clustering1m = as.factor(sub_df2$cell_clustering1m), 
                    heatmapcluster = T, 
                    filename =paste0("./output/", sub_celltype, "_subclustering_heatmap_annotated.pdf"))



unique(output$cell_clustering2m)

# overwrite output with cellclustering_2m (specific cell type annotation)
saveRDS(output, 'backup/backup_output.rds')

df_output<- data.frame(fsApply(output$fcs1, exprs)[, c("X_position", "Y_position")], sample_id = output$sample_ids, 
                       cluster = output$cell_clustering2m, cell_clustering1m=output$cell_clustering1m)
saveRDS(df_output, './data/df_output.rds')