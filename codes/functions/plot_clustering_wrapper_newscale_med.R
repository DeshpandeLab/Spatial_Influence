plot_clustering_heatmap_wrapper_newscale_med <- function(expr, cell_clustering, nclusters=40,
                                            color_clusters, cluster_merging = NULL, 
                                            cluster_by=output$cluster_by,
                                            clusterMergeFile=NULL,
                                            fileName = 'clusteringheatmap.pdf',lower_quantile=0.33,upper_quantile=0.99){
  require(matrixStats);require(dplyr);require(RColorBrewer);require(pheatmap);require(readxl);require(flowCore);require(scales);require(pals);require(ComplexHeatmap)
  ## Will output the heatmap object and print it 
  
  #get expression
  #expr <- fsApply(fcs, exprs);
  expr <-expr[,cluster_by]
  ## Scale expression of all markers to values between 0 and 1
  #rng <- colQuantiles(expr, probs = c(0.01, 0.99))
  expr01 <- sapply(colnames(expr), function(x) uniform_quantile_scale(expr[, x],lower_quantile = lower_quantile, upper_quantile = upper_quantile))
  expr01 <- expr01[,cluster_by]
  
  
  ## Calculate the median expression##################################################
  
  expr_median <- data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  
  ## Calculate cluster frequencies
  
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  names(clustering_prop)<- names(table(cell_clustering))
  ## Sort the cell clusters with hierarchical clustering
  
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  ## Annotation for the original clusters
  
  annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  
  ## Annotation for the merged clusters
  
  if(!is.null(clusterMergeFile)){
    ifelse(grepl(clusterMergeFile,pattern='.xlsx'),cluster_merging <- read_excel(clusterMergeFile),cluster_merging <- read.csv(clusterMergeFile,header = TRUE))
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Merged <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Merged <- color_clusters2
  }
  
  ## Colors for the heatmap
  
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  
  clusternames<-sort(unique(cell_clustering))
  
  colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(clusternames))
  
  names(colorassigned)<-clusternames
  
  color_list = list(clusters=colorassigned)
  
  color_list_byoriginal = colorassigned[match(cell_clustering,names(colorassigned))]
  
  cp<-rowAnnotation(clusters=names(colorassigned),
                    col=color_list,
                    gp = gpar(col = "white", lwd = .5),
                    prop=anno_barplot(
                      clustering_prop[rownames(expr_heat)], 
                      gp = gpar(fill=colorassigned, col=F),
                      border = F,
                      bar_width = 0.75, 
                      width = unit(2,"cm")))
  
  q <- Heatmap(expr_heat, name="scaled",
               col=rev(brewer.rdbu(100)),
               row_order = names(colorassigned),
               cluster_columns = T,
               cluster_rows = T,
               border = NA,
               rect_gp = gpar(col = "white", lwd = .5),
               right_annotation = cp,
               show_row_names = T,
               row_names_gp = gpar(fontsize=7),
               column_names_gp = gpar(fontsize=10),
               heatmap_legend_param = list(at=seq(from = 0, to = 1, by = 0.2)),
               width = unit(10, "cm"))
  
  print('Colors:')
  print(color_clusters)
  
  pdf(fileName, width=9, height=9) 
  
  grid.draw(q$gtable)
  
  dev.off() 
  
}
