#----------------------------------------#
#           visualize heatmap            #
#----------------------------------------#

plot_scaled_heatmap<-function(expr = expr, # include annotation 
                              markers = celltype_markers, 
                              colorassigned = colorassigned, 
                              heatmapcluster = FALSE, 
                              cell_clustering1m = factor(output$cell_clustering1m, levels=clusterlevels), 
                              filename ="./output/scaled_heatmap.pdf"
                              ){
  require(ComplexHeatmap); require(pals); 
  expr0 <-expr[,markers] # subset markers 
  expr1 <- sapply(colnames(expr0), function(x) uniform_quantile_scale(expr0[, x],
                                                                      lower_quantile = 0.33, 
                                                                      upper_quantile = 0.99))
  
  # calculate median 
  expr01_median <- data.frame(expr1, cell_clustering =cell_clustering1m, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr1)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  color_list = list(clusters=colorassigned)
  
  if (heatmapcluster){
    cp<-rowAnnotation(col=color_list,
                      gp = gpar(col = "white", lwd = .5),
                      counts= anno_barplot(
                        as.vector(table(cell_clustering1m)[rownames(expr_heat)]),
                        gp = gpar(fill=color_list$clusters),
                        border = F,
                        bar_width = 0.75, 
                        width = unit(2,"cm")))
    
    pdf(filename, height = 10, width=10)
    p<-Heatmap(as.matrix(expr_heat),
               name = "scaled",
               col = kovesi.diverging_bwr_40_95_c42(100),
               clustering_distance_rows = "euclidean",
               clustering_distance_columns = "euclidean",
               cluster_columns = T,
               cluster_rows = T,
               border = "white",
               rect_gp = gpar(col = "white", lwd = 0.5),
               right_annotation = cp,
               show_row_names = TRUE,
               row_names_gp = gpar(fontsize = 10),
               column_names_gp = gpar(fontsize = 10),
               width = ncol(expr_heat) * unit(4, "mm"), 
               height = nrow(expr_heat) * unit(4, "mm"))
    print(p)
    dev.off() 
  } else { 
    pdf(filename, height = 10, width=10)
    p<-Heatmap(as.matrix(expr_heat),
               name = "scaled",
               col = kovesi.diverging_bwr_40_95_c42(100),
               clustering_distance_rows = "euclidean",
               clustering_distance_columns = "euclidean",
               cluster_columns = F,
               cluster_rows = F,
               border = "white",
               rect_gp = gpar(col = "white", lwd = 0.5),
               show_row_names = TRUE,
               row_names_gp = gpar(fontsize = 10),
               column_names_gp = gpar(fontsize = 10),
               width = ncol(expr_heat) * unit(4, "mm"), 
               height = nrow(expr_heat) * unit(4, "mm"))
    print(p)
    dev.off() 
  }
  
  }
