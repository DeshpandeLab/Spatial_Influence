plot_marker_level <- function(coord_data = coord, 
                              expression_data= new_expr, 
                              uniqID = unique_sampleID,
                              marker ="PDL1",
                              filenameprefix="Tumor_PDL1", 
                              height_size = 5, 
                              width_size =18,
                              ncols=5) {
  require(ggplot2); require(dplyr); require(pheatmap); require(scales); 
  
  
  data<- cbind(coord_data, expression_data[, marker])
  colnames(data)[length(data)]<-marker # name the last column (added from expr data) -> marker name 
  data_filtered <- data[data$sample_id%in% uniqID, ]
  
  pdf(paste0('./output/', filenameprefix, '_level.pdf'), height=height_size, width=width_size)
  p1<-ggplot(data_filtered,aes(x=X_position, y=Y_position, color=get(marker))) + 
    geom_point(size=0.5) + 
    scale_color_distiller(palette = "Spectral", direction = -1)+
    theme_minimal()+
    theme(
      panel.grid = element_blank(),         
      axis.text = element_blank(),          
      axis.ticks = element_blank(),         
      axis.title = element_blank(),         
    )+
    scale_y_reverse()+ 
    facet_wrap(~sample_id, ncol=ncols)+
    labs(color = marker)
  print(p1)
  dev.off()
  
  data_filtered <- data_filtered %>%
    group_by(sample_id) %>%
    mutate(scaled = scales::rescale(get(marker), to = c(0, 1))) %>% # Scale to 0-1 within each group
    ungroup()
  
  pdf(paste0('./output/', filenameprefix, '_scaled.pdf'), height=height_size, width=width_size)
  p2<- ggplot(data_filtered,aes(x=X_position, y=Y_position, color=scaled)) + 
    geom_point(size=0.5) + 
    scale_color_distiller(palette = "Spectral", direction = -1)+
    theme_minimal()+
    theme(
      panel.grid = element_blank(),         
      axis.text = element_blank(),          
      axis.ticks = element_blank(),         
      axis.title = element_blank(),         
    )+
    scale_y_reverse()+ 
    facet_wrap(~sample_id, , ncol=ncols)+ 
    labs(color = paste(marker, "scaled"))
  print(p2)
  dev.off()
}

