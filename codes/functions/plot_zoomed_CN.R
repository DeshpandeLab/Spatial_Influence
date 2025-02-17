#-------------------------------------#
#        plot CN cropped area         #
#-------------------------------------#

plot_zoomed_CN <- function(sample_id = 26, 
                           coord= coord, 
                           colorassigned = colorassigned,
                           show = 4,
                           seed=seed,
                           range=50,
                           CNnum = 3,
                           ncols=2){
  
  pickcell<- coord[coord$sample_id==sample_id, ]
  set.seed(seed)
  pickcell_CN <- pickcell[pickcell$CN==CNnum, ][sample(nrow(pickcell[pickcell$CN==CNnum, ]), show),]  # random pick n number of ref cells

  range =range
  # Filter data for each sample_id based on the reference from pickcell
  data_for_plot_list <- lapply(1:nrow(pickcell_CN), function(i) {
    
    reference <- pickcell_CN[i, ]
    # Filter coordinates based on the range
    data_for_plot<- coord[coord$sample_id==sample_id, ]
    data_for_plot$is_reference <- FALSE 
    data_for_plot$is_reference[rownames(data_for_plot)==rownames(reference)]<-TRUE
    data_for_plot_filtered<- data_for_plot %>%
      dplyr::filter(
          X_position >= reference$X_position - range & 
          X_position <= reference$X_position + range &
          Y_position >= reference$Y_position - range & 
          Y_position <= reference$Y_position + range)%>% 
      mutate(show = i)
  })
  
  data_for_plot_combined <- bind_rows(data_for_plot_list)
  
  pdf(paste0('./output/core',sample_id, "_CN",CNnum,'.pdf'), height = 5, width=4)
  p<- ggplot(data_for_plot_combined, aes(x = X_position, y = Y_position)) + 
    geom_point(aes(color = cluster), size = 2) +
    geom_point(
      data = data_for_plot_combined %>% dplyr::filter(is_reference == TRUE),
      aes(x = X_position, y = Y_position),
      shape = 21, fill = "white", color = "black", size = 4, stroke = 1.2
    ) +
    theme_bw() +
    xlab("") + ylab("")+
    theme(strip.background = element_blank(), 
          strip.text = element_blank(), 
          aspect.ratio = 1)+
    scale_color_manual(values = colorassigned, name="cellTypes") +
    facet_wrap(~show, scales = "free", ncol=ncols)
  print(p)
  dev.off()
}
