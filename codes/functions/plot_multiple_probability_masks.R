#===========================================#
###      Plot multiple probability Mask        ###
#===========================================#

plot_multiple_probability_masks <- function(coord_data = coord, 
                                            uniqID = unique_id, 
                                            sampINFO = sampInfo, 
                                            probability_mask_folder = 'Probability_masks/',
                                            colorassigned = colorassigned,
                                            ncols = 5, 
                                            filenameprefix= "selected_sampID_cellmask") {
  require(raster); require(sp); require(tiff); require(pals); require(reshape2);
  require(stringr);require(scales);require(flowCore);require(ggplot2);require(gridExtra);
  
  plot_list <- list()
  
  for(sampIDs in uniqID){
    filename= paste0(probability_mask_folder, sampIDs, "_Probabilities_mask.tiff")
    if (file.exists(filename)) {
      patient_id = sampINFO[sampINFO$sample_id==sampIDs, ]$Patient
      Site = sampINFO[sampINFO$sample_id==sampIDs, ]$Site
      
      ##load segmentation mask
      img<-readTIFF(filename)
      
      img_rev<-img*(2^16-1) #16bit tiff image that numbered each cell with an integer number that goes from 0 (saved as divided by the max number)
      uniquetabs_rev<-table(img_rev)
      img_rast<-raster(img_rev)
      
      img_rast_df <- as.data.frame(rasterToPoints(img_rast))
      colnames(img_rast_df) <- c("x", "y", "value")
      
      colorlist <- colorassigned[as.character(coord_data[coord_data$sample_id==sampIDs, ]$cluster)]
      
      p<- ggplot() +
        geom_raster(data = img_rast_df, aes(x = x, y = y, fill = value)) +
        scale_fill_gradientn(colors = c("white", colorlist)) +
        coord_fixed() +  
        theme_void() +
        theme(legend.position = "none", 
              plot.title=element_text(hjust=0.5, size=8)) + 
        ggtitle(paste("Core:", sampIDs,"| Site:", Site, "| Patient:", patient_id)) 
      
      plot_list[[as.character(sampIDs)]] <- p
    }
  }
  
  size = length(uniqID)*2.5 
  
  pdf(paste0('./output/', filenameprefix, 'cellmask.pdf'), height=size, width=size)
  print(grid.arrange(grobs = plot_list, ncol = ncols)) 
  dev.off()
  
}
