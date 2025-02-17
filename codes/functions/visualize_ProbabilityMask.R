#--------------------------------#
#           cell Mask            #
#--------------------------------#


visualize_ProbabilityMask<- function(expr0 = expr0, # should contain sampleID, X,Y position, celltype annotation 
                                     samp_id = 26, 
                                     colorassigned = colorassigned
                                     )
  {
  require(raster); require(sp); require(tiff); require(pals); require(reshape2);
  require(stringr); require(scales); require(flowCore); 
  filename= paste0('Probability_masks/', samp_id, "_Probabilities_mask.tiff")

  ##load segmentation mask
  img<-readTIFF(filename)
    
  img_rev<-img*(2^16-1) #16bit tiff image that numbered each cell with an integer number that goes from 0 (saved as divided by the max number)
  uniquetabs_rev<-table(img_rev)
  img_rast<-raster(img_rev)
    
  img_rast_df <- as.data.frame(rasterToPoints(img_rast))
  colnames(img_rast_df) <- c("x", "y", "value")
    
  colorlist <- colorassigned[as.character(expr0[expr0$sample_id==samp_id, ]$cluster)]
    
  pdf(paste0('./output/', samp_id, "_cellmask.pdf"), width=5, height=5)
  p<- ggplot() +
      geom_raster(data = img_rast_df, aes(x = x, y = y, fill = value)) +
      scale_fill_gradientn(colors = c("white", colorlist)) +
      coord_fixed() +  
      theme_void() +
      theme(legend.position = "none") 
  print(p)
  dev.off()
  }

  
