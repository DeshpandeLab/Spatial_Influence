#-------------------------------------#
#           plot CN cell mask         #
#-------------------------------------#


plot_CN_cellmask2 <- function(sample_id = 26, 
                              coord= coord, 
                              numClusters = 10, 
                              path="./Probability_masks/"){
  require(tiff); require(raster); require(dittoSeq);
  ## prepare data 
  exprwanted <- coord[coord$sample_id==sample_id,]
  ## load segmentation mask
  filename <- paste0(sample_id, "_Probabilities_mask.tiff")
  img<-readTIFF(paste0(path, filename))
  img_rev<-img*(2^16-1) 
  uniquetabs_rev<-table(img_rev)
  img_rast<-raster(img_rev)
  
  ## color by CN 
  ## make other celltypes all grey, color by CN colors 
  CNcolors<- dittoColors(reps = 1, get.names = FALSE)[1:numClusters]
  names(CNcolors)<- c(1:numClusters)
  CNcolors<-c(CNcolors, "#d3d3d3")
  names(CNcolors)[length(CNcolors)]<-"others"
  
  colorlist <- CNcolors[exprwanted$CN]
  
  pdf(paste0('./output/', sample_id, "_CN_cellmask.pdf"), 
      height=5, width=5) 
  par(mar=c(0,0,0,0))
  image(img_rast, 
        col = c("white", colorlist), 
        axes = FALSE, 
        main="CN clusters",
        frame.plot = FALSE)
  dev.off()
}
