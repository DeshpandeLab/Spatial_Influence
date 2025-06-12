###############################################
#  do_Cellular_Influence_(CI)_quantification  #
###############################################

do_CI_quantification <- function(expr = df_rm,
                                 sampleID = sampleID,
                                 x_col = "X_position",
                                 y_col = "Y_position",
                                 celltype_col = "celltype",
                                 sample_col = "sample_id",
                                 kernels = "gaussian",
                                 clusterlevels = clusterlevels,
                                 sigma = 10) {
  require(spatstat);
  
  spatwt<- c()
  for (i in sampleID){
    expr_k <- expr[expr[[sample_col]] == i, ]
    coord <- expr_k[, c(x_col, y_col)]
    
    rownames(coord) <- paste0(rownames(expr_k), ":", expr_k[[celltype_col]])
    
    df<- data.frame(matrix(0, nrow = nrow(coord), ncol = length(clusterlevels)))
    colnames(df) <- clusterlevels
    
    coord <- cbind(coord, df)
    
    sapply(1:nrow(coord), function(r)
      coord[r, as.character(lapply(strsplit(as.character(rownames(coord)[r]), split = ":"),"[[",2))] <<- 1)
    
    # calculation of spatial weight 
    mywin <- owin(xrange = range(coord[[x_col]]), yrange = range(coord[[y_col]]))
    spObj <- ppp(x = coord[[x_col]],
                 y = coord[[y_col]],
                 window = mywin,
                 marks = coord[, clusterlevels])
    
    # adjust sigma 
    # leaveoneout = T -> removing the self
    K <- Smooth(spObj, kernel = kernels, at = "points", leaveoneout = TRUE, sigma = sigma)
    rownames(K) <- rownames(coord)
    K <- as.data.frame(K)
    K[[sample_col]] <- i
    
    spatwt<- rbind(spatwt, K)
  } 
  
  return (spatwt) 
}
  
