##################################
#        Cluster using SOM       #
##################################

clusterfcs_dim <- function(fcs=output$fcs,
                           cluster_by = output$cluster_by,
                           xdims=35, ydims=35, 
                           seed=1234,plottitle='consensus_plots',
                           scaleoption=F,
                           numclusters=50){
  ## Cell population identification with FlowSOM and ConsensusClusterPlus
  require(dplyr);require(FlowSOM);require(ConsensusClusterPlus)
  set.seed(seed)
  som <- ReadInput(fcs, transform = FALSE, scale = scaleoption)
  som <- BuildSOM(som, colsToUse = cluster_by,xdim=xdims,ydim=ydims,distf=1,outlierMAD=2)
  
  ## Get the cell clustering into 100 SOM codes
  cell_clustering_som <- som$map$mapping[,1]
  
  ## Metaclustering into numclusters with ConsensusClusterPlus
  codes <- som$map$medianValues
  mc <- ConsensusClusterPlus(t(codes), maxK = numclusters, reps = 100,
                             pItem = 0.9, pFeature = 1, title = plottitle, 
                             plot = "png", clusterAlg = "hc", 
                             innerLinkage = "average", finalLinkage = "average",
                             distance = "euclidean", seed = 1234)
  
  ## Get cluster ids for each cell
  code_clustering <- mc[[numclusters]]$consensusClass #metaclusters consensus
  cell_clustering <- code_clustering[cell_clustering_som] #cell clustering from som
  return(list('code_clustering'=code_clustering,'cell_clustering'=cell_clustering,'metaclusters'=mc))
}


