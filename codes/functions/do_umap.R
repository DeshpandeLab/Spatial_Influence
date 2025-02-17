#######################################
#             do UMAP                 #
#######################################

# umap with sampling 

#separate UMAP also created in Giotto
do_umap <- function(expr,subtype_markers,sample_ids,cell_clustering,
                    clusterMergeFile,
                    seed = 1234, ncells=5000,sample_subset=NULL){
  require(umap);require(flowCore);require(readxl)
  expr <-expr[,subtype_markers]
  
  ## Create vector to later find and skip duplicates
  dups <- duplicated(expr[, subtype_markers])
  dups <- which(!(dups))## Find and skip duplicates
  
  
  ## Data subsampling: create indices by sample
  inds <- split(1:length(sample_ids), sample_ids) #to get original indexes belonging to each cluster
  samplenames <- names(inds) #create a name vector of the files
  custom.settings = umap.defaults
  custom.settings$seed = seed
  
  ifelse(is.null(sample_subset),
         umap_ncells <- pmin(table(sample_ids), ncells),
         umap_ncells <- pmin(table(sample_ids), ncells)[sample_subset]
  )
  if(!is.null(sample_subset)){inds <- inds[sample_subset]}
  umap_inds <- lapply(names(inds), function(i){
    s <- sample(inds[[i]], umap_ncells[i], replace = FALSE)
    intersect(s, dups)
  })
  set.seed(seed)
  umap_inds <- unlist(umap_inds)
  umap_out <- umap(expr[umap_inds, subtype_markers], config = custom.settings, method = 'naive')
  umapRes2D = data.frame(umap1 = umap_out$layout[, 1], umap2 = umap_out$layout[, 2], 
                         expr[umap_inds, subtype_markers],
                         sample_id = sample_ids[umap_inds], cell_clustering = factor(cell_clustering[umap_inds]), check.names = FALSE)
  return(umapRes2D)
}
