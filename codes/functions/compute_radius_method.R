# -------------------------------
# compute radius method
# -------------------------------

compute_radius_method <- function(df_core, celltype_col, clusterlevels, reference_type = "Tumor", rmax_pairs = 50) {
  require(data.table)
  # df_core: data.frame/data.table for one sample_id with columns: X_position, Y_position, celltype
  DT <- as.data.frame(df_core)
  rownames(DT)<- paste0(rownames(DT), ":", DT[[celltype_col]])
  # build ppp with celltype  
  pp <- with(DT, ppp(X_position, Y_position,
                     window = owin(range(X_position), range(Y_position)),
                     marks  = DT[[celltype_col]]))
  
  # pick reference points
  m  <- as.character(marks(pp))
  ref_idx <- grepl(reference_type, m)
  ref_pp  <- pp[ref_idx]
  if (ref_pp$n == 0) return(NULL)
  
  # all pairs from reference -> all points within rmax_pairs
  xcp <- crosspairs(ref_pp, pp, rmax = rmax_pairs, what = "all")
  
  # keep <= radius_keep and > 0 (no self-self)
  keep <- xcp$d > 0
  if (!any(keep)) return(NULL)
  
  from_ids <- rownames(DT)[ref_idx][xcp$i[keep]] 
  to_ids   <- marks(pp)[xcp$j[keep]]

  # proportions by (from -> neighbor type)
  W <- data.table(from = from_ids, to_type = to_ids)[, .N, by = .(from, to_type)]
  TOT <- W[, .(total = sum(N)), by = from]
  W <- merge(W, TOT, by = "from")
  W[, prop := N / total]
  
  # wide (rows = from ids; cols = neighbor types)
  rad_wt <- data.frame(matrix(0, nrow = length(unique(W$from)), ncol = length(clusterlevels)),
                       check.names = FALSE)
  
  rownames(rad_wt)<- unique(W$from)
  names(rad_wt)<- clusterlevels
  
  r_idx <- match(W$from, rownames(rad_wt))
  c_idx <- match(W$to_type, names(rad_wt))
  
  rad_wt[cbind(r_idx, c_idx)] <- W$prop

  return(as.data.frame(rad_wt))
}
