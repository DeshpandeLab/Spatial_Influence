### post-processing of fcs file 
## remove benign cells

library(flowCore)

# core 59 
benign <- read.FCS('59_benign.fcs')
benign_expr<- as.data.frame(benign@exprs)

markerlabels <- benign@parameters@data[["name"]]
markers<- benign@parameters@data[["desc"]]
name_mapping <- setNames(markers, markerlabels)
colnames(benign_expr) <- name_mapping[colnames(benign_expr)]
index<- benign_expr$CellId # extract benign cellId 

c59<- read.FCS('59.fcs')
c59_rm <- c59
c59_rm@exprs <- c59@exprs[-index, ] # remove benign cellId

# save into new fcs file
dataDir<- "Data"
outFile <- file.path(dataDir, "c59_rmBenign.fcs")
write.FCS(c59_rm, outFile)
