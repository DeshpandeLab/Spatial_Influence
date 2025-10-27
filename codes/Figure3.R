#============================#
###       Figure 3         ###
#============================#

## compute FunCN spatial weight ====

# Spatial weight calculations for each core
sampleID <- unique(df_output$sample_id)

df_rm <- df_output%>% dplyr::select(-cluster)
colnames(df_rm)[4]<- "celltype"

df_rm$celltype <- factor(df_rm$celltype, levels=unique(df_output$cell_clustering1m))

# CI quantification for broader celltypes
spatwt_1m <- compute_FunCN(expr = df_rm,
                           x_col = "X_position", # x_position column name
                           y_col = "Y_position", # x_position column name
                           celltype_col = "celltype", # cell type column name
                           multisample = TRUE, # TRUE = multisamples, FALSE = single sample
                           sample_col = "sample_id", # sample id column name
                           kernels = "gaussian", # pick the kernel style 
                           sigma = 10) # adjust the sigma for the region of interest
spatwt_1m$sample_id<- factor(spatwt_1m$sample_id, levels=samplevels)


## compute KNN ====
compute_knn<- function(data, cell_col, clusterlevels, x_col, y_col, k) {
  # data frame initialization
  knn_weights<- data.frame(matrix(0, nrow = nrow(data), ncol = length(clusterlevels)))
  colnames(knn_weights) <- clusterlevels
  
  nn_indices<- data.frame(matrix(0, nrow = nrow(data), ncol = k))
  colnames(nn_indices) <- paste0("N", 1:k)
  
  # compute knn
  knn_result <- get.knn(data[, c(x_col, y_col)], k = k)
  # knn_result$nn.index: index of the cell (10 columns with increasing dist neighbors for 10-nn)
  # knn_result$nn.dist : distance of the cell from the reference cell 
  
  # add k-neighbors to the matrix
  for (i in 1:nrow(knn_result$nn.index)) {
    neighbors <- knn_result$nn.index[i, ]
    
    # nearest neighbor indices
    nn_indices[i, ] <- as.character(rownames(data)[neighbors])
    
    for (neighbor in neighbors) {
      celltype <- as.character(data[neighbor, cell_col])
      knn_weights[i, celltype] <- knn_weights[i, celltype] + 1
    }
  }
  
  # calculate fraction of neighbors 
  knn_weights<- knn_weights/k  # divide by the k number 
  rownames(knn_weights)<- rownames(data)
  rownames(nn_indices) <- rownames(data)
  
  return(list(knn_weights = knn_weights, nn_indices = nn_indices))
}


knnwt<- c()
nn_list<- c()
for (i in sampleID){   # by the core
  expr_k<- df_output[df_output$sample_id==i, ]
  dd <- expr_k[, c("cell_clustering1m", "X_position", "Y_position")]
  rownames(dd)<- paste0(rownames(expr_k), ":", expr_k$cell_clustering1m)
  
  result<- compute_knn(data = dd, 
                       cell_col = "cell_clustering1m", 
                       clusterlevels= unique(df_output$cell_clustering1m), 
                       x_col = "X_position", 
                       y_col = "Y_position", 
                       k = 10) # 10-nearest neighbors 
  K <- result$knn_weights 
  K$sample_id <- i 
  knnwt<- rbind(knnwt, K)
  
  N <- result$nn_indices 
  N$sample_id <- i
  nn_list<- rbind(nn_list, N) 
}


## Contour plot for validation ====
# choose a core to test the result 
select_core =51

# FunCN gaussian influence of T cells on Tumor 
id_1 <- which(grepl("Tumor", rownames(spatwt_1m)))

spatwt_filtered<- spatwt_1m[id_1, ]
coord_filtered <- df_output[id_1, ]

spatwt_Tcell_on_tum<- data.frame(Tcell = spatwt_filtered$CD8T+spatwt_filtered$CD4T, 
                                 sample_id = spatwt_filtered$sample_id, 
                                 X_position= coord_filtered$X_position,
                                 Y_position= coord_filtered$Y_position)

id_2 <- which(grepl("CD8T|CD4T", rownames(spatwt_1m)))

coord_filtered_2<- df_output[id_2, ]
max_range<- max(df_output$X_position, df_output$Y_position)

flipped_y <- max_range - coord_filtered_2[coord_filtered_2$sample_id == select_core, ]$Y_position

w2 <- ppp(coord_filtered_2[coord_filtered_2$sample_id==select_core, ]$X_position,
          flipped_y, # flip y-axis
          window = owin(c(0,max_range), c(0,max_range)))

densityplot <- density(w2, sigma=7)
plot(densityplot, main="T cell density")


pdf('./output/Figure3A.pdf', width=10, height=5)
par(mfrow = c(1, 2))
plot(densityplot)

# Adjust margins to make space for the legend on the right
#par(mar = c(bottom, left, top, right))
par(mar = c(6, 5, 5, 6.5))  # Add space on the right for the legend (last value = 8)

colors <- colorRampPalette(c("grey", "red"), alpha = TRUE)(100)[as.numeric(cut(spatwt_Tcell_on_tum[spatwt_Tcell_on_tum$sample_id==select_core, ]$Tcell, 100))]
plot(
  spatwt_Tcell_on_tum[spatwt_Tcell_on_tum$sample_id==select_core, ]$X_position, 
  -spatwt_Tcell_on_tum[spatwt_Tcell_on_tum$sample_id==select_core, ]$Y_position,  # flip y-axis
  col = colors,         
  pch = 20, cex = 0.6,             
  xlab = "",ylab = "",
  xaxt = "n",yaxt = "n",  
  main = "Tumors with T-cell influence"
)

# Add the legend to the right
image.plot(
  zlim = range(spatwt_Tcell_on_tum$Tcell),  
  legend.only = TRUE,               
  col = colorRampPalette(c("grey", "red"), alpha = TRUE)(100),               
  add = TRUE,                       
  legend.mar = 5,  
  legend.width = 2, 
  legend.args = list(text="", side = 3, line = 2.5, cex = 5)
)
dev.off()


## Compare scatter plot FunCN vs KNN ====
# Subset data for the specified core
spatwt_subset <- spatwt_1m %>% dplyr::filter(sample_id == select_core)
knnwt_subset <- knnwt %>% dplyr::filter(sample_id == select_core)


knnwt_subset<- knnwt_subset[, colnames(spatwt_subset)]
identical(rownames(spatwt_subset), rownames(knnwt_subset))
identical(colnames(spatwt_subset), colnames(knnwt_subset))

spatwt_subset_melted <- reshape::melt(subset(spatwt_subset, select = -sample_id))
knn_wt_subset_melted<- reshape::melt(subset(knnwt_subset, select = -sample_id))


colnames(spatwt_subset_melted)[2]<- "FunCN"
colnames(knn_wt_subset_melted)[2]<-"knn"

df<- cbind(spatwt_subset_melted, knn_wt_subset_melted[,"knn"])
names(df)<-c("neighbor", "FunCN", "knn")


## plot spatwt vs knn weights 
coi <- c("CD8T", "Myeloid", "Stroma")
df_sub <- df[df$neighbor%in%coi,]
df_sub$neighbor<- factor(df_sub$neighbor, levels=coi)


pdf('./output/Figure3B.pdf', height=2.7,width=7)
p3B<- ggplot(df_sub, aes(x= knn, y=FunCN, color=neighbor) ) + 
  geom_point() + 
  labs(x = "knn proportion", y = "FunCN influence") +
  theme_minimal()+ 
  theme(
    plot.title = element_text(size=13,color = "black"),
    axis.title.x = element_text(size=12,color = "black"),
    axis.title.y = element_text(size=12,color = "black"),
    axis.text.x = element_text(size=9,color = "black"),
    axis.text.y = element_text(size=9,color = "black"),
    strip.text = element_text(size=12,color = "black"), 
    legend.position = "none",
  )+
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~neighbor, scales="free")
print(p3B)
dev.off()


## Evaluate difference in KNN vs spatwt_1m(FunCN method)  ====
generate_neighborweights <- function(spatwt, knnwt, nn_type) {
  # data.frame of gaussian and knn weights for each cell
  wts <- lapply(
    rownames(spatwt),
    function(rowname) {
      coi <- str_replace(nn_list[rowname, nn_type], ".*:", "")
      
      data.frame(
        cellID = rowname,
        weights = c(spatwt[rowname, coi], knnwt[rowname, coi]),
        type = c("FunCN", "knn"),
        neighbor = coi,
        nn = nn_type,
        stringsAsFactors = FALSE
      )
    }
  ) %>% bind_rows()
  return(wts)
}

nn_types <- paste0("N", 1:4)
data <- nn_types %>%
  lapply(function(nn_type) generate_neighborweights(spatwt_subset, knnwt_subset, nn_type)) %>%
  bind_rows()

# Add a column for referece celltype
data <- data %>%
  mutate(
    ref_celltype = str_replace(cellID, ".*:", ""),
    across(c(nn, type, neighbor), as.factor)
  )


pdf('./output/Figure3C.pdf', height=3, width=2.5)
for (i in unique(data$neighbor)){
  p3C<- ggplot(data[data$neighbor%in%i,], aes(x = nn, y = weights, color = nn)) +
    geom_violin()  +
    geom_boxplot(width=0.2)+ 
    labs(x = "Nearest Neighbors", y = i) +
    theme_minimal()+
    theme(legend.position = "none", 
          plot.title =element_text(size=12,color = "black", hjust = 0.5) )+
    facet_wrap(~type)
  print(p3C)
}
dev.off()

## Extract the most different TME and visualize ====
identical(rownames(spatwt_1m), rownames(knnwt))
knnwt_filtered<- knnwt[id_1, ]

compare_CD8T<- as.data.frame(cbind(spatwt_filtered[,"CD8T"],knnwt_filtered[,c("CD8T", "sample_id")]))
names(compare_CD8T)<-c("spatwt", "knnwt", "sample_id")
compare_CD8T$diff <- compare_CD8T$spatwt - compare_CD8T$knnwt

# order by diff 
compare_CD8T<- compare_CD8T[order(compare_CD8T$diff, decreasing=T), ]
compare_CD8T<- compare_CD8T[!duplicated(compare_CD8T$sample_id), ]
compare_CD8T$refID<- gsub(":.*", "", rownames(compare_CD8T))

colorassigned<- c("#8c564b",#Immune_mix
                  "#1f77b4",#CD8T
                  "#ff7f0e",#CD4T
                  "#FEE480",#NK
                  "#d62728",#Myeloid
                  "#9C4DF5",#Neutrophil 
                  "#BBE5E9",#Stroma
                  "#d3d3d3",#Tumor
                  "#5c6068")#UAs
names(colorassigned)<-c("Immune_Mix", "CD8T", "CD4T", "NK","Myeloid", "Neutrophil","Stroma", "Tumor", "UA")

# pick top 5 cases where the CD8T cell influence value on Tumor 
pickcell<- compare_CD8T[1:5, ]
print(pickcell)
range =50
# Filter data for each sample_id based on the reference from pickcell
data_for_plot_list <- lapply(1:nrow(pickcell), function(i) {
  # Extract the sample ID and reference coordinates for the current pickcell
  sampID <- as.numeric(pickcell$sample_id[i])
  reference <- df_output[pickcell$refID[i],]
  
  # Filter coordinates based on the range
  data_for_plot<- df_output
  data_for_plot$is_reference <- FALSE 
  data_for_plot$is_reference[as.numeric(pickcell$refID[i])]<-TRUE
  data_for_plot_filtered<- data_for_plot %>%
    dplyr::filter(
      sample_id == sampID & 
        X_position >= reference$X_position - range & 
        X_position <= reference$X_position + range &
        Y_position >= reference$Y_position - range & 
        Y_position <= reference$Y_position + range)
})

data_for_plot_combined <- bind_rows(data_for_plot_list)

# plot Figure 3D
data_for_plot_combined$sample_id<-paste0("core:", data_for_plot_combined$sample_id)
data_for_plot_combined$sample_id<- as.factor(data_for_plot_combined$sample_id)

pdf('./output/Figure3D.pdf', height = 3, width=15)
p3D<- ggplot(data_for_plot_combined, aes(x = X_position, y = Y_position)) + 
  geom_point(aes(color = cell_clustering1m), size = 4) +
  geom_point(
    data = data_for_plot_combined %>% dplyr::filter(is_reference == TRUE),
    aes(x = X_position, y = Y_position),
    shape = 21, fill = "white", color = "black", size = 4, stroke = 1.2
  ) +
  theme_bw() +
  scale_y_reverse() + # flip y-axis to match MCD image 
  scale_color_manual(values = colorassigned, name="cellTypes") +
  facet_wrap(~sample_id, scales = "free", ncol=5) + 
  labs(x = expression("X Position (" ~ mu ~ "m)"), y = "Y Position (" ~ mu ~ "m)")
print(p3D)
dev.off()

## Compare different sigmas and different methods (FunCN vs KNN vs Radius) ====
# test different sigma ##
sigmas <- c(5, 10, 30, 50)
# KNN method 
df_core <- subset(df_output, sample_id == select_core)

knnwt_core<- c()
for (nn in sigmas) {
  
  knn_res <- compute_knn(
    data = df_core,
    cell_col = "cell_clustering1m",
    clusterlevels = unique(df_core$cell_clustering1m),
    x_col = "X_position", y_col = "Y_position",
    k = nn
  )
  K <- knn_res$knn_weights 
  K$cell_Id<- paste0(rownames(df_core), ":", df_core$cell_clustering1m)
  K$sigma <- nn
  knnwt_core<- rbind(knnwt_core, K)
}

knnwt_core$method <- "KNN"
head(knnwt_core)
reorder <- c("cell_Id",unique(df_core$cell_clustering1m) ,"sigma", "method")
knnwt_core<- knnwt_core[, reorder]
head(knnwt_core)

# FunCN method

funCN_core <- c()
for (si in sigmas) {
  funCN_res <- compute_FunCN(expr = df_core,
                             x_col = "X_position",
                             y_col = "Y_position",
                             celltype_col = "cell_clustering1m",
                             multisample = FALSE, 
                             sample_col = NA,
                             kernels = "gaussian",
                             sigma = si) 
  funCN_res$sigma <- si
  funCN_res<- funCN_res%>% rownames_to_column("cell_Id")
  funCN_core<- rbind(funCN_core, funCN_res)
}
funCN_core$method <- "FunCN"
head(funCN_core)

# radius method ---- calculate distances from tumor cells
radii <- c()
for (radius in sigmas) {
  radii_res <- compute_radius_method(df_core = df_core,
                                     celltype_col = "cell_clustering1m",
                                     clusterlevels = unique(df_output$cell_clustering1m),
                                     reference_type = "Tumor", 
                                     rmax_pairs = radius)
  # skip if NULL or has zero rows
  if (!is.null(radii_res) && NROW(radii_res) > 0) {
    radii_res$sigma <- radius
    radii_res<- radii_res%>% rownames_to_column("cell_Id")
    radii <- rbind(radii, radii_res)}
}
    
radii$method <- "Radius"
reorder <- c("cell_Id",unique(df_output$cell_clustering1m) ,"sigma", "method")
radii<- radii[, reorder]
head(radii)

## plot CD8T influence on Tumor cells using different methods ====
combined_data <- rbind(knnwt_core, funCN_core, radii)
df_core$cell_Id <- paste0(rownames(df_core), ":", df_core$cell_clustering1m)

combined_data <- combined_data %>%
  dplyr::left_join(df_core %>% dplyr::select(cell_Id, X_position, Y_position),
                   by = "cell_Id") 

cols_to_plot<- c("CD8T", "X_position", "Y_position", "sigma", "method")

## magnify zone ====
y_rng <- c(0,200)
x_rng <- c(650, 850)  

# filter x,y-coordinates
data_plot_filter <- combined_data%>% 
  dplyr::filter(between(X_position, x_rng[1], x_rng[2]),
                between(Y_position, y_rng[1], y_rng[2]))

cd8_coords_filter <- df_core %>%
  dplyr::filter(cell_clustering1m == "CD8T") %>%
  dplyr::select(X_position, Y_position) %>% 
  dplyr::filter(between(X_position, x_rng[1], x_rng[2]),
                between(Y_position, y_rng[1], y_rng[2]))


p_range<- ggplot(data_plot_filter[grep("Tumor",data_plot_filter$cell_Id), cols_to_plot], 
                 aes(x = X_position, y = Y_position)) + 
  geom_point(aes(color = CD8T), size = 6) +
  geom_point(data = cd8_coords_filter,
             aes(x = X_position, y = Y_position),
             inherit.aes = FALSE, shape = 3, 
             color = "red", alpha = 0.8, size=6) +
  theme_bw() +
  xlab("") + ylab("")+
  scale_y_reverse()+
  coord_equal()+
  scale_color_viridis_c(name = "CD8T influence", option = "viridis", limits = c(0, 1), direction = -1)+
  facet_grid(sigma~method)+ 
  theme(legend.position = "right", 
        axis.text = element_text(size=18, color="black"), 
        strip.text = element_text(size=25, face="bold"), 
        strip.background = element_blank(),
        legend.text = element_text(size = 18),      
        legend.title = element_text(size = 20, face = "bold"))

pdf('./output/FigS2C.pdf', width = 20, height=30)
print(p_range)
dev.off()

# check all CD8T cell influence on all tumor cells in core 51 
DT <- as.data.frame(df_core)
rownames(DT)<- paste0(rownames(DT), ":", DT[["cell_clustering1m"]])
# build ppp with celltype  
pp <- with(DT, ppp(X_position, Y_position,
                   window = owin(range(X_position), range(Y_position)),
                   marks  = DT[["cell_clustering1m"]]))

# pick reference points
m  <- as.character(marks(pp))
ref_idx <- grepl("Tumor", m)
ref_pp  <- pp[ref_idx]

# all pairs from reference -> all points within 30um
xcp <- crosspairs(ref_pp, pp, rmax = 30, what = "all")

# keep <= radius_keep and > 0 (no self-self)
keep <- xcp$d > 0

from_ids <- rownames(DT)[ref_idx][xcp$i[keep]] 
to_ids   <- marks(pp)[xcp$j[keep]]
dist_ids <- xcp$d[keep]
# proportions by (from -> neighbor type)
W <- data.table(from = from_ids, to_type = to_ids, dist = dist_ids)
CD8T_table <- W[W$to_type=="CD8T", ]

CD8T_table_summarized <- CD8T_table%>%
  group_by(from) %>% 
  dplyr::summarize(
    total_n = n(),
    avg_dist = mean(dist, na.rm=T))

core_51 <- combined_data[combined_data$sigma==10&combined_data$method%in%c("KNN", "FunCN")|
                           combined_data$sigma==30&combined_data$method%in%c("Radius"), ]

CD8T_table_summarized<- as.data.frame(CD8T_table_summarized)
colnames(CD8T_table_summarized)[1]<-"cell_Id"
CD8T_table_summarized<- left_join(CD8T_table_summarized, core_51[,c("cell_Id", "CD8T","method")], by="cell_Id")


pS2D_1<- ggplot(CD8T_table_summarized, aes(x=total_n, y=avg_dist, color=CD8T))+
  geom_point()+ 
  #scale_x_continuous(breaks = breaks_x, minor_breaks = NULL,
  #                   expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_reverse()+ theme_bw() + 
  scale_color_viridis_c(option = "viridis", limits = c(0, 1), direction=-1)+
  ylab("Average Distance (Tumor - CD8T)") + xlab("Number of CD8T (within 30um distance)")+
  theme(legend.position = "none", 
        axis.title = element_text(size=14, color="black"), 
        axis.text = element_text(size=12, color="black"), 
        strip.text = element_text(size=14, face="bold"), 
        strip.background = element_blank())+ 
  facet_wrap(~method)


pS2D_2 <- ggplot(CD8T_table_summarized[CD8T_table_summarized$total_n%in%c(1,2,3), ], 
       aes(x=avg_dist, y=CD8T, color=method))+
  geom_col()+ theme_bw() + 
  xlab("Average Distance (Tumor - CD8T)") + ylab("CD8T influence")+
  theme(legend.position = "none", 
        axis.title = element_text(size=14, color="black"), 
        axis.text = element_text(size=10, color="black"), 
        strip.text = element_text(size=14, face="bold"), 
        strip.background = element_blank())+ 
  facet_grid(total_n~method)

pdf('./output/FigureS2D_1.pdf', height=6.5, width=5)
print(pS2D_1)
dev.off()

pdf('./output/FigureS2D_2.pdf', height=4, width=5.5)
print(pS2D_2)
dev.off()


## Apply FunCN & compare TME between PRI pancr and MET liver ====
# Dot plot for pairwise Cell-cell interaction in Pancreas & Liver
df_rm_2m <- df_output%>% dplyr::select(-cell_clustering1m)
colnames(df_rm_2m)[4]<- "celltype"

df_rm_2m$celltype <- factor(df_rm_2m$celltype, levels=clusterlevels)

spatwt <- compute_FunCN(expr = df_rm_2m,  # should contain sample_id, X, Y position, celltype in the df
                        x_col = "X_position",
                        y_col = "Y_position",
                        celltype_col = "celltype",
                        multisample = TRUE,
                        sample_col = "sample_id", 
                        kernels = "gaussian",
                        sigma = 10) 
spatwt$sample_id<- factor(spatwt$sample_id, levels=samplevels)
# save spatwt 
# saveRDS(spatwt, "./backup/spatwt.rds")

df_cci <- spatwt%>% rownames_to_column("ref")
df_cci$ref <- gsub(".*:", "", df_cci$ref)
df_cci$ref<- as.factor(df_cci$ref)
df_cci<- left_join(df_cci, Specimen_designation[,c("sample_id", "Site")], by="sample_id")
df_cci<- reshape::melt(df_cci)
colnames(df_cci)[4]<-"influence"

# scale
df_cci_summary <- df_cci %>%
  dplyr::group_by(ref, influence, Site) %>%
  dplyr::summarise(mean_value = mean(value), .groups = "drop") %>%
  dplyr::filter(ref != influence) %>% 
  dplyr::group_by(Site, influence)%>%
  dplyr::mutate(scaled = as.numeric(scale(mean_value)))%>% 
  ungroup()


df_cci_summary$ref<- factor(df_cci_summary$ref, levels=rev(clusterlevels))
df_cci_summary$influence<- factor(df_cci_summary$influence, levels=clusterlevels)
df_cci_summary$Site<- factor(df_cci_summary$Site, levels=sitelevels)

df_cci_summary<- df_cci_summary%>% 
  dplyr::filter(mean_value>0.01) # removing pairwise interaction mean_value < 0.01

## filter 
pdf('./output/FigureS3A.pdf', height=8, width=12)
p2_1<-ggplot(df_cci_summary, aes(x=influence, y=ref)) + 
  geom_point(pch = 21, stroke = 0.5, col="grey", aes(size = mean_value, fill = scaled)) + 
  theme_classic() +
  ylab("") + 
  xlab("")+
  theme(axis.line = element_blank(), 
        axis.text= element_text(colour ="black", size = 12),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  
        panel.grid.major = element_line(size = 0.5, linetype = 'solid'),  
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid'),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7), 
        legend.key.size = unit(0.3, "cm"),
        legend.position="left") + 
  #scale_color_discrete("grey")+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       name="Relative scale\nfor reference\n(scaled down\nthe column)")+
  scale_y_discrete(position = "right") + 
  scale_x_discrete(
    guide = guide_axis(angle = 45),position = "top") + 
  scale_size_continuous(breaks = seq(0, max(df_cci_summary$mean_value), by = 0.1))+
  facet_wrap(~Site)
print(p2_1)
dev.off()


## overall CD8T+ TME comparsion between Pancreas and Liver (radial plot) ====
# Define the index
idx <- list(
  which(df_output$cluster == "CD8T" & df_output$sample_id != 43),
  which(df_output$cluster == "Tumor" & df_output$sample_id != 43)
)

# sample info dataframe
unique_id <- unique(df_output$sample_ids)
sampInfo <- data.frame(sample_id = unique_id)
sampInfo$Site <- as.factor(df_output$Site[match(sampInfo$sample_id, df_output$sample_id)])
sampInfo$Patient <- as.factor(df_output$Patient[match(sampInfo$sample_id, df_output$sample_id)])

majorCelltypes <- list(c("CD8T", "CD4T", "Treg", "M_I", "M_II", "Neutrophil", "Str_I", "Str_II", "Tumor"), 
                       c("CD8T","CD4T","Treg","NK","M_I", "M_II", "Neutrophil", "Str_I", "Str_II", "UA"))
cellTypes<- c("CD8T", "Tumor")

tme_colors <- c(
  "#1f77b4",  # CD8T
  "#ff7f0e",  # CD4T
  "#ff9f0e",  # Treg
  "#FEE480",  # NK
  "#2ca02c",  # M_I
  "#d62728",  # M_II
  "#9C4DF5",  # Neutrophil
  "#BBE5E9",  # Str_I
  "#0099B5",  # Str_II
  "#d3d3d3",  # Tumor
  "#5c6068"   # UA
)
names(tme_colors) <- c("CD8T", "CD4T","Treg","NK", "M_I", "M_II", "Neutrophil", "Str_I", "Str_II", "Tumor", "UA")

### TME comparison between Pancreas and Liver (CD8T and Tumor)
cell_types <- c("Tumor", "CD8T")

all_results <- list()

for (cell_type in cell_types) {
  idx <- which(df_output$cluster == cell_type & 
                 df_output$sample_id != 43)
  
  # Subset data
  spatwt_df <- spatwt[idx, ]
  spatwt_df$Site <- factor(df_output$Site[match(spatwt_df$sample_id, df_output$sample_id)], levels = sitelevels)
  spatwt_df$Patient <- factor(df_output$Patient[match(spatwt_df$sample_id, df_output$sample_id)], levels = unique(df_output$Patient))
  
  # Filter data
  filtered_df <- spatwt_df[, clusterlevels]
  
  # Identify pancreas and liver indices
  PanIdx <- which(spatwt_df$Site == "Pancreas")
  LivIdx <- which(spatwt_df$Site == "Liver")
  
  # Initialize results dataframe
  results <- data.frame(CellType = character(), 
                        PValue = numeric(), 
                        AdjustedPValue = numeric(), 
                        MeanPancreas = numeric(),
                        MeanLiver = numeric(),
                        log2FC = numeric())
  
  for (cell in clusterlevels) {
    pancreas_vals <- filtered_df[PanIdx, cell]
    liver_vals <- filtered_df[LivIdx, cell]
    
    # Run Wilcoxon test
    test_result <- wilcox.test(pancreas_vals, liver_vals, exact = FALSE)
    
    # Calculate means
    mean_pancreas <- mean(pancreas_vals, na.rm = TRUE)
    mean_liver <- mean(liver_vals, na.rm = TRUE)
    log2fc <- log2(mean_liver / mean_pancreas)
    
    results <- rbind(results, data.frame(CellType = cell,
                                         PValue = test_result$p.value,
                                         AdjustedPValue = NA,
                                         MeanPancreas = mean_pancreas,
                                         MeanLiver = mean_liver,
                                         log2FC = log2fc))
  }
  
  # Adjust p-values using Bonferroni correction
  results$AdjustedPValue <- p.adjust(results$PValue, method = "bonferroni")
  results$logP <- -log10(results$AdjustedPValue)
  
  # Replace Inf with maximum finite value
  results$logP[is.infinite(results$logP)] <- max(results$logP[!is.infinite(results$logP)], na.rm = TRUE)
  
  all_results[[cell_type]] <- results
  
  pdf(paste0('./output/FigureS3B_', cell_type, '.pdf'), height = 5, width = 5)
  pS3B<- ggplot(results, aes(CellType, log2FC, fill = logP)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    scale_y_continuous(breaks = scales::breaks_width(1)) +
    scale_fill_gradient(low = "yellow", high = "red", name = "-log10(adj.pvalue)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    coord_polar()
  print(pS3B)
  dev.off()
}