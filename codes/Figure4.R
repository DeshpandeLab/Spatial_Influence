#============================#
###       Figure 4         ###
#============================#

## Tumor TME analysis ====
spatwt<- readRDS('./backup/spatwt.rds')

# filter all tumors (remove core 43)
index <- which(output$cell_clustering2m=="Tumor"&
                 output$sample_ids!=43)

spatwt_df_tum <- spatwt[index, ]

# subset bulk tumor & tumor boundary 

cutoff = 0.7
subset <- which(spatwt_df_tum$Tumor<cutoff)


# visualize cellmask 
colorassigned<- c("#8c564b",#Immune_mix
                  "#1f77b4",#CD8T
                  "#ff7f0e",#CD4T
                  "#ff9f0e",#Treg
                  "#FEE480",#NK
                  rep("#d62728",6), 
                  "#9C4DF5",#Neutrophil 
                  rep("#BBE5E9",7), 
                  "#d3d3d3",#Tumor
                  "#5c6068")#UAs
names(colorassigned)<-clusterlevels

mask_data <- data.frame(sample_id = output$sample_ids, cluster = output$cell_clustering2m)
visualize_ProbabilityMask(expr0 = mask_data, 
                          samp_id = 18, 
                          colorassigned = colorassigned)
visualize_ProbabilityMask(expr0 = mask_data, 
                          samp_id = 25, 
                          colorassigned = colorassigned)

# visualize tumor classification
coord <- fsApply(output$fcs1, exprs)[, c("X_position", "Y_position")]
coord <-as.data.frame(coord)
coord$sample_id <- output$sample_ids
coord$cluster <- output$cell_clustering2m
coord$Y_position<- 1000-coord$Y_position


dd <- coord[index,] 
dd$type<- "BK"
dd[subset, "type"]<- "BD"
dd$type<- as.factor(dd$type)

pdf('./output/Figure4B.pdf', height=6, width=3)
p4B<-ggplot(dd[dd$sample_id%in%c("18", "25"),], aes(x=X_position, y=Y_position, color=type))+
  geom_point(size = 0.5)+ 
  scale_color_manual(values = c("BD" = "black", "BK" = "#DF536C")) + 
  theme_void() + 
  theme(legend.position = "none",)+
  facet_wrap(~sample_id, ncol=1)
print(p4B)
dev.off()

# proportion of CNs 
expr_tum<- as.data.frame(new_expr[index,])
expr_tum$sample_id <- factor(output$sample_ids[index], levels=samplevels)
expr_tum$Patient <-factor(output$Patient[match(expr_tum$sample_id,  output$meta_data$sample_id)], levels=unique(output$Patient))
expr_tum$Site <- factor(output$Site[match(expr_tum$sample_id,  output$meta_data$sample_id)], levels=sitelevels)

tum_markers<- c("KI67", "CK", "CD86", "VISTA", "PDL1")
markerExpr <- expr_tum[, c(tum_markers, "Site", "sample_id")]
markerExpr$type <- "BK"
markerExpr[subset, "type"]<- "BD"


CNprop <- as.data.frame(table(markerExpr$type, markerExpr$Site))
names(CNprop)<-c("CN", "Site", "Freq")
CNprop<- CNprop%>% 
  group_by(Site)%>% 
  mutate(total_counts = sum(Freq), 
         prop=Freq/total_counts*100)

CNprop_1 <- as.data.frame(CNprop)
CNprop_1$CN <- as.character(CNprop_1$CN)
CNprop_1$CN[CNprop_1$CN!="BK"]<-"BD"

CNprop_1 <- CNprop_1 %>% 
  group_by(Site, CN) %>% 
  dplyr::summarise(counts = sum(Freq), .groups = "drop") %>%  # Calculate counts per Site and CN
  group_by(Site) %>% 
  mutate(
    total_counts = sum(counts),    # Total counts per Site
    prop = (counts / total_counts) * 100  # Proportion in percentage
  ) %>% 
  ungroup()


pdf('./output/Figure4C.pdf', height=4, width=2.5)
p4C<- ggplot(CNprop_1, aes(x=Site, y=prop, fill=CN))+
  scale_fill_manual(values = c("BK" = "#DF536C", "BD" = "black"))+
  geom_col(position="fill", width = 0.7) +
  theme_minimal()+
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x * 100)) +
  labs(title = "", x = "", y = "", fill="Tumor\nClassification")+ 
  theme(axis.text=element_text(size=14,color = "black"),
        axis.text.x = element_text(angle=45, hjust =1),
        legend.position="right", 
        legend.text = element_text(size=13)) 
print(p4C)
dev.off()



#====

data_for_violin <- melt(cbind(spatwt_df_tum[, clusterlevels], markerExpr))
data_for_violin <- data_for_violin %>%
  dplyr::mutate(type = ifelse(type == "BK", "BK", "BD"))

summary_df_2 <- data_for_violin %>%
  group_by(variable, sample_id,type, Site)%>% 
  dplyr::summarise(
    mean_value = mean(value))


selected_markers <- c("CK","CD86")


# Tumor Marker expression level between Sites
pdf('./output/Figure4D.pdf', height=2.5, width=3)
p4D<- ggplot(summary_df_2[summary_df_2$variable%in%selected_markers, ], aes(x = type, y = mean_value, fill = Site)) +
  geom_violin(alpha = 0.6, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 0.8), alpha = 0.5) +
  stat_compare_means(
    method = "wilcox.test",  
    paired = F, 
    label = "p.signif", 
    hide.ns = T, 
    size=5
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  facet_wrap(~variable, scales = "free", ncol = 2) +  
  theme_bw() + 
  ylab("Mean Marker Expression") +
  xlab("Tumor type")+ 
  theme(legend.position = "top", 
        axis.text = element_text(size=12, color="black"), 
        axis.title=element_text(size=12, color="black"))
print(p4D)
dev.off()


selected_celltypes <- c("Treg","M_I", "Str_I", "Str_II","Tumor")


# Cellular influences on Tumors (Compare between Sites)
pdf('./output/Figure4E.pdf', height=2.5, width=6)
p4E<- ggplot(summary_df_2[summary_df_2$variable%in%selected_celltypes, ], aes(x = type, y = mean_value, fill = Site)) +
  geom_violin(alpha = 0.6, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 0.8), alpha = 0.5) +
  stat_compare_means(
    method = "wilcox.test",  
    paired = F, 
    label = "p.signif", 
    hide.ns = T, 
    size=5
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  facet_wrap(~variable, scales = "free", ncol = 5) +  
  theme_bw() + 
  ylab("Mean Cellular Influences") +
  xlab("Tumor type")+ 
  theme(legend.position = "top", 
        axis.text = element_text(size=12, color="black"), 
        axis.title=element_text(size=12, color="black"))
print(p4E)
dev.off()

pdf('./output/tumor_BK.pdf', height=10, width=10)
p<- ggplot(summary_df_2[summary_df_2$type=="BK", ], aes(x = Site, y = mean_value, color = Site)) +
  geom_violin(alpha = 0.6, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 0.8), alpha = 0.5) +
  stat_compare_means(
    method = "wilcox.test",  
    comparisons = site_comparison,  
    label = "p.format") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  facet_wrap(~variable, scales = "free", ncol = 5) +  
  theme_bw() 
print(p)
dev.off()



## spearman correlation between weight vs markers ==== 

data_for_heatmap <- melt(cbind(spatwt_df_tum[, clusterlevels], markerExpr))
data_for_heatmap <- data_for_heatmap %>%
  dplyr::mutate(type = ifelse(type == "BK", "BK", "BD"))


selected_celltypes <- c("Immune_Mix","CD8T","CD4T","Treg", "NK", 
                        "M_I","M_II","M_III","Neutrophil", "Str_I", "Str_II", 
                        "Str_III", "Str_IV", "Tumor", "UA")

selected_markers<- c("CK", "KI67", "PDL1", "CD86", "VISTA")

summary_df <- data_for_heatmap %>%
  group_by(variable, sample_id,type, Site)%>% 
  dplyr::summarise(
    mean = mean(value))

summary_df_selected<- summary_df%>% 
  dplyr::filter(variable%in%selected_celltypes)

data_melted<- summary_df %>%
  dplyr::filter(variable %in% selected_markers) %>%  # Keep only selected variables
  pivot_wider(
    names_from = variable,  # Use the selected variables as new column names
    values_from = mean      # Populate the new columns with 'mean'
  )


data_for_plot_2<- left_join(summary_df_selected, data_melted, by=c("sample_id", "type", "Site"))

pdf('./output/CK_spearman_BK.pdf', height=10, width=10)
p2<- ggplot(data_for_plot_2[data_for_plot_2$type=="BK", ], 
           aes(x=CK, y=mean, color=Site)) +
  geom_point(alpha=0.5) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))+
  stat_cor(method='spearman', 
           label.y.npc="top", 
           label.x.npc = "left",
           label.sep='\n',
           size=4) + 
  facet_wrap(~variable, scales = 'free') + 
  geom_smooth(method=lm, aes(fill=Site), alpha=0.2)+
  ylab('Mean Cellular Influences') +
  xlab("Mean CK Expression")
print(p2)
dev.off()


pdf('./output/CK_spearman_BD.pdf', height=10, width=10)
p3<- ggplot(data_for_plot_2[data_for_plot_2$type=="BD", ], 
       aes(x=CK, y=mean, color=Site)) +
  geom_point(alpha=0.5) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))+
  stat_cor(method='spearman', 
           label.y.npc="top", 
           label.x.npc = "left",
           label.sep='\n',
           size=4) + 
  facet_wrap(~variable, scales = 'free') + 
  geom_smooth(method=lm, aes(fill=Site), alpha=0.2)+
  ylab('Mean Cellular Influences') +
  xlab("Mean CK Expression")
print(p3)
dev.off()



# table for spearman correlation 

# Function to calculate p-values and rho values
calculate_stats <- function(data, selected_markers) {
  p_res <- list()
  rho_res <- list()
  
  for (i in selected_markers) {
    for (site in unique(data$Site)) {
      site_data <- data[data$Site == site, ]
      
      # Calculate p-values and rho for each marker
      stats <- site_data %>%
        dplyr::group_by(variable) %>%
        dplyr::summarise(
          p_value = cor.test(.data[[i]], mean, method = "spearman")$p.value,
          rho = cor(.data[[i]], mean, method = "spearman"),
          .groups = "drop"
        ) %>%
        dplyr::mutate(marker = i, Site = site)
      
      # Append results to lists
      p_res[[paste(i, site, sep = "_")]] <- stats %>%
        dplyr::select(variable, marker, Site, p_value)
      rho_res[[paste(i, site, sep = "_")]] <- stats %>%
        dplyr::select(variable, marker, Site, rho)
    }
  }
  
  # Combine results
  p_combined <- do.call(rbind, p_res)
  rho_combined <- do.call(rbind, rho_res)
  
  return(list(p_values = p_combined, rho_values = rho_combined))
}


selected_markers <- c("CK", "KI67", "PDL1", "CD86", "VISTA")  # Replace with your markers
data_types <- c("BK", "BD")
sites <- c("Pancreas", "Liver")


selected_markers <- c("CK", "PDL1", "VISTA", "CD86", "KI67")  # Replace with your markers

for (data_type in data_types) {
    type_data <- data_for_plot_2[data_for_plot_2$type == data_type, ]
    
    # Calculate stats
    stats <- calculate_stats(type_data, selected_markers)
    
    rho_data = stats$rho_values
    p_data = stats$p_values
    Site = stats$site
    data_type = data_type
    
    
    rho_data$label<- paste0(rho_data$Site, "_", rho_data$marker)
    rho_filtered <- rho_data %>%
      dplyr::select(variable, label, rho)%>%
      tidyr::pivot_wider(names_from = label, values_from = rho)
    
    p_data$label<- paste0(p_data$Site, "_", p_data$marker)
    p_filtered <- p_data %>%
      dplyr::select(variable, label, p_value)%>%
      tidyr::pivot_wider(names_from = label, values_from = p_value)
    
    # Replace NA values with 0 for rho and 1 for p-values to avoid issues
    rho_filtered[is.na(rho_filtered)] <- 0
    p_filtered[is.na(p_filtered)] <- 1
    
    # Extract row and column names
    rho_filtered<- rho_filtered%>% column_to_rownames("variable")
    p_filtered<- p_filtered%>% column_to_rownames("variable")
    
    
    # Convert to matrix
    rho_matrix <- as.matrix(rho_filtered)  # Exclude the variable column
    
    p_matrix <- as.matrix(p_filtered)  # Exclude the variable column
    p_matrix<- t(p_matrix)
    
    cell_function <- function(j, i, x, y, width, height, fill) {
      # Define significance level
      significance <- ""
      if (p_matrix[i, j] < 0.0001) {
        significance <- "****"
      } else if (p_matrix[i, j] < 0.001) {
        significance <- "***"
      } else if (p_matrix[i, j] < 0.01) {
        significance <- "**"
      } else if (p_matrix[i, j] < 0.05) {
        significance <- "*"
      }
      
      # Add asterisks to the cell
      if (significance != "") {
        grid.text(
          significance,
          x, y,
          gp = gpar(fontsize = 10, col = "black")
        )
      }
    }
    
    site_colors <- setNames(
      ifelse(grepl("Liver", colnames(rho_matrix)), "#00BFC4", "#F8766D"),  
      colnames(rho_matrix)  
    )
    
    
    site_annotation <- rowAnnotation(
      Site = colnames(rho_matrix),
      col = list(Site = site_colors),
      simple_anno_size = unit(0.2, "cm")
    )

    pdf(paste0('./output/Figure4F_',data_type, '_tumor.pdf'), height=6, width=7)
    ht <- Heatmap(
      as.matrix(t(rho_matrix)),
      name = "Spearman\nRho",
      col = colorRamp2(c(-1, 0, 1), c("#1f77b4", "white", "#ff7f0e")),  
      cluster_columns = T, 
      cell_fun = cell_function,
      row_names_side = "right", 
      column_names_side = "bottom", 
      column_names_rot = 45, 
      width = unit(nrow(rho_matrix) * 0.6, "cm"),  
      height = unit(ncol(rho_matrix) * 0.4, "cm"),
      column_title = paste0(data_type," Tumor"), 
      rect_gp = gpar(col = "black", lwd = 1.5), 
      right_annotation = site_annotation
    )
    draw(ht)
    dev.off()
  }



## Representative Images ==== 
select_sampID<-  c("17", "20", "53", "56")

colorassigned<- c("#8c564b",#Immune_mix
                  "#1f77b4",#CD8T
                  "#ff7f0e",#CD4T
                  "#2ca02c",#Treg
                  "#FEE480",#NK
                  rep("#d62728",6), 
                  "#9C4DF5",#Neutrophil 
                  rep("#BBE5E9",7), 
                  "#d3d3d3",#Tumor
                  "#5c6068")#UAs

names(colorassigned)<-clusterlevels

unique_id<- unique(output$sample_ids)


sampInfo <- data.frame(sample_id = unique_id)
sampInfo$Site<- as.factor(output$Site[match(sampInfo$sample_id,output$meta_data$sample_id)])
sampInfo$Patient <- as.factor(output$Patient[match(sampInfo$sample_id,output$meta_data$sample_id)])

mask_data <- data.frame(sample_id = output$sample_ids, cluster = output$cell_clustering2m)
plot_multiple_probability_masks(coord_data = mask_data,
                                uniqID =select_sampID,
                                sampINFO = sampInfo,
                                probability_mask_folder = 'Probability_masks/',
                                colorassigned = colorassigned,
                                ncols = 4, 
                                filenameprefix= "cellmask")


index <- which(output$cell_clustering2m=="Tumor"&
                     output$sample_ids%in%select_sampID)

plot_marker_level(coord_data = coord[index, ],
                  expression_data= new_expr[index, ], 
                  uniqID = select_sampID,
                  marker ="CK", 
                  filenameprefix="Tumor_CK",
                  ncols=4)

# Supplementary Figure 4
## Matched Patients ==== 
# matched patients 
matchedPA<- intersect(unique(output$Patient[output$Site=="Pancreas"]), 
                      unique(output$Patient[output$Site=="Liver"]))

matched_sampID<- unique(output$meta_data$sample_id[output$Patient%in%matchedPA&
                                                     output$meta_data$sample_id!=43])

## CD8T/Treg ratio 
Tcell_ratio <-as.data.frame.matrix(table(output$sample_ids, output$cell_clustering2m)[,c("CD8T", "Treg")])
Tcell_ratio <- Tcell_ratio%>% rownames_to_column("sample_id")
Tcell_ratio$Patient <-  factor(output$Patient[match(Tcell_ratio$sample_id,output$meta_data$sample_id)], levels=unique(output$Patient))
Tcell_ratio$Site <-  factor(output$Site[match(Tcell_ratio$sample_id,output$meta_data$sample_id)], levels=sitelevels)

Tcell_ratio$ratio <- Tcell_ratio$CD8T/ Tcell_ratio$Treg

# remove when both CD8T and Treg count is zero
Tcell_ratio<- Tcell_ratio[!is.na(Tcell_ratio$ratio),]

# select matched sample ids 
Tcell_ratio_PA <- Tcell_ratio[Tcell_ratio$sample_id%in%matched_sampID, ]
site_comparison<- list(c("PRI", "MET"))

Tcell_ratio_PA$Site <- as.character(Tcell_ratio_PA$Site)

Tcell_ratio_PA[Tcell_ratio_PA$Site == "Pancreas", "Site"] <- "PRI"
Tcell_ratio_PA[Tcell_ratio_PA$Site == "Liver", "Site"] <- "MET"
Tcell_ratio_PA$Site <- factor(Tcell_ratio_PA$Site, levels=c("PRI", "MET"))


pdf('./output/FigureS4A.pdf', height=3, width = 3)
pS4A<- ggplot(Tcell_ratio_PA, aes(x = Site, y = ratio, fill = Site)) +
  geom_violin(alpha=0.6)+
  geom_boxplot(width=0.3, fill = NA, color = "black", alpha = 0.5)+
  geom_point(aes(shape = Patient),  # Map Patient to shape
             position = position_jitter(width = 0.2, height = 0), # Add jitter for clarity
             size = 3, alpha = 0.8) +  # Adjust size and alpha for visibility
  stat_compare_means(
    method = "wilcox.test",
    comparisons = site_comparison,
    paired = FALSE, 
    hide.ns = FALSE,
    label="p.format", 
    size=5)+ 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+ 
  theme_bw()+ 
  xlab("") + ylab("CD8+T/Treg count ratio")+
  theme(axis.text.x = element_text(size=14, angle=45,  hjust =1), 
        axis.text = element_text(size=14, color="black"),
        axis.title = element_text(size=14, color="black"))
print(pS4A)
dev.off()


# Cellular Influence on Tumors in matched patients 
idx_tum <- which(output$cell_clustering2m=="Tumor"& 
                   output$Patient%in%matchedPA& 
                   output$sample_ids!=43)

idx_cd8T<- which(output$cell_clustering2m=="CD8T"& 
                   output$Patient%in%matchedPA& 
                   output$sample_ids!=43)


# Define the list of `spatwt` objects 
spatwt_list <- list(spatwt[idx_tum,] , spatwt[idx_cd8T,])  
output_names <- c("tumor", "cd8T")
celltypes<- list(c("M_I", "M_II", "M_III",
                   "Str_I", "Str_II", "Tumor", "Neutrophil", "UA"), 
                 c("CD8T","CD4T", "Treg", "M_I", "M_II", "Str_I", "Str_II", "Tumor"))  

for (idx in seq_along(spatwt_list)) {
  spatwt_PA <- spatwt_list[[idx]]  
  celltype <- output_names[idx]  
  majorCelltypes<- celltypes[[idx]]
  
  spatwt_PA$Site <- factor(output$Site[match(spatwt_PA$sample_id, output$meta_data$sample_id)], levels = sitelevels)
  spatwt_PA$Patient <- factor(output$Patient[match(spatwt_PA$sample_id, output$meta_data$sample_id)], levels = unique(output$Patient))
  
  # Initialize empty matrices
  spatwt_pan <- c()
  spatwt_liv <- c()
  
  # Calculate means for Pancreas
  for (i in matchedPA) {
    tmp <- apply(spatwt_PA[spatwt_PA$Site == "Pancreas" &
                             spatwt_PA$Patient == i, clusterlevels], 2, mean)
    spatwt_pan <- cbind(spatwt_pan, tmp)
  }
  colnames(spatwt_pan) <- paste0(matchedPA, "_P")
  
  # Calculate means for Liver
  for (i in matchedPA) {
    tmp <- apply(spatwt_PA[spatwt_PA$Site == "Liver" &
                             spatwt_PA$Patient == i, clusterlevels], 2, mean)
    spatwt_liv <- cbind(spatwt_liv, tmp)
  }
  colnames(spatwt_liv) <- paste0(matchedPA, "_L")
  
  # Compute log2 fold-change
  FC_result <- log2(spatwt_liv / spatwt_pan)
  FC_result[FC_result == 'Inf'] <- NA
  FC_result[FC_result == '-Inf'] <- NA
  FC_result[is.na(FC_result)] <- NA
  colnames(FC_result) <- matchedPA
  
  # Define palette and breaks
  my_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  max_abs <- max(abs(FC_result[majorCelltypes, ]), na.rm = TRUE)
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  
  ht <- ComplexHeatmap::pheatmap(as.matrix(FC_result[majorCelltypes, ]),
                                 name = "Cellular\nInfluence\nlog2FC",
                                 color = my_palette,
                                 breaks = breaks,
                                 cluster_rows = TRUE,
                                 cluster_cols = FALSE,
                                 show_colnames = TRUE,
                                 show_rownames = TRUE,
                                 scale = "none",
                                 border_color = "black",
                                 cellwidth = 10,
                                 cellheight = 10,
                                 treeheight_row = 20,
                                 treeheight_col = 20)
  
  # Generate heatmap and save to PDF
  pdf(paste0('./output/log2FC_spatwt_', celltype, '.pdf'), height = 5, width = 5)
  draw(ht)
  dev.off()
}

# Define the list of marker Expr objects 
new_expr_id<- as.data.frame(new_expr)
new_expr_id$sample_id <- output$sample_ids
expr_list<- list(new_expr_id[idx_tum,] , new_expr_id[idx_cd8T,])  
output_names <- c("tumor", "cd8T")
markers<- list(c("CD86","CK", "KI67", "PDL1", "VISTA"),
               c("ICOS", "GranzB", "PD1", "HLADR", 
                 "TIGIT", "LAG3", "PTPN22", "pSTAT3", "TIM3"))  

for (idx in seq_along(expr_list)) {
  expr_PA <- expr_list[[idx]]  
  colnames(expr_PA)<- gsub("pSTAT", "pSTAT3", colnames(expr_PA))
  colnames(expr_PA)<- gsub("Granzyme", "GranzB", colnames(expr_PA))
  
  expr_PA$Site <- factor(output$Site[match(expr_PA$sample_id, output$meta_data$sample_id)], levels = sitelevels)
  expr_PA$Patient <- factor(output$Patient[match(expr_PA$sample_id, output$meta_data$sample_id)], levels = unique(output$Patient))
  
  
  celltype <- output_names[idx]  
  majorMarkers<- markers[[idx]]
  
  # Initialize empty matrices
  expr_pan <- c()
  expr_liv <- c()
  
  # Calculate means for Pancreas
  for (i in matchedPA) {
    tmp <- apply(expr_PA[expr_PA$Site == "Pancreas" &
                           expr_PA$Patient == i, majorMarkers], 2, mean)
    expr_pan <- cbind(expr_pan, tmp)
  }
  colnames(expr_pan) <- paste0(matchedPA, "_P")
  
  # Calculate means for Liver
  for (i in matchedPA) {
    tmp <- apply(expr_PA[expr_PA$Site == "Liver" &
                           expr_PA$Patient == i, majorMarkers], 2, mean)
    expr_liv <- cbind(expr_liv, tmp)
  }
  colnames(expr_liv) <- paste0(matchedPA, "_L")
  
  # Compute log2 fold-change
  FC_result <- log2(expr_liv / expr_pan)
  FC_result[FC_result == 'Inf'] <- NA
  FC_result[FC_result == '-Inf'] <- NA
  FC_result[is.na(FC_result)] <- NA
  colnames(FC_result) <- matchedPA
  
  # Define palette and breaks
  my_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  max_abs <- max(abs(FC_result[majorMarkers, ]), na.rm = TRUE)
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  
  ht <- ComplexHeatmap::pheatmap(as.matrix(FC_result[majorMarkers, ]),
                                 name = "Marker\nExpression\nlog2FC",
                                 color = my_palette,
                                 breaks = breaks,
                                 cluster_rows = TRUE,
                                 cluster_cols = FALSE,
                                 show_colnames = TRUE,
                                 show_rownames = TRUE,
                                 scale = "none",
                                 border_color = "black",
                                 cellwidth = 10,
                                 cellheight = 10,
                                 treeheight_row = 20,
                                 treeheight_col = 20)
  
  # Generate heatmap and save to PDF
  pdf(paste0('./output/log2FC_expr_', celltype, '.pdf'), height = 5, width = 5)
  draw(ht)
  dev.off()
}



