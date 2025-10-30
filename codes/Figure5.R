#============================#
###       Figure 5         ###
#============================#

# Tumor TME analysis ====
spatwt<- readRDS('./data/spatwt.rds')

# filter all tumors (remove core 43)
index <- which(df_output$cluster=="Tumor"&
                 df_output$sample_id!=43)

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

visualize_ProbabilityMask(expr0 = df_output, 
                          samp_id = 18, 
                          colorassigned = colorassigned)
visualize_ProbabilityMask(expr0 = df_output, 
                          samp_id = 25, 
                          colorassigned = colorassigned)

# visualize tumor classification
dd <- df_output[index,] 
dd$type<- "BK"
dd[subset, "type"]<- "BD"
dd$type<- as.factor(dd$type)

pdf('./output/Figure5B.pdf', height=6, width=3)
p5B<-ggplot(dd[dd$sample_id%in%c("18", "25"),], aes(x=X_position, y=Y_position, color=type))+
  geom_point(size = 0.5)+ 
  scale_color_manual(values = c("BD" = "black", "BK" = "#DF536C")) + 
  theme_void() + 
  scale_y_reverse()+
  theme(legend.position = "none",)+
  facet_wrap(~sample_id, ncol=1)
print(p5B)
dev.off()

# proportion of CNs 
expr_tum<- as.data.frame(new_expr[index,])
expr_tum$sample_id <- factor(df_output$sample_id[index], levels=samplevels)
expr_tum$Patient <-factor(df_output$Patient[match(expr_tum$sample_id,  df_output$sample_id)], levels=unique(df_output$Patient))
expr_tum$Site <- factor(df_output$Site[match(expr_tum$sample_id,  df_output$sample_id)], levels=sitelevels)

tum_markers<- c("KI67", "CK", "CD86", "VISTA", "PDL1")
markerExpr <- expr_tum[, c(tum_markers, "Site", "sample_id")]
markerExpr$type <- "BK"
markerExpr[subset, "type"]<- "BD"


CNprop <- as.data.frame(table(markerExpr$type, markerExpr$Site))
names(CNprop)<-c("CN", "Site", "Freq")
CNprop<- CNprop%>% 
  dplyr::group_by(Site)%>% 
  dplyr::mutate(total_counts = sum(Freq), 
         prop=Freq/total_counts*100)

CNprop_1 <- as.data.frame(CNprop)
CNprop_1$CN <- as.character(CNprop_1$CN)
CNprop_1$CN[CNprop_1$CN!="BK"]<-"BD"

CNprop_1 <- CNprop_1 %>% 
  dplyr::group_by(Site, CN) %>% 
  dplyr::summarise(counts = sum(Freq), .groups = "drop") %>%  # Calculate counts per Site and CN
  dplyr::group_by(Site) %>% 
  dplyr::mutate(
    total_counts = sum(counts),    # Total counts per Site
    prop = (counts / total_counts) * 100  # Proportion in percentage
  ) %>% 
  ungroup()


CNprop_1$Site<- factor(CNprop_1$Site, levels=c("Liver", "Pancreas"))
pdf('./output/Figure5C.pdf', height=2, width=4.5)
p5C<- ggplot(CNprop_1, aes(x=Site, y=prop, fill=CN))+
  scale_fill_manual(values = c("BK" = "#DF536C", "BD" = "black"))+
  geom_col(position="fill", width = 0.7) +
  theme_minimal()+
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x * 100)) +
  labs(title = "", x = "", y = "", fill="Tumor\nClassification")+ 
  theme(axis.text=element_text(size=14,color = "black"),
        legend.position="right", 
        legend.text = element_text(size=13))+ 
  coord_flip()
print(p5C)
dev.off()



#====

data_for_violin <- reshape::melt(cbind(spatwt_df_tum[, clusterlevels], markerExpr))
data_for_violin <- data_for_violin %>%
  dplyr::mutate(type = ifelse(type == "BK", "BK", "BD"))

summary_df_2 <- data_for_violin %>%
  dplyr::group_by(variable, sample_id,type, Site)%>% 
  dplyr::summarise(
    mean_value = mean(value))

selected_celltypes <- c("M_I", "Treg", "Str_I", "Str_II")
subset_df <- subset(summary_df_2, variable%in%selected_celltypes)
subset_df$variable<- factor(subset_df$variable, levels=selected_celltypes)
# Cellular influences on Tumors (Compare between Sites)
pdf('./output/Figure5D.pdf', height=5, width=3)
p5D<- ggplot(subset_df, aes(x = type, y = mean_value, fill = Site)) +
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
  ylab("Mean Cellular Influences") +
  xlab("Tumor type")+ 
  theme(legend.position = "top", 
        axis.text = element_text(size=12, color="black"), 
        axis.title=element_text(size=12, color="black"))
print(p5D)
dev.off()


selected_markers <- c("CK","CD86")

# Tumor Marker expression level between Sites
summary_df_2$Tumor_type <- paste0(summary_df_2$Site, "_", summary_df_2$type)
summary_df_2$Tumor_type<- factor(summary_df_2$Tumor_type, levels = c("Pancreas_BK", "Pancreas_BD", "Liver_BK", "Liver_BD"))
comparisons <- list(c("Pancreas_BK", "Pancreas_BD"), c("Pancreas_BK", "Liver_BK"), c("Pancreas_BK","Liver_BD"), 
                    c("Pancreas_BD", "Liver_BD"), c("Liver_BK", "Liver_BD"))

pdf('./output/Figure5E.pdf', height=3.5, width=3)
p5E<- ggplot(summary_df_2[summary_df_2$variable%in%selected_markers, ], aes(x = Tumor_type, y = mean_value, fill = Tumor_type)) +
  geom_violin(alpha = 0.6, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 0.8), alpha = 0.5) +
  stat_compare_means(
    method = "wilcox.test", 
    comparisons = comparisons, 
    paired = FALSE, 
    label = "p.signif", 
    hide.ns = FALSE, 
    size=3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~variable, scales = "free", ncol = 2) +  
  theme_bw() + 
  ylab("Mean Marker Expression") +
  xlab("")+ 
  theme(legend.position = "none", 
        axis.text = element_text(size=12, color="black"), 
        axis.title=element_text(size=12, color="black"), 
        axis.text.x = element_text(angle=45, hjust=1))
print(p5E)
dev.off()


site_comparison<- list(c("Pancreas", "Liver"))
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

data_for_heatmap <- reshape::melt(cbind(spatwt_df_tum[, clusterlevels], markerExpr))
data_for_heatmap <- data_for_heatmap %>%
  dplyr::mutate(type = ifelse(type == "BK", "BK", "BD"))


selected_celltypes <- c("Immune_Mix","CD8T","CD4T","Treg", "NK", 
                        "M_I","M_II","M_III","Neutrophil", "Str_I", "Str_II", 
                        "Str_III", "Str_IV", "Tumor", "UA")

selected_markers<- c("CK", "KI67", "PDL1", "CD86", "VISTA")

summary_df <- data_for_heatmap %>%
  dplyr::group_by(variable, sample_id,type, Site)%>% 
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

    pdf(paste0('./output/Figure5G_',data_type, '_tumor.pdf'), height=6, width=7)
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
    ComplexHeatmap::draw(ht)
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

unique_id<- unique(df_output$sample_id)


sampInfo <- data.frame(sample_id = unique_id)
sampInfo$Site<- as.factor(df_output$Site[match(sampInfo$sample_id,df_output$sample_id)])
sampInfo$Patient <- as.factor(df_output$Patient[match(sampInfo$sample_id,df_output$sample_id)])

plot_multiple_probability_masks(coord_data = df_output,
                                uniqID =select_sampID,
                                sampINFO = sampInfo,
                                probability_mask_folder = 'Probability_masks/',
                                colorassigned = colorassigned,
                                ncols = 4, 
                                filenameprefix= "Figure5F_")

plot_5F_data <- df_output[df_output$sample_id%in%select_sampID, ]
plot_5F_data$cluster <- factor(plot_5F_data$cluster,
                               levels = c("Tumor","M_I","M_II","Str_I","Str_II","Neutrophil", "Treg"))

cols <- c(
  Tumor   = "#d3d3d3",
  M_I = "#d62728",
  M_II = "#d62728",
  Str_I   = "#BBE5E9",
  Str_II = "#1f77b4",
  Neutrophil = "#9C4DF5",
  Treg   = "#2ca02c", 
  others ="white"
)

p5F<- ggplot(plot_5F_data, aes(x=X_position, y=Y_position, color=cluster))+
  geom_point(size=1)+
  facet_wrap(~sample_id, nrow=1) +
  scale_y_reverse()+
  scale_color_manual(values = cols, drop = FALSE, na.value = "white") + 
  coord_equal()+
  theme_minimal()+
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())

pdf('./output/Fig5F.pdf', width=20, height=8)
print(p5F)
dev.off()

index <- which(df_output$cluster=="Tumor"&
                 df_output$sample_id%in%select_sampID)

plot_marker_level(coord_data = df_output[index, ],
                  expression_data= new_expr[index, ], 
                  uniqID = select_sampID,
                  marker ="CK", 
                  filenameprefix="Tumor_CK",
                  height_size = 5, 
                  width_size = 18, 
                  ncols=4)

# Supplementary Figure 4
df <- cbind(spatwt_df_tum[, clusterlevels], markerExpr)

selected_celltypes <- c("Immune_Mix","CD8T","CD4T","Treg", "NK", 
                        "M_I","M_II","M_III","Neutrophil", "Str_I", "Str_II", 
                        "Str_III", "Str_IV", "Tumor", "UA")

selected_markers<- c("CK", "KI67", "PDL1", "CD86", "VISTA")

celltype_summary <- df %>%
  dplyr::select(all_of(selected_celltypes), sample_id, Site) %>%
  dplyr::group_by(sample_id, Site) %>%
  dplyr::summarise(across(all_of(selected_celltypes), mean, na.rm = TRUE), .groups = "drop")

marker_summary <- df %>%
  dplyr::select(all_of(selected_markers), sample_id, Site) %>%
  dplyr::group_by(sample_id,Site) %>%
  dplyr::summarise(across(all_of(selected_markers), mean, na.rm = TRUE), .groups = "drop")

merged_summary <- left_join(celltype_summary, marker_summary, by = c("sample_id", "Site"))


CD8T_med <- median(merged_summary$CD8T)
PDL1_med <- median(merged_summary$PDL1)


merged_summary$type <- ifelse(
  merged_summary$CD8T >= CD8T_med & merged_summary$PDL1 >= PDL1_med, "CD8T+PDL1+",
  ifelse(merged_summary$CD8T <= CD8T_med & merged_summary$PDL1 >= PDL1_med, "CD8T-PDL1+",
         ifelse(merged_summary$CD8T >= CD8T_med & merged_summary$PDL1 <= PDL1_med, "CD8T+PDL1-", "CD8T-PDL1-")
  )
)


merged_summary$type<- factor(merged_summary$type, 
                             levels = c("CD8T+PDL1+", 
                                        "CD8T-PDL1-", 
                                        "CD8T-PDL1+",
                                        "CD8T+PDL1-"))

# PDL1 expr for tumor types (individual tumor cells)
tumor_types<- merged_summary[, c("sample_id", "type")]

pdf('./output/FigureS5A.pdf', height=3, width=3.5)
pS5A<- ggplot(merged_summary,
              aes(x=type, y=PDL1, color=type)) +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  geom_jitter(width=0.2, alpha=0.6, size=1)+ 
  theme_bw()+ 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15),
        axis.text.x = element_text(angle=45, color="black", hjust=1),
        plot.title = element_text(hjust = 0.5))+
  facet_wrap(~Site) + 
  xlab("type")
print(pS5A)
dev.off()

ddf <- as.data.frame(table(merged_summary$Site, merged_summary$type))

pdf('./output/FigS5B.pdf', width=2.5, height=4)
pS5B<-ggplot(ddf, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col(position="fill") +
  scale_y_continuous(expand = c(0, 0),labels = function(x) paste0(x * 100)) +
  labs(title = "", x = "", y = "", fill="type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text= element_text(colour ="black", size = 12),
        legend.key.size = unit(0.4, "cm"))
print(pS5B)
dev.off()

pdf('./output/FigureS5C.pdf', height=4, width=4)
pS5C<- ggplot(merged_summary,
       aes(x=type, y=VISTA, color=type)) +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  geom_jitter(width=0.2, alpha=0.6, size=1)+ 
  theme_bw()+ 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15),
        axis.text.x = element_text(angle=45, color="black", hjust=1),
        plot.title = element_text(hjust = 0.5))+
  facet_wrap(~Site)
print(pS5C)
dev.off()


ggplot(merged_summary,
       aes(x=VISTA, y=Neutrophil, color=type)) +
  geom_point(alpha=0.5) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))+
  stat_cor(aes(color=type), method='spearman', 
           label.y.npc="top", 
           label.x.npc = "left",
           label.sep='\n',
           size=4) + 
  geom_smooth(method=lm, aes(fill=type), alpha=0.2)+ 
  facet_wrap(~Site)


pdf('./output/FigureS5D.pdf', height=4, width=5)
pS5D<- ggplot(merged_summary,
       aes(x= VISTA, y=Neutrophil, color=type)) +
  geom_point(alpha=0.5) +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15),
        axis.text= element_text(color="black"),
        plot.title = element_text(hjust = 0.5))+
  #stat_cor(aes(color=type), method='spearman', 
  #         label.y.npc="top", 
  #         label.x.npc = "left",
  #         size=3) + 
  geom_smooth(method=lm, aes(fill=type), alpha=0.2)+ 
  facet_wrap(~Site)
print(pS5D)
dev.off()



## Matched Patients ==== 
# matched patients 
matchedPA<- intersect(unique(df_output$Patient[df_output$Site=="Pancreas"]), 
                      unique(df_output$Patient[df_output$Site=="Liver"]))

matched_sampID<- unique(df_output$meta_data$sample_id[df_output$Patient%in%matchedPA&
                                                        df_output$meta_data$sample_id!=43])

## Figure S5F-G
# observe cellular influence and marker expr level change by Patient (average by patient)
df<- left_join(df, Specimen_designation[, c("sample_id", "Patient")], by="sample_id")

df_filtered<- df%>% 
  dplyr::filter(Patient%in%c("PA2", "PA5", "PA6", "PA8", "PA10", "PA11"))


celltype_summary_2 <- df_filtered %>%
  dplyr::select(all_of(selected_celltypes), Patient, Site) %>%
  dplyr::group_by(Patient, Site) %>%
  dplyr::summarise(across(all_of(selected_celltypes), mean, na.rm = TRUE), .groups = "drop")

marker_summary_2 <- df_filtered %>%
  dplyr::select(all_of(selected_markers), Patient, Site) %>%
  dplyr::group_by(Patient, Site) %>%
  dplyr::summarise(across(all_of(selected_markers), mean, na.rm = TRUE), .groups = "drop")

merged_summary2 <- left_join(celltype_summary_2, marker_summary_2, by = c("Patient",  "Site"))

merged_summary2$Patient<- factor(merged_summary2$Patient, 
                                 levels=c("PA2", "PA5", "PA6", "PA8", "PA10", "PA11"))


merged_summary2_melt <- reshape::melt(as.data.frame(merged_summary2))

merged_summary2_celltype <- merged_summary2_melt%>% 
  dplyr:: filter(variable%in%c("CD8T", "CD4T", "Treg", "NK", 
                               "M_I", "M_II", "M_III", "Neutrophil", 
                               "Str_I", "Str_II"))
merged_summary2_melt$variable <- factor(merged_summary2_melt$variable)

pdf('./output/FigureS5F.pdf', height=4, width=8)
pS5F<- ggplot(merged_summary2_celltype, aes(x=Site, y=value, group=Patient))+
  geom_line(aes(color=Patient), size = 1) + 
  geom_point(aes(shape=Patient))+ 
  facet_wrap(~variable, scales = "free_y", ncol = 5)+
  theme_bw() +
  labs(
    title = "Cellular influence comparison at patient level",
    x = "",
    y = "Average influence"
  ) +
  scale_x_discrete(labels = c("PRI", "MET")) + 
  theme(axis.text.x = element_text(size=14, angle=45,  hjust =1))
#guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
print(pS5F)
dev.off()

merged_summary2_marker <- merged_summary2_melt%>% 
  dplyr:: filter(variable%in%c("CK", "KI67","PDL1", "CD86", "VISTA"))

pdf('./output/FigureS5G.pdf', height=4, width=4.5)
pS5G<- ggplot(merged_summary2_marker, aes(x=Site, y=value, group=Patient))+
  geom_line(aes(color=Patient), size = 1) + 
  geom_point(aes(shape=Patient))+ 
  facet_wrap(~variable, scales = "free_y", ncol = 3)+
  theme_bw() +
  labs(
    title = "Marker Expression comparison at patient level",
    x = "",
    y = "Average marker expression"
  ) +
  scale_x_discrete(labels = c("PRI", "MET")) + 
  theme(axis.text.x = element_text(size=14, angle=45,  hjust =1))
#guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
print(pS5G)
dev.off()
