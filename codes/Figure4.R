#============================#
###        Figure 4        ###
#============================#

## CD8+T TME analysis====
spatwt<- readRDS('./data/spatwt.rds')

# remove core43 
idx <- which(df_output$cluster=="CD8T"&
               df_output$sample_id!=43)

spatwt_df_cd8T_2m<- spatwt[idx, ]
spatwt_df_cd8T_2m$Site<- factor(df_output$Site[match(spatwt_df_cd8T_2m$sample_id,df_output$sample_id)], levels=sitelevels)
spatwt_df_cd8T_2m$sample_id<- factor(spatwt_df_cd8T_2m$sample_id, levels=unique(spatwt_df_cd8T_2m$sample_id))
spatwt_df_cd8T_2m$Patient <- factor(df_output$Patient[match(spatwt_df_cd8T_2m$sample_id,df_output$sample_id)], levels=unique(df_output$Patient))

expr_cd8T<- as.data.frame(new_expr[idx,])
expr_cd8T$sample_id<- spatwt_df_cd8T_2m$sample_id
expr_cd8T$Patient<- spatwt_df_cd8T_2m$Patient
expr_cd8T$Site <- spatwt_df_cd8T_2m$Site



colnames(expr_cd8T)<- gsub("pSTAT", "pSTAT3", colnames(expr_cd8T))
colnames(expr_cd8T)<- gsub("Granzyme", "GZMB", colnames(expr_cd8T))


funcMarkers<- c("GZMB", "ICOS", "LAG3", "PD1", "TIGIT", "TIM3", "HLADR", "PTPN22", "pSTAT3")

CD8T_cn_2m<-do_CN_analysis(spatwt_df = spatwt_df_cd8T_2m,
                           expr0 = expr_cd8T,
                           funcMarkers = funcMarkers,
                           clusterby = clusterlevels,
                           numClusters =8, 
                           scaleOption = F,
                           rmCN = FALSE,  # default = FALSE 
                           filenameprefix = "Figure4_")


# CN cell mask 
mask_data <- data.frame(sample_id = df_output$sample_id, cluster = df_output$cluster)
mask_data$CN<-"others"
mask_data <- mask_data[mask_data$sample_id!=43, ] 
mask_data[mask_data$cluster=="CD8T", "CN"]<- CD8T_cn_2m

# Figure 3C
# core 51 (pancreas)
plot_CN_cellmask2(sample_id = 51, 
                  coord= mask_data,
                  numClusters = 8,
                  path="./Probability_masks/")
# core 26 (Liver)
plot_CN_cellmask2(sample_id = 26, 
                  coord= mask_data,
                  numClusters = 8,
                  path="./Probability_masks/")

# plot zoomed CNs - validation of CN
colorassigned<- c("#1f77b4",#CD8T
                  "#ff6f0e",#CD4T
                  "#ff9f0e",#Treg
                  "#2ca02c",#M_I
                  "#d62728",#M_II
                  "#BBE5E9",#Str_I
                  "#FEE480",#Str_II
                  "#9C4DF5",#Neutrophil 
                  "#d3d3d3",#Tumor
                  rep("#5c6068",12))#others


majorCelltypes<- c("CD8T", "CD4T", "Treg", "M_I", "M_II", 
                   "Str_I", "Str_II", "Neutrophil", "Tumor")

names(colorassigned)<-c(majorCelltypes,
                        clusterlevels[!clusterlevels %in% majorCelltypes])

# flip y-axis to match cell mask & MCD image
coord <- df_output
coord$CN<-"others"
coord <- coord[coord$sample_id!=43, ] 
coord[coord$cluster=="CD8T", "CN"]<- CD8T_cn_2m

# Figure 3D (range +/-50)
for(i in 1:8){
  plot_zoomed_CN(sample_id = 51,
                 coord= coord,
                 colorassigned = colorassigned,
                 show = 1,
                 seed=1234,
                 range=50,
                 CNnum = i,
                 ncols=1)
}

for(i in 1:8){
  plot_zoomed_CN(sample_id = 26,
                 coord= coord,
                 colorassigned = colorassigned,
                 show = 1,
                 seed=1234,
                 range=50,
                 CNnum = i,
                 ncols=1)
}


# abundance plot of CD8T CNs for each core
cd8T_table <- as.data.frame(table(spatwt_df_cd8T_2m$sample_id, CD8T_cn_2m))
names(cd8T_table)<-c("sample_id", "CN", "Freq")

cd8T_table$Site<-factor(df_output$Site[match(cd8T_table$sample_id, df_output$sample_id)], levels=sitelevels)

CN_colors <- dittoColors()[c(1:length(unique(cd8T_table$CN)))]

names(CN_colors)<-c(1:length(unique(cd8T_table$CN)))

#CN_abundance_by_Cores 
ggplot(cd8T_table, aes(x=sample_id, y=Freq, fill=CN))+
  geom_col(position="fill", width = 0.8) +
  scale_fill_manual(values = CN_colors) + 
  theme_minimal()  + 
  scale_y_continuous(expand = c(0, 0),labels = function(x) paste0(x * 100)) +
  labs(x = "", y = "")+ 
  theme(plot.title = element_text(size = 10,  hjust = 0.5), 
        axis.text=element_text(size=10, color = "black"),
        axis.text.x = element_text(angle=90, vjust =0.5, size = 6), 
        legend.position="right", 
        strip.background = element_blank(), 
        strip.text = element_text(size=10, color = "black"))+ 
  guides(fill = guide_legend(ncol = 1)) +  
  facet_wrap(~Site, ncol=2, scales = "free_x", drop = TRUE)

matched_patients<-c("2","3","6","8","13","14","16","17","18","25","39","46","47","53","54","56","59")
# CN_abundance_by_Cores for matched patients only 
ggplot(cd8T_table[cd8T_table$sample_id%in%matched_patients,], aes(x=sample_id, y=Freq, fill=CN))+
  geom_col(position="fill", width = 0.8) +
  scale_fill_manual(values = CN_colors) + 
  theme_minimal()  + 
  scale_y_continuous(expand = c(0, 0),labels = function(x) paste0(x * 100)) +
  labs(x = "", y = "")+ 
  theme(plot.title = element_text(size = 10,  hjust = 0.5), 
        axis.text=element_text(size=10, color = "black"),
        axis.text.x = element_text(angle=90, vjust =0.5, size = 6), 
        legend.position="right", 
        strip.background = element_blank(), 
        strip.text = element_text(size=10, color = "black"))+ 
  guides(fill = guide_legend(ncol = 1)) +  
  facet_wrap(~Site, ncol=2, scales = "free_x", drop = TRUE)



# CN functional characteristic (CN relative)

markerExpr_cd8T <- expr_cd8T[, c(funcMarkers, "Site")]
markerExpr_cd8T$CN <- as.factor(paste0("CN",CD8T_cn_2m))
data_for_plot <- melt(markerExpr_cd8T)

summary_df <- data_for_plot %>%
  group_by(variable, CN, Site)%>% 
  dplyr::summarise(
    median_value = median(value))%>%
  group_by(variable)%>% 
  dplyr::mutate(scaled_median = as.numeric(scale(median_value)))

# scaled Marker expression (across CN)
p<- ggplot(summary_df, aes(variable, scaled_median, fill = Site)) +
  geom_bar(width = 1, stat = "identity", color = "black",  alpha = 0.7)  + 
  scale_y_continuous(breaks = scales::breaks_width(1))+
  geom_hline(yintercept = 0, color = "blue", size = 0.7) +
  #scale_fill_manual(values = tme_colors) +  # Custom colors
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  coord_polar() +
  facet_wrap(~CN, ncol=3)

pdf('./output/Figure4E.pdf', width=6, height = 8)
print(p)
dev.off()


# Supplementary Figure 3
## Tumor TME analysis====
# filter all tumors (remove core 43)
index <- which(df_output$cluster=="Tumor"&
                 df_output$sample_id!=43)

spatwt_df_tum <- spatwt[index, ]

# subset bulk tumor & tumor boundary 

cutoff = 0.7
subset <- which(spatwt_df_tum$Tumor<cutoff)


# do CN analysis and get CN classification
data4s <- spatwt_df_tum[subset, ]
data4s$Site<- factor(df_output$Site[match(data4s$sample_id,df_output$sample_id)], levels=sitelevels)
data4s$Patient <- factor(df_output$Patient[match(data4s$sample_id,df_output$sample_id)], levels=unique(df_output$Patient))



expr_tum<- as.data.frame(new_expr[index,])
expr_tum$sample_id <- factor(df_output$sample_id[index], levels=samplevels)
expr_tum$Site <- factor(df_output$Site[match(expr_tum$sample_id,df_output$sample_id)], levels=sitelevels)
expr_tum$Patient <- factor(df_output$Patient[match(expr_tum$sample_id,df_output$sample_id)], levels=unique(df_output$Patient))

expr_4s<- expr_tum[subset, ]


tum_markers<- c("CK", "CD86", "VISTA", "KI67", "PDL1")

numClust=12
tumor_cn<- do_CN_analysis(spatwt_df = data4s,
                          expr0 = expr_4s,
                          funcMarkers = tum_markers,
                          clusterby = clusterlevels,
                          numClusters = numClust, 
                          scaleOption = F, 
                          rmCN = TRUE, 
                          filenameprefix = "FigureS3_")


#saveRDS(markerExpr, './backup/markerExpr.rds')
markerExpr_tum <- expr_tum[, c(tum_markers, "Site", "sample_id")]
markerExpr_tum$type <- "BK"
markerExpr_tum[subset, "type"]<- tumor_cn

CNprop <- as.data.frame(table(markerExpr_tum$type, markerExpr_tum$sample_id))
names(CNprop)<-c("CN", "sample_id", "Freq")
CNprop$Patient <- factor(df_output$Patient[match(CNprop$sample_id,df_output$sample_id)], levels=unique(df_output$Patient))
CNprop$Site <- factor(df_output$Site[match(CNprop$sample_id,df_output$sample_id)], levels=sitelevels)

# Compare Tumor CN composition between Sites
CNprop_3 <- as.data.frame(CNprop)%>% 
  dplyr::filter(CN!="BK")
CNprop_3$CN <- factor(paste0("CN", CNprop_3$CN), levels=c(paste0("CN", 1:numClust)))

CNprop_3 <- CNprop_3 %>% 
  group_by(Site, CN, Patient, sample_id) %>% 
  summarise(counts = sum(Freq))%>% 
  ungroup()

CNprop_3 <- CNprop_3 %>% 
  group_by(sample_id) %>% 
  mutate(total_counts = sum(counts))%>% 
  ungroup()

CNprop_3$prop<- CNprop_3$counts/CNprop_3$total_counts*100
CNprop_3$Site<- factor(CNprop_3$Site, levels=sitelevels)

site_comparison<-list(c("Pancreas", "Liver"))

pdf('./output/FigureS3C.pdf', height=6, width =5)
pS3C<- ggplot(CNprop_3, aes(x = Site, y = prop, fill = Site)) +
  geom_violin(alpha=0.6)+
  geom_boxplot(width=0.3, fill = NA, color = "black", alpha = 0.5)+
  facet_wrap(~CN, scales = "free_y") + 
  stat_compare_means(
    method = "wilcox.test",
    comparisons = site_comparison,
    hide.ns = FALSE,
    label="p.format")+ 
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+ 
  theme_bw()+ 
  xlab("") + ylab("Proportion of CN")+
  theme(strip.text = element_text(size=12, color="black"), 
        axis.text.x = element_text(angle=45,  hjust =1), 
        axis.text = element_text(size=12, color="black"), 
        legend.position = "none")
print(pS3C)
dev.off()


# CN functional characteristic ====
data_for_plot <- melt(markerExpr_tum)
data_for_plot<- data_for_plot %>%
  dplyr::mutate(type = factor(
    ifelse(type == "BK", "BK", paste0("CN", type)), 
    levels = c("BK", paste0("CN", c(1:numClust)))
  ))


summary_df <- data_for_plot %>%
  group_by(variable, type, Site)%>% 
  dplyr::summarise(
    mean_value = mean(value), 
    median_value = median(value))%>%
  group_by(variable)%>% 
  dplyr::mutate(scaled_median = as.numeric(scale(median_value)))

# scaled Marker expression (across CN)
p<- ggplot(summary_df, aes(variable, scaled_median, fill = Site)) +
  geom_bar(width = 1, stat = "identity", color = "black",  alpha = 0.7)  + 
  scale_y_continuous(breaks = scales::breaks_width(1))+
  geom_hline(yintercept = 0, color = "blue", size = 0.7) +
  #scale_fill_manual(values = tme_colors) +  # Custom colors
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  coord_polar() +
  facet_wrap(~type, ncol=5)

pdf('./output/FigureS3D.pdf', width=8, height = 8)
print(p)
dev.off()

## line plot including bulk ====
data_melted<- melt(markerExpr_tum)
data_melted$type<- factor(data_melted$type, levels=c("BK", 1:numClust))

data_for_plot_line<- data_melted

summary_df_line <- data_for_plot_line %>%
  group_by(Site, variable)%>% 
  dplyr::mutate(baselevel = median(value))%>%
  group_by(type, Site, variable, baselevel) %>%
  dplyr::summarise(
    median_value = median(value),
    mean_value = mean(value),
    se = sd(value) / sqrt(n()),  # Standard error
    margin_of_error = 2 * se,  # Margin of error for 95% CI (z=1.96)
    lower_ci = median_value - margin_of_error,  # Lower bound of 95% CI
    upper_ci = median_value + margin_of_error,   # Upper bound of 95% CI
    first_quantile = quantile(value, 0.25), 
    third_quantile =quantile(value, 0.75) 
  )

# order manual
summary_df_line$type <- as.character(summary_df_line$type)

summary_df_line$type <- ifelse(summary_df_line$type != "BK",
                               paste0("CN", summary_df_line$type),
                               summary_df_line$type)


summary_df_line$type<- factor(summary_df_line$type, 
                              levels=c("BK", paste0("CN", c(2,1,3,4,7,8, 9, 10, 5, 12,6, 11))))

summary_df_line$variable<- factor(summary_df_line$variable, levels=tum_markers)

pdf('./output/FigureS3E.pdf', height=7, width=6)
pS3E<- ggplot(summary_df_line, aes(x = type, y = median_value, color = Site, group = Site)) +
  geom_line(size = 1) +  # Line for the mean value
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Site), alpha = 0.2, color = NA) +
  facet_wrap(~ variable, scales = "free_y", ncol=2) +  # Facet by Marker
  theme_minimal() + 
  labs(
    x = "TME",
    y = "Scaled Expression",
    color = "Site",
    fill = "Site"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle=45,  hjust =1),
    axis.title = element_text(size = 12)
  )+ 
  geom_hline(data= summary_df_line, 
             aes(yintercept=baselevel, color=Site), 
             linetype="dashed")
print(pS3E)
dev.off()


# Visualize UA cells 
colorassigned<- c("#FEE480",#Immune_mix
                  "#FEE480",#CD8T
                  "#FEE480",#CD4T
                  "#FEE480",#Treg
                  "#FEE480",#NK
                  rep("#d62728",6), 
                  "#d62728",#Neutrophil 
                  rep("#BBE5E9",7), 
                  "#d3d3d3",#Tumor
                  "#5c6068")#UAs
names(colorassigned)<-clusterlevels

visualize_ProbabilityMask(expr0 = df_output, 
                          samp_id = 49, 
                          colorassigned = colorassigned)
visualize_ProbabilityMask(expr0 = df_output, 
                          samp_id = 48, 
                          colorassigned = colorassigned)
