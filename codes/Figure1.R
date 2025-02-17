#============================#
###        Figure 1        ###
#============================#

colorassigned<- c("#8c564b",#Immune_mix,
                  "#1f77b4",#CD8T
                  "#ff7f0e",#CD4T
                  "#ff7f0e",#Treg
                  "#FEE480",#NK
                  rep("#d62728",6), 
                  "#9C4DF5",#Neutrophil 
                  rep("#BBE5E9",7), 
                  "#d3d3d3",#Tumor
                  "#5c6068")#UAs
names(colorassigned)<-clusterlevels

# Figure 1B - cell phenotype heatmap
plot_scaled_heatmap(expr = new_expr,
                    markers = c("CD3","CD8", "CD4","FOXP3","CD45RA","CD57", 
                                "CD86","CD68", "CD163", "CD206", "ARG1", "HLADR",
                                "DCSIGN", "CD15",
                                "Collagen","CD74","SMAVIM", "Podoplanin",
                                "CK"), 
                    colorassigned = colorassigned, 
                    heatmapcluster = F, 
                    cell_clustering1m = factor(output$cell_clustering2m, levels = clusterlevels), 
                    filename ="./output/Figure1B.pdf") 



# Figure 1C - cell mask
coord <- fsApply(output$fcs1, exprs)[, c("X_position", "Y_position")]
coord <-as.data.frame(coord)
coord$sample_id <- output$sample_ids
coord$cluster <- output$cell_clustering2m


visualize_ProbabilityMask(expr0 = coord, # should contain sampleID, X,Y position, celltype annotation 
                          samp_id = 51, 
                          colorassigned = colorassigned)

visualize_ProbabilityMask(expr0 = coord, # should contain sampleID, X,Y position, celltype annotation 
                          samp_id = 26, 
                          colorassigned = colorassigned)


# abundance plots 
counts_table <- table(output$cell_clustering2m, output$sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)
areas <- read_xlsx('./Config/new_areas.xlsx')
densities <- t(t(counts)/areas$TotalArea)
densities<- densities[clusterlevels, ]

# save data 
#write.csv(counts,'./output/Results_counts.csv')
#write.csv(props,'./output/Results_props.csv')
#write.csv(densities, './output/Results_densities.csv')

# celltype abundance by Site =====
counts_frac<- counts%>% rownames_to_column("cellTypes")
counts_frac <- melt(counts_frac)
colnames(counts_frac)[2]<-'sample_id'

counts_frac$Site<- factor(output$Site[match(counts_frac$sample_id,output$meta_data$sample_id)], levels=sitelevels)
counts_frac$Patient <- factor(output$Patient[match(counts_frac$sample_id,output$meta_data$sample_id)], levels=unique(output$Patient))
counts_frac$cellTypes <- factor(counts_frac$cellTypes, levels=clusterlevels)


# total number of cells for pancreas and liver 
total_counts<- counts_frac%>% 
  group_by(Site)%>% 
  summarise(count = sum(value))

pct_total_counts<- melt(total_counts)%>% 
  mutate(pct = value/sum(value)*100)
pct_total_counts$variable<-"total fraction"

# Site fraction for all cells
p0<- ggplot(pct_total_counts, aes(x = variable, y = pct, fill = Site)) +
  geom_col(position="fill", width = 0.7) +
  labs(title = "", x = "", y = "") +
  theme_minimal() + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.text.x = element_text(size=10,
                                   color = "black",
                                   angle=45, 
                                   hjust =1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "top", 
        legend.title = element_blank(),  
        legend.text = element_blank(),    
        legend.key = element_blank(), 
        legend.background = element_blank()) 

# group by cell type and site, then take average of values 
table1<- counts_frac %>% 
  group_by(cellTypes,Site) %>% 
  mutate(aggr_counts = sum(value)) 

table1<-table1[!duplicated(table1[c("cellTypes","Site")]),] 
dim(table1) # 46 X 7 (23 cellTypes * 2 Sites)

# calculate Site fractions for each cellTypes
table1<- table1%>% 
  group_by(cellTypes)%>% 
  mutate(pct = aggr_counts/sum(aggr_counts)*100)

# Reorder levels of the Celltype factor based on pct values (sort with Liver)
order <- (table1[table1$Site=="Liver",]$cellTypes)[order(table1[table1$Site=="Liver",]$pct)]

table1$cellTypes <- factor(table1$cellTypes, levels=order)

table1 <- table1 %>% 
  arrange(cellTypes, Site)

p1<- ggplot(table1, aes(x=cellTypes, y=pct, fill=Site))+
  geom_col(position="fill", width = 0.7) +
  theme_minimal() + 
  scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x * 100)) +
  labs(title = "", x = "", y = "")+ 
  theme(axis.text=element_text(size=10,color = "black"),
        axis.text.x = element_text(angle=45, hjust =1),
        legend.position="top") 

# by cores 
# by patients 
# by cellTypes

celltype_colors <- hcl.colors(length(clusterlevels), palette ="Spectral")
# celltype abundance for each core 

table2<- counts_frac%>% 
  group_by(sample_id)%>% 
  mutate(pct = value/sum(value)*100)

p2<-ggplot(table2, aes(x=sample_id, y=pct, fill=cellTypes))+
  geom_col(position="fill", width = 0.8) +
  scale_fill_manual(values = celltype_colors) + 
  theme_minimal() + 
  scale_y_continuous(expand = c(0, 0),labels = function(x) paste0(x * 100)) +
  labs(x = "", y = "")+ 
  theme(plot.title = element_text(size = 10,  hjust = 0.5), 
        axis.text=element_text(size=10, color = "black"),
        axis.text.x = element_text(angle=90, vjust =0.5, size = 6), 
        legend.position="right", 
        strip.background = element_blank(), 
        strip.text = element_text(size=10, color = "black"))+ 
  guides(fill = guide_legend(ncol = 2)) +  
  facet_wrap(~Site, ncol=2, scales = "free_x", drop = TRUE)

# celltype abundance for each Patient =====

table3<- counts_frac %>% 
  group_by(cellTypes, Patient, Site) %>% 
  mutate(aggr_counts = sum(value)) 

table3<-table3[!duplicated(table3[c("cellTypes","Patient", "Site")]),] 
dim(table3) 

# calculate Site fractions for each cellTypes
table3<- table3%>% 
  group_by(Patient)%>% 
  mutate(pct = aggr_counts/sum(aggr_counts)*100)


p3<- ggplot(table3, aes(x=Patient, y=pct, fill=cellTypes))+
  geom_col(position="fill", width = 0.8) +
  scale_fill_manual(values = celltype_colors) + 
  theme_minimal() + 
  scale_y_continuous(expand = c(0, 0),labels = function(x) paste0(x * 100)) +
  labs(x = "", y = "")+ 
  theme(plot.title = element_text(size = 10,  hjust = 0.5), 
        axis.text=element_text(size=10, color = "black"),
        axis.text.x = element_text(angle=45, hjust =1), 
        legend.position="right", 
        strip.background = element_blank(), 
        strip.text = element_text(size=10, color = "black"))+ 
  guides(fill = guide_legend(ncol = 2)) +  
  facet_wrap(~Site, ncol=2, scales = "free_x", drop = TRUE)

# patient proportion for each celltype =====
Patientcolors <- hcl.colors(length(unique(table3$Patient)), palette = "Spectral")
names(Patientcolors)<- unique(table3$Patient)
p4<- ggplot(table3, aes(x=cellTypes, y=pct, fill=Patient))+
  geom_col(position="fill", width = 0.8) +
  scale_fill_manual(values = Patientcolors) + 
  theme_minimal() + 
  scale_y_continuous(expand = c(0, 0),labels = function(x) paste0(x * 100)) +
  labs(x = "", y = "")+ 
  theme(plot.title = element_text(size = 10,  hjust = 0.5), 
        axis.text=element_text(size=10, color = "black"),
        axis.text.x = element_text(angle=45, hjust =1), 
        legend.position="right", 
        strip.background = element_blank(), 
        strip.text = element_text(size=10, color = "black"))+ 
  guides(fill = guide_legend(ncol = 2)) +  
  facet_wrap(~Site, ncol=2, scales = "free_x", drop = TRUE)

pdf("./output/Figure1D.pdf", width=7, height = 4)
print(plot_grid(p1, p0, rel_widths=c(7,1)))
dev.off()

pdf("./output/AbundancePlots.pdf", width=9, height = 3.5)
print(p2)
print(p3)
print(p4)
dev.off()


# density plot =====
# plot for selected cellTypes
selected <- c("M_I", "M_IV", "M_V", 
              "Str_II","Str_VI", "Tumor")
densities_2 <- densities[selected,]
pdf('./output/Figure1E.pdf',width=4.5, height=4)
density_cores<- density_byCore(densities = densities_2, 
                               specimen = Specimen_designation, 
                               coi = selected, 
                               site1_ref = "Pancreas", 
                               site2="Liver", 
                               ncols = 3)
print(density_cores)
dev.off()


# proportion of cell types (compare between pancreas and liver)
props_df<- melt(props%>% rownames_to_column("cellTypes"))
colnames(props_df)[2]<- "sample_id"
props_df<- left_join(props_df, Specimen_designation[,c("sample_id", "Site")], by="sample_id")

props_df <- props_df %>%
  group_by(cellTypes, Site) %>%
  mutate(median_value = median(value, na.rm = TRUE)) %>%
  ungroup()

props_df<- props_df[order(props_df$median_value, decreasing=T), ]
props_df$cellTypes<- factor(props_df$cellTypes, levels=rev(unique(props_df$cellTypes)))
props_df$Site<- factor(props_df$Site, levels=sitelevels)
pdf('./output/FigureS1G_1.pdf', height=5, width=5)
prop1<-ggplot(props_df, 
       aes(x=value, y=cellTypes, color=fct_rev(Site))) +
  geom_boxplot()+
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text=element_text(size=12, color = "black"),
        legend.position = 'right')+
  scale_x_continuous(expand = expansion(mult = c(0,0.3))) +
  xlim(0,75) +
  ylab('All cellTypes')+
  xlab("proportion of all cells")+
  scale_color_manual(name = "Site", 
                     values = scales::hue_pal()(length(unique(props_df$Site))),
                     breaks = levels(props_df$Site))
print(prop1)
dev.off()


# remove Tumor and Str -> compare the proportion of Lymp, and Mac cellTypes 
select_cellTypes <- c("CD8T","CD4T", "Treg","NK","Immune_Mix","Neutrophil",paste0("M_", c("I", "II", "III", "IV", "V", "VI")))
counts_table_rm <- counts_table[select_cellTypes, ]
props_table_rm <- t(t(counts_table_rm) / colSums(counts_table_rm)) * 100
props_rm <- as.data.frame.matrix(props_table_rm)

props_df_rm<- melt(props_rm%>% rownames_to_column("cellTypes"))
colnames(props_df_rm)[2]<- "sample_id"
props_df_rm<- left_join(props_df_rm, Specimen_designation[,c("sample_id", "Site")], by="sample_id")

props_df_rm <- props_df_rm %>%
  group_by(cellTypes, Site) %>%
  mutate(median_value = median(value, na.rm = TRUE)) %>%
  ungroup()

props_df_rm<- props_df_rm[order(props_df_rm$median_value, decreasing=T), ]
props_df_rm$cellTypes<- factor(props_df_rm$cellTypes, levels=rev(unique(props_df$cellTypes)))
props_df_rm$Site<- factor(props_df_rm$Site, levels=sitelevels)

pdf('./output/FigureS1G_2.pdf', height=3, width=5)
prop2<-ggplot(props_df_rm, 
       aes(x=value, y=cellTypes, color=fct_rev(Site))) +
  theme_bw() + 
  geom_boxplot()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text=element_text(size=12, color = "black"),
        legend.position = 'right')+
  scale_x_continuous(expand = expansion(mult = c(0,0.3))) +
  xlim(0,75) +
  ylab('Immune Cells')+
  xlab("proportion of Immune cells")+
  scale_color_manual(name = "Site", 
                     values = scales::hue_pal()(length(unique(props_df_rm$Site))),
                     breaks = levels(props_df_rm$Site))
print(prop2)
dev.off()
