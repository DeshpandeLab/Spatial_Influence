## Density plot (by cores)
density_byCore <- function(densities, 
                           specimen, 
                           coi, 
                           site1_ref, 
                           site2, 
                           ncols){
  require(ggplot2); require(ggpubr);
  
  
  melted <- melt(densities)
  colnames(melted)[1] <- "celltypes"
  colnames(melted)[2]<-"sample_id"
  melted$sample_id<- gsub("X","", melted$sample_id) 
  melted$sample_id<- factor(melted$sample_id, levels = levels(specimen$sample_id))
  
  melted<- left_join(melted, specimen, by="sample_id")
  melted$Site <- factor(melted$Site, levels = c(site1_ref,site2))
  
  # take mean by cores
  table<- melted %>% 
    group_by(celltypes, sample_id, Site) %>% 
    mutate(avg= mean(value))
  table<-table[!duplicated(table[c("celltypes","sample_id","Site")]),]
  
  
  df = table[table$celltypes%in%coi,]
  df$celltypes<- factor(df$celltypes, levels= coi)
  

  p<- ggplot(df, aes(Site, avg, fill=Site)) + 
    geom_point(size = 2, shape=22, color="black", stroke=0.8, 
               position = position_jitter(width = 0.1, height = 0),
               alpha=0.8)+
    xlab("") +
    ylab(expression(No.~of~cells/mm^2))+ 
    theme_bw() + 
    theme(axis.text.x=element_text(size=12,color="black", angle=45, vjust=0.5, hjust=0.5), 
          axis.text.y=element_text(size=12, color="black"),
          axis.title=element_text(size=12, color="black"),
          legend.position="none", 
          strip.background = element_blank(), 
          strip.text = element_text(size=12, color = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_rect(color = "black")) +
    stat_summary(fun.y = median, geom = "errorbar", 
                 aes(ymax = ..y.., ymin = ..y.., group =factor(Site)),
                 width = 0.6)  +
    stat_summary(geom = "errorbar", fun.data = mean_sd, position = "dodge", width=0.3) + 
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    stat_compare_means(method="wilcox.test",
                       label="p.signif", 
                       hide.ns = F, 
                       paired=F, 
                       label.y.npc =0.9,
                       label.x.npc = 0.5, 
                       size=5) +
    facet_wrap(~celltypes, scales = "free_y", ncol = ncols)
  return(p)
}

