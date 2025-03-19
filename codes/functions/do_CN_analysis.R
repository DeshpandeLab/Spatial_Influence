#--------------------------------#
#          CN analysis           #
#--------------------------------#


do_CN_analysis <- function(spatwt_df = spatwt_df_tum, #cell of interest filtered, should contain patient, sample_id
                           expr0 = expr0_tum, 
                           funcMarkers = funcMarkers, 
                           clusterby = clusterlevels, #should match the independent variables of spatwt
                           numClusters = 10, 
                           scaleOption = F, 
                           rmCN = F,  # default = NA 
                           filenameprefix = "TUM_")
  {
  require(reshape); require(ggplot2); require(pheatmap);
  require(dplyr); require(Hmisc);require(grid); require(RColorBrewer); 
  require(dittoSeq); require(tidyverse); require(ggpubr);
  
  # used clusterfunc using mean metric 
  somClustering <- clusterfcs_dim(fcs=as.matrix(spatwt_df[, clusterby]),
                              cluster_by = clusterby,
                              xdims=10, ydims=10, 
                              numclusters=numClusters, 
                              scaleoption = scaleOption)
  
  spatwt_df$CN <- somClustering$cell_clustering
  print(table(spatwt_df$CN, spatwt_df$sample_id))
  
  CNcolors<- dittoColors(reps = 1, get.names = FALSE)[1:numClusters]
  
  ### Dot plot for CN clusters 
  chTME<-c()
  # mean of neighborhood weights 
  for (i in 1:numClusters){
    tmp<- apply(spatwt_df[spatwt_df$CN==i,clusterby],2, median)
    chTME<- cbind(chTME, tmp)
    colnames(chTME)[i]<-paste0('CN', i)
  }
  scaled_chTME <- t(apply(chTME, 1, scale))
  colnames(scaled_chTME)<- colnames(chTME)
  
  df1<- reshape2::melt(scaled_chTME)
  names(df1)<-c("celltype",'tme', "zscore")
  df2<- reshape2::melt(chTME)
  names(df2)<-c("celltype",'tme', "median")
  
  df_joined <- left_join(df1, df2, by = c("celltype", "tme"))
  
  df_joined$tme <- factor(df_joined$tme, levels = rev(levels(df_joined$tme)))

  # plot only Median weights > 0.01 
  df_filtered <- df_joined %>% 
    dplyr::filter(!is.na(zscore))%>% 
    dplyr::filter(median>0.01) # greater than 1% influence 
  
  ## dot plot for CN clustering ===
  p2<-ggplot(df_filtered, aes(x=celltype, y=tme)) + 
    geom_point(pch = 21, stroke = 0.5, col="black", aes(size = median, fill = zscore)) + 
    theme_classic() +
    labs(size = "Median\nCellular Influences", fill="z-score")+ 
    ylab("") + 
    xlab("")+
    theme(axis.line = element_blank(), 
          axis.text.x= element_text(colour ="black", size = 9),
          axis.text.y= element_text(colour ="black", size = 12),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.7),  
          panel.grid.major = element_line(size = 0.5, linetype = 'solid'),  
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid'), 
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7), 
          legend.key.size = unit(0.2, "cm"),
          legend.position="right") + 
    scale_color_discrete("black")+
    scale_fill_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0)+
    scale_y_discrete(position = "left") + 
    scale_x_discrete(
      guide = guide_axis(angle = 45),position = "top") + 
    scale_size_continuous(breaks = seq(0, max(df_joined$median, na.rm=TRUE), by = 0.1))
  
  legend <- get_legend(p2)
  
  
  pdf(paste0('./output/', filenameprefix, 'CNdotplot.pdf'), width=5, height=3)
  print(p2+theme(legend.position = "none"))
  dev.off()
  
  
  
  ## stacked bar of CN abundance ===
  df<- as.data.frame(table(spatwt_df$CN, spatwt_df$Site))
  names(df)<- c("TME", "Site", "Freq")
  df$TME <- paste0("CN", df$TME)
  df$TME <- factor(df$TME, levels=paste0("CN", 1:numClusters))
  df$prop <- with(df, Freq / ave(Freq, Site, FUN = sum))
  
  
  print("CN frequency")
  print(df)
  
  pdf(paste0('./output/', filenameprefix, 'CNstackedbar.pdf'), width=2, height=4)
  p3<-ggplot(df, aes(x = Site, y = prop, fill = TME)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(expand = c(0, 0),labels = function(x) paste0(x * 100)) +
    labs(title = "", x = "", y = "") +
    theme_minimal() +
    scale_fill_manual(values = CNcolors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text= element_text(colour ="black", size = 12),
          legend.key.size = unit(0.4, "cm"))
  print(p3)
  dev.off()
  
  ## Patient composition of CN 
  df0<- as.data.frame(table(spatwt_df$CN, paste0(spatwt_df$Patient, "_",spatwt_df$Site)))
  names(df0)<- c("TME", "Patient", "Freq")
  df0$TME <- paste0("CN", df0$TME)
  df0$TME <- factor(df0$TME, levels=rev(paste0("CN", 1:numClusters)))
  df0$prop <- with(df0, Freq / ave(Freq, TME, FUN = sum))
  
  patientColors <- hcl.colors(length(unique(df0$Patient)), palette = "roma")
  p4<-ggplot(df0, aes(x = prop, y = TME, fill = Patient)) +
    geom_col(width = 0.8)  +
    labs(title = "", x = "", y = "", fill="Samples") +
    theme_bw() +
    scale_fill_manual(values = patientColors) +
    theme(axis.text.x  = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.6),  
          legend.key.size = unit(0.4, "cm"))+ 
    scale_x_continuous(expand = expansion(mult = c(0,0)))+ 
    scale_y_discrete(expand = expansion(mult = c(0,0)))
  legend_2<- get_legend(p4)
  
  pdf(paste0('./output/', filenameprefix, 'CN_composition.pdf'), width=2, height=3)
  print(p4+theme(legend.position = "none"))
  dev.off()

  pdf(paste0('./output/', filenameprefix, 'CNdotplot_legend.pdf'), width=5, height=4)
  plot(legend)
  plot(legend_2)
  dev.off()
  
  if (rmCN) {
    rm<- df[df$Freq == max(df$Freq), ]$TME
    idxrm <- as.numeric(gsub("CN","", rm))
    pdf(paste0('./output/', filenameprefix, 'CNstackedbar_selected.pdf'), width=2, height=4)
    p5<- ggplot(df[df$TME%nin%rm, ], aes(x = Site, y = prop, fill = TME)) +
      geom_bar(stat = "identity") +
      scale_y_continuous(expand = c(0, 0),labels = function(x) paste0(x * 100)) +
      labs(title = "", x = "", y = "") +
      theme_minimal() +
      scale_fill_manual(values = CNcolors[-idxrm]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text= element_text(colour ="black", size = 12),
            legend.key.size = unit(0.4, "cm"))
    print(p5)
    dev.off()
    
    df1<- as.data.frame(table(spatwt_df$CN,spatwt_df$Site, spatwt_df$Patient))
    names(df1)<- c("TME","Site" ,"Patient", "Freq")
    df1$TME <- paste0("CN", df1$TME)
    df1$TME <- factor(df1$TME, levels=paste0("CN", 1:numClusters))
    df1$prop <- with(df1, Freq / ave(Freq, Patient, Site, FUN = sum))
    
    pdf(paste0('./output/', filenameprefix, 'CNstackedbar_patients.pdf'), width=4, height=4)
    p6<- ggplot(df1[df1$TME%nin%rm, ], aes(x = Patient, y = prop, fill = TME)) +
      geom_bar(stat = "identity") +
      scale_y_continuous(expand = c(0, 0),labels = function(x) paste0(x * 100)) +
      labs(title = "", x = "", y = "") +
      theme_minimal() +
      scale_fill_manual(values = CNcolors[-idxrm]) +
      theme(axis.text.x = element_text(angle =90, vjust = 0.5, hjust=1, size = 8),
            axis.text= element_text(colour ="black", size = 11),
            legend.key.size = unit(0.4, "cm"))+ 
      facet_wrap(~ Site)
    print(p6)
    dev.off()
    
  }
    
  
  ## line plot for functional markers ===
  expr0$TME <- somClustering$cell_clustering
  
  expr0$TME<- paste0("CN",expr0$TME)
  expr0$TME<- factor(expr0$TME, levels=paste0("CN",1:numClusters))
  
  ## Overall functional state comparison between Sites 
  data <- melt(expr0)
  data_filtered <- data[data$variable%in%funcMarkers, ]

  data_for_plot<- data_filtered

  summary_df <- data_for_plot %>%
    group_by(Site, variable)%>% 
    dplyr::mutate(baselevel = median(value))%>%
    group_by(TME, Site, variable, baselevel) %>%
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
  
  summary_df$TME<- factor(summary_df$TME, levels=paste0("CN", c(1:numClusters)))
  
  pdf(paste0('./output/', filenameprefix, 'median_expr_level_1_3_quantile.pdf'), height=6, width=8)
  p7<-ggplot(summary_df, aes(x = TME, y = median_value, color = Site, group = Site)) +
    geom_line(size = 1) +  # Line for the mean value
    geom_ribbon(aes(ymin = first_quantile, ymax = third_quantile, fill = Site), alpha = 0.2, color = NA) +
    facet_wrap(~ variable, scales = "free_y") +  # Facet by Marker
    theme_minimal() + 
    labs(
      x = "TME",
      y = "Scaled Expression",
      color = "Site",
      fill = "Site"
    ) +
    theme(
      legend.position = "right",
      strip.text = element_text(size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
      
    ) + 
    geom_hline(data= summary_df, 
               aes(yintercept=baselevel, color=Site), 
               linetype="dashed")
  print(p7)
  dev.off()
  
  pdf(paste0('./output/', filenameprefix, 'median_expr_level_95conf.pdf'), height=6, width=8)
  p8<-ggplot(summary_df, aes(x = TME, y = median_value, color = Site, group = Site)) +
    geom_line(size = 1) +  # Line for the mean value
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Site), alpha = 0.2, color = NA) +
    facet_wrap(~ variable, scales = "free_y") +  # Facet by Marker
    theme_minimal() + 
    labs(
      x = "TME",
      y = "Scaled Expression",
      color = "Site",
      fill = "Site"
    ) +
    theme(
      legend.position = "right",
      strip.text = element_text(size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
      
    )+ 
    geom_hline(data= summary_df, 
               aes(yintercept=baselevel, color=Site), 
               linetype="dashed")
  print(p8)
  dev.off()
  
  return(somClustering$cell_clustering)
  
}
