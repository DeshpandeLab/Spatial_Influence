#######################################
#             plot umap               #
#######################################


plotUmap <- function(umapRes,seed=1234,neighbors=10, 
                     midpoint,
                     color_clusters=colorassigned,
                     subtype_markers=NULL)
{require(umap);require(ggplot2);require(viridis);require(ggrepel); require(tidyverse)
  custom.settings = umap.defaults
  custom.settings$seed = seed
  custom.settings$n.neighbors = neighbors
  
  label.df <- umapRes %>% 
    dplyr::group_by(cell_clustering) %>% 
    dplyr::summarize(umap1 = mean(umap1), umap2 = mean(umap2))
  colnames(label.df)[1]<-'label'
  
  ggp1 <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = cell_clustering)) +
    geom_point(size = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    scale_color_manual(values = color_clusters, name="CLUSTERS") +
    ggrepel::geom_label_repel(data = label.df, aes(x = umap1, y = umap2, label = label, color = label),
                              show.legend = FALSE)+
    guides(color = guide_legend(override.aes = list(shape = 16, size = 3), ncol = 1))
    
    
  print(ggp1)
 
  ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = cell_clustering)) +
    geom_point(size = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    scale_color_manual(values = color_clusters, name="CLUSTERS") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
    
  print(ggp + facet_wrap(~ Site, ncol=4)+ggtitle('Site'))
  
  ggp3 <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = Site)) +
    geom_point(size = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 2))
  
  print(ggp3)
  
  ggp4 <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = Patient)) +
    geom_point(size = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
  
  print(ggp4)
  
  ggp2 <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = sample_id)) +
    geom_point(size = 1) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
  
  print(ggp2)
  

  
  #can specify which markers to display
  if(!is.null(subtype_markers)){
    for(i in subtype_markers)
    {
      ggp <- ggplot(umapRes,  aes(x = umap1, y = umap2, color = umapRes[,i])) +
        geom_point(size = 1) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        scale_color_gradient2(i, low="dark blue",mid="white",high="dark red", midpoint = mean(unlist(umapRes[,i])))
      print(ggp)
    }
  }
}
