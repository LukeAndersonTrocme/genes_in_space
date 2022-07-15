pacman::p_load(dplyr)
pacman::p_load(ggplot2)

make_enrichment_plot <- function(watershed_enrichment){
  
  supp.labs <- c("Migration", "Identity by Descent", "Cross Coalescence")
  names(supp.labs) <- c("1_migration_rate", "2_ibd_sharing", "3_coal_sharing")
  
  ggplot(watershed_enrichment) +
    geom_hline(yintercept = 1, size = 0.5, color = alpha("black",0.1)) +
    geom_errorbar(aes(x = bin, y = mean, 
                      ymin = mean - se, 
                      ymax = mean + se), 
                  width = 0, size = 0.3) +
    geom_point(aes(x = bin, y = mean), size = 1.5) +
    labs(y = "Enrichment within watershed", 
         x = "Distance (Km)") + 
    theme_bw() + 
    scale_color_grey() + 
    scale_x_discrete(labels = c(20,40,60,80,100,120,140,160,180,200)) +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          strip.background = element_blank(),
          strip.text = element_text(size = 20),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank()) +
    facet_wrap(metric~., ncol = 1, scales = "free_y", 
               labeller = labeller(metric = supp.labs))
  
}


make_supplementary_enrichment_plot <- function(watershed_enrichment){
  
  supp.labs <- c("Migration", "Identity by Descent", "Cross Coalescence")
  names(supp.labs) <- c("1_migration_rate", "2_ibd_sharing", "3_coal_sharing")
  
  ggplot(watershed_enrichment) +
    geom_hline(yintercept = 1, size = 0.5, color = "grey60") +
    geom_errorbar(aes(x = bin, y = mean, 
                      ymin = mean - se, 
                      ymax = mean + se), 
                  width = 0, size = 0.3) +
    geom_point(aes(x = bin, y = mean),size = 1.5) +
    labs(y = "Enrichment within watershed", x = "Distance (Km)") + 
    theme_bw() + 
    scale_color_grey() + 
    scale_x_discrete(labels = c(20,40,60,80,100,120,140,160,180,200)) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          strip.background = element_blank(),
          strip.text = element_text(size = 20),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank()) +
    facet_grid(region~metric, labeller = labeller(metric = supp.labs), scales = "free_y")
  
}