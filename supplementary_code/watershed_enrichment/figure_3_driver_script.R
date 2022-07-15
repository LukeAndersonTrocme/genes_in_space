#  __  __       _                 
# |  \/  |     | |                
# | \  / | __ _| | _____          
# | |\/| |/ _` | |/ / _ \         
# | |  | | (_| |   <  __/         
# |_|__|_|\__,_|_|\_\___|         
# |  ____(_)                      
# | |__   _  __ _ _   _ _ __ ___  
# |  __| | |/ _` | | | | '__/ _ \ 
# | |    | | (_| | |_| | | |  __/ 
# |_|____|_|\__, |\__,_|_|  \___| 
# |__   __| |__/ |                
#    | |  | |___/_ __ ___  ___    
#    | |  | '_ \| '__/ _ \/ _ \   
#    | |  | | | | | |  __/  __/   
#    |_|  |_| |_|_|  \___|\___|   
                                 
pacman::p_load(dplyr, data.table,forcats, ggplot2, cowplot)
setwd("~/Documents/genes_in_space/supplementary_code/watershed_enrichment/")

source("sum_over_sets_of_towns.R")

# compute omega for three metrics
source("omega_migration_rate_og.R")
omega_migration_rate <- omega_migration_rate_og()
source("omega_identity_by_descent.R")
omega_identity_by_descent <- omega_identity_by_descent()
source("omega_cross_coalescence.R")
omega_cross_coalescence <- omega_cross_coalescence()

mig_enrch <- 
  omega_migration_rate %>% mutate(weight = 1, omega = omega * 1) %>%
  sum_over_sets_of_towns() %>% mutate(metric = "1_migration_rate")
  
ibd_enrch <-
  omega_identity_by_descent %>% mutate(weight = 1, omega = omega * 1) %>%
  sum_over_sets_of_towns() %>% mutate(metric = "2_ibd_sharing")

coal_enrch <- 
  omega_cross_coalescence %>% mutate(weight = 1, omega = omega * 1) %>%
  sum_over_sets_of_towns() %>% mutate(metric = "3_coal_sharing")

watershed_enrichment <-
  bind_rows(mig_enrch, ibd_enrch, coal_enrch) %>%   
  filter(ref_area != "St. Lawrence", 
         #n_between > 1, 
         #n_within > 1
         ) %>%
  group_by(bin, metric) %>%
  summarise_at(c("e_b"),
               funs(count = n(), sd, mean,
                    se = sd(.) / sqrt( n() ) 
               )) %>%
  filter(count > 10)

source("make_enrichment_plot.R")
enrichment_plot <- 
  make_enrichment_plot(watershed_enrichment)
  
source("where_the_wind_blows.R")
field_plot <- where_the_wind_blows()

fig3<-cowplot::plot_grid(field_plot, enrichment_plot, nrow = 1, rel_widths = c(1,0.25), labels = c("A","B"),  label_size = 30)#, align = "v", axis = "b" )
ggsave(fig3, file = "~/Documents/Genizon/Genizon_Scripts/Latex/Figures/fig3.8.jpg", 
       height = 6.5, width = 17.75, dpi = 400, bg = "white")


#   _____                   _                           _                    
#  / ____|                 | |                         | |                  
# | (___  _   _ _ __  _ __ | | ___ _ __ ___   ___ _ __ | |_ __ _ _ __ _   _ 
#  \___ \| | | | '_ \| '_ \| |/ _ \ '_ ` _ \ / _ \ '_ \| __/ _` | '__| | | |
#  ____) | |_| | |_) | |_) | |  __/ | | | | |  __/ | | | || (_| | |  | |_| |
# |_____/ \__,_| .__/| .__/|_|\___|_| |_| |_|\___|_| |_|\__\__,_|_|   \__, |
# |  ____(_)   | |   | |                                               __/ |
# | |__   _  __|_|_  |_|_ __ ___  ___                                 |___/ 
# |  __| | |/ _` | | | | '__/ _ \/ __|                                      
# | |    | | (_| | |_| | | |  __/\__ \                                      
# |_|    |_|\__, |\__,_|_|  \___||___/                                      
#           __/ |                                                          
#          |___/                                                           

administrative_regions <- 
  fread("~/Documents/Genizon/Data/qc_administrative_regions_lieu.csv") %>%
  dplyr::rename(combo_name = ref_town) %>%
  dplyr::select(combo_name, region, yax, rgb) %>%
  mutate(region = forcats::fct_reorder(region, yax))

watershed_enrichment_admin <-
  bind_rows(mig_enrch, ibd_enrch, coal_enrch) %>%   
  mutate(combo_name = paste(ref_town, ref_area)) %>%
  #filter(ref_area != "St. Lawrence") %>%
  inner_join(administrative_regions) %>%
  group_by(bin, metric, region) %>%
  summarise_at(c("e_b"),
               funs(count = n(), sd, mean,
                    se = sd(.) / sqrt( n() ) 
               )) %>%
  filter(count > 5)

admin <- make_supplementary_enrichment_plot(watershed_enrichment_admin)
ggsave(admin, height = 18, width = 16, filename = "~/Documents/Genizon/Genizon_Scripts/Latex/Figures/supplementary_enrichment_administrative.jpg")

watershed_enrichment_hydro <-
  bind_rows(mig_enrch, ibd_enrch, coal_enrch) %>%   
  mutate(region = case_when(
    ref_area == "Ottawa River" ~ "Outaouais",
    ref_area == "Nothern Qebec" ~ "Abitibi\nT\'emiscamingue",
    TRUE ~ ref_area
  )) %>%
  group_by(bin, metric, region) %>%
  summarise_at(c("e_b"),
               funs(count = n(), sd, mean,
                    se = sd(.) / sqrt( n() ) 
               )) %>%
  filter(count > 5) 

hydro <- make_supplementary_enrichment_plot(watershed_enrichment_hydro)
ggsave(hydro, height = 18, width = 12, filename = "~/Documents/Genizon/Genizon_Scripts/Latex/Figures/supplementary_enrichment_wts_a.jpg")

## Decay over space
all_dist <-
  bind_rows(omega_migration_rate %>% mutate(metric = "Migration Rate (log scale)"), 
            omega_identity_by_descent %>% mutate(metric = "Identity By Descent (cM) (log scale)"), 
            omega_cross_coalescence %>% mutate(metric = "Relative cross-coalescence (%)")) %>%
  mutate(combo_name = paste(ref_town, ref_area),
         watershed = case_when(w == "same_t" ~ "Same\nwatershed",
                               TRUE ~ "Different\nwatershed")) %>%
  inner_join(administrative_regions) %>%
  distinct(omega,region,rgb,d,watershed, ref_town, distal_town, metric)
saveRDS(all_dist, file = "~/Documents/Genizon/Data/all_dist.RDS")


