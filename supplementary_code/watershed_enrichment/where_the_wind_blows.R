pacman::p_load(dplyr)
pacman::p_load(tidyr)
pacman::p_load(ggplot2)
pacman::p_load(ggspatial)
pacman::p_load(sf)

source("~/Documents/Genizon/Genizon_Scripts/cleaner/get_loc_pairs.R")

where_the_wind_blows <- function(loc_pairs_fn = "/Users/luke/Documents/Genizon/Data/loc_pairs_may2022.RDS", 
                                 wts_loc_fn = "/Users/luke/Documents/Genizon/Data/watershed_locations_feb2022.csv",
                                 pedigree_fn = "/Users/luke/Documents/Genizon/BALSAC/Balsac_aout_2021_v2/tout_balsac.csv",
                                 stats_fn = "/Users/luke/Documents/Genizon/BALSAC/Balsac_aout_2021_v2/balsac_summary_stats.csv",
                                 political_processed_fn = "~/Documents/geographic_data/political_processed.RDS",
                                 wts_processed_fn = "~/Documents/geographic_data/wts_processed.RDS",
                                 simplified_water_area_fn = "~/Documents/geographic_data/simplified_water_area.RDS",
                                 elevation_fn = "/Users/luke/Documents/can_dem/quebec_raster_may2022_0.001_1e7_new_crs.RDS",
                                 crs_string = "+proj=omerc +lat_0=46.8560266 +lonc=-71.6218555 +alpha=0 +k_0=.7 +datum=WGS84 +units=m +no_defs +gamma=35",
                                 rerun = F) {
  
  if(rerun){
    bw <- 20  # bin width (km)
    dl <- 200 # dist limit (km)
    loc_pairs <- get_pairwise_locations(fread(wts_loc_fn), bin_width = bw, distance_limit = dl)
    saveRDS(loc_pairs, loc_pairs_fn)
  }
  
  loc_pairs <- readRDS(loc_pairs_fn)
  
  
  loc_m <- fread(wts_loc_fn) %>% 
    dplyr::select(lieu, a, t) %>% 
    dplyr::rename(lieum = lieu, 
                  a_2 = a, t_2 = t) %>% 
    mutate(name_2 = paste(t_2, a_2))
  
  loc_mp <- loc_m %>% 
    dplyr::rename(name_1 = name_2,
                  lieump = lieum, 
                  a_1 = a_2, t_1 = t_2)
  
  political_processed <- readRDS(political_processed_fn)
  wts_processed <- readRDS(wts_processed_fn)
  simplified_water_area <- readRDS(simplified_water_area_fn)
  m_gplot <- readRDS(file = elevation_fn) %>% filter(value > 0) %>% mutate(value = ifelse(value >= 1000, 999,value))
  
  lon_1 <- fread(wts_loc_fn) %>% 
    dplyr::rename(a_1 = a, t_1 = t, Lon_1 = Lon, Lat_1 = Lat) %>% distinct(a_1, t_1, Lon_1, Lat_1)
  lon_2 <- fread(wts_loc_fn) %>% 
    dplyr::rename(a_2 = a, t_2 = t, Lon_2 = Lon, Lat_2 = Lat) %>% distinct(a_2, t_2, Lon_2, Lat_2)
  
  space_time_contributions <- 
    fread(stats_fn) %>%
    left_join(fread(pedigree_fn), by = "ind") %>%
    left_join(loc_m) %>% left_join(loc_mp) %>%
    group_by(a_1, t_1, a_2, t_2) %>%
    dplyr::summarise(migrants = n(),
                     genomes = sum(expected_contribution, na.rm = T),
                     gdate = weighted.mean(datem, w = expected_contribution, na.rm = T)) %>%
    left_join(loc_pairs)
  
  rate_summary <- 
    space_time_contributions %>%
    left_join(lon_1) %>% 
    left_join(lon_2, by = c("a_2", "t_2")) %>%
    filter(d > 10, genomes > 10, a_2 != "Northern Quebec", !is.na(Lon_1), !is.na(Lon_2)) %>%
    group_by(a_2, t_2) %>% 
    dplyr::mutate(rank = rank(-genomes, ties.method = "first")) %>%
    filter(rank == 1) %>%
    ungroup() %>%
    mutate(discrete_size = case_when(genomes > 10000 ~ 10000,
                                     genomes > 1000 ~ 1000,
                                     TRUE ~ 100),
           genome_bins = cut(genomes, c(1,1e3,1e4,Inf), right = F))
  
  out <-
    ggplot()+
    geom_contour_filled(data = m_gplot, aes(x = x, y = y, z = value), bins = 12, alpha = 0.9, color = "NA") +
    scale_fill_grey(name = "Elevation (m)", start=0.9, end=0.2) + 
    geom_sf(data=political_processed %>% filter(ctry_en %in% c("Canada"), !juri_en %in% c("Quebec")), fill = "gray95", color = NA) +
    geom_sf(data=political_processed %>% filter(!ctry_en %in% c("Canada","Ocean")), fill = "gray95", color = NA)+
    geom_sf(data=wts_processed, fill = NA, color = "black", size = 0.15) +
    geom_sf(data = simplified_water_area, size = 0.15, color = "aliceblue", fill = "aliceblue") +
    ggspatial::geom_spatial_segment(data = rate_summary %>% filter(genome_bins == "[1,1e+03)"),
                                    aes(x= Lon_1, 
                                        y = Lat_1, 
                                        xend = Lon_2, 
                                        yend = Lat_2,
                                        color = gdate,
                                        size = genome_bins
                                    ),
                                    lineend = "butt", linejoin = "mitre",
                                    arrow = arrow(length = unit(0.075,"cm"), type = "closed"),
                                    alpha = 1
    ) +
    ggspatial::geom_spatial_segment(data = rate_summary %>% filter(genome_bins == "[1e+03,1e+04)"),
                                    aes(x= Lon_1, 
                                        y = Lat_1, 
                                        xend = Lon_2, 
                                        yend = Lat_2,
                                        color = gdate,
                                        size = genome_bins
                                    ),
                                    lineend = "butt", linejoin = "mitre",
                                    arrow = arrow(length = unit(0.075,"cm"), type = "closed"),
                                    alpha = 1
    ) +
    ggspatial::geom_spatial_segment(data = rate_summary %>% filter(genome_bins == "[1e+04,Inf)"),
                                    aes(x= Lon_1, 
                                        y = Lat_1, 
                                        xend = Lon_2, 
                                        yend = Lat_2,
                                        color = gdate,
                                        size = genome_bins
                                    ),
                                    lineend = "butt", linejoin = "mitre",
                                    arrow = arrow(length = unit(0.075,"cm"), type = "closed"),
                                    alpha = 1
    ) +  
    scale_size_manual( values = c(0.2, 0.8, 2), name = "Genomes Contributed", labels = c("< 1k", "< 10k", "> 10k")) +
    scale_color_viridis_c(name = "Time",
                          option = "B", end = 0.90,
                          breaks = c(1700, 1825, 1950),
                          limits = c(1650, 2000)) +
    guides(color = guide_colorbar(direction = "horizontal",title.position="left", barwidth = 8.5),
           size = guide_legend(direction = "horizontal",title.position="top", override.aes = aes(lineend = "butt", linejoin = "mitre")),
           fill = "none") +
    theme_bw() +
    coord_sf(crs = crs_string,
             xlim = c(-249000,382000),
             ylim = c(-118000,160000)) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks= element_blank(),
          legend.box = "horizontal",
          legend.background = element_rect(fill=alpha("white", 0.5)),
          legend.key = element_rect(fill = 'transparent'),
          legend.position = c(0.79, 0.93),
          legend.text = element_text(size=16),
          legend.title = element_text(size=20),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "aliceblue"))+  
    ggspatial::annotation_scale(location = "bl", width_hint = 0.2) +
    ggspatial::annotation_north_arrow(location = "bl", which_north = "true", 
                                      pad_x = unit(0.1, "in"), pad_y = unit(0.3, "in"),
                                      style = north_arrow_fancy_orienteering)
  return(out)
}