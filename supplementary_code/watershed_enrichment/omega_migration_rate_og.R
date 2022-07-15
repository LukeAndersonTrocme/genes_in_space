pacman::p_load(dplyr)
pacman::p_load(tidyr)

omega_migration_rate_og <- function(loc_pairs_fn = "/Users/luke/Documents/Genizon/Data/loc_pairs_may2022.RDS", 
                                 wts_loc_fn = "/Users/luke/Documents/Genizon/Data/watershed_locations_feb2022.csv",
                                 pedigree_fn = "/Users/luke/Documents/Genizon/BALSAC/Balsac_aout_2021_v2/tout_balsac.csv",
                                 rerun = F) {
  
  if(rerun){
    bw <- 20  # bin width (km)
    dl <- 200 # dist limit (km)
    loc_pairs <- get_pairwise_locations(fread(wts_loc_fn), bin_width = bw, distance_limit = dl)
    saveRDS(loc_pairs, loc_pairs_fn)
  }
  
  loc_pairs <- readRDS(loc_pairs_fn)
  
  # reference town
  ref_loc <- fread(wts_loc_fn) %>% 
    dplyr::select(lieu, a, t) %>% 
    dplyr::rename(lieum = lieu, 
                  a_2 = a, t_2 = t) %>% 
    mutate(ref_town = t_2,
           ref_area = a_2)
  
  # distal town
  distal_loc <- fread(wts_loc_fn) %>% 
    dplyr::select(lieu, a, t) %>% 
    dplyr::rename(lieump = lieu, 
                  a_1 = a, t_1 = t) %>% 
    mutate(distal_town = t_1,
           distal_area = a_1)
  
  total_pedigree <-  
    fread(pedigree_fn) %>%
    filter(!is.na(lieum))
  
  pedigree_migrations <-
    total_pedigree %>%
    left_join(distal_loc, by = "lieump") %>% 
    left_join(ref_loc, by = "lieum") %>% 
    group_by(distal_town, distal_area, 
             ref_town, ref_area,
             a_1, t_1, a_2, t_2) %>%  
    tally(name = "migrants") %>%
    ungroup()

  distal_town_size <- 
    total_pedigree %>% 
    left_join(distal_loc, by = "lieump") %>% 
    group_by(distal_town, distal_area) %>% 
    tally(name = "distal_town_size") %>%
    ungroup()
  
  omega_migration_rate <-                                     # compute out-bound migration rate
    left_join(pedigree_migrations,                            # number of migrants from `A` to `B`
              distal_town_size, by = c("distal_town","distal_area")) %>%       # total population size of town `A`
    mutate(emigration_rate = migrants / distal_town_size) %>% # m_A_B = migrants_A_B / pop_size_A
    left_join(loc_pairs, by = c("a_1", "t_1",
                                "a_2", "t_2")) %>%
    filter(bin != "") %>%
    mutate(omega = emigration_rate) %>%
    distinct(ref_town, ref_area, omega, distal_town_size,
             distal_town, distal_area, bin, w, d)
  
  return(omega_migration_rate)
}