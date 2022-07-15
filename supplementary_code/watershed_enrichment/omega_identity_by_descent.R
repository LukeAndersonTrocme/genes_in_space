pacman::p_load(dplyr)
pacman::p_load(tidyr)
source("~/Documents/Genizon/Genizon_Scripts/cleaner/get_ibd_in_quebec.R")
source("~/Documents/Genizon/Genizon_Scripts/cleaner/pairwise_ibd_sharing.R")

omega_identity_by_descent <- function(loc_pairs_fn = "/Users/luke/Documents/Genizon/Data/loc_pairs_may2022.RDS",
                                      wts_loc_fn = "/Users/luke/Documents/Genizon/Data/watershed_locations_feb2022.csv",
                                      pedigree_fn = "/Users/luke/Documents/Genizon/BALSAC/Balsac_aout_2021_v2/tout_balsac.csv",
                                      raw_ibd_fn = "/Users/luke/Documents/Genizon/Data/RDS/hap-ibd_2022-02-08_wide_sizeIID.RDS",
                                      temp_ibd_fn = "/Users/luke/Documents/Genizon/Data/individual_pairwise_ibd_june2022.RDS",
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
    mutate(name_2 = paste(t, a),
           a_2 = a, t_2 = t) %>% 
    dplyr::rename(lieum = lieu, 
                  ref_area = a, 
                  ref_town = t)
  
  # distal town
  distal_loc <- fread(wts_loc_fn) %>% 
    dplyr::select(lieu, a, t) %>% 
    mutate(name_1 = paste(t, a),
           a_1 = a, t_1 = t) %>% 
    dplyr::rename(lieump = lieu, 
                  distal_area = a, 
                  distal_town = t)
  
  if(rerun){ # 8 minute total run time
    individual_ibd <- readRDS(raw_ibd_fn)
    individuals_in_quebec <- get_ibd_in_quebec(individual_ibd)
    individual_pairwise_ibd <- pairwise_ibd_sharing(individual_ibd, 
                                                    individuals_in_quebec, 
                                                    loc_pairs)
    saveRDS(individual_pairwise_ibd, temp_ibd_fn)
  }
  # IID2, name_1, IID1, name_2, ibd
  individual_pairwise_ibd <- readRDS(temp_ibd_fn) 
  
  ibd_w <- 
    individual_pairwise_ibd %>%
    distinct(IID2, name_1) %>%
    group_by(name_1) %>% 
    tally(name = "distal_town_size") %>% 
    ungroup() %>%
    left_join(distal_loc, by = "name_1") %>%
    dplyr::select(distal_area, distal_town, distal_town_size)
  
  omega_identity_by_descent <-
    individual_pairwise_ibd %>%
    group_by(name_1, name_2) %>%
    dplyr::summarise(mean_ibd = mean(ibd)) %>%
    ungroup() %>%
    left_join(ref_loc, by = "name_2") %>% 
    left_join(distal_loc, by = "name_1") %>% 
    left_join(loc_pairs, by = c("a_1","t_1","a_2","t_2")) %>%
    left_join(ibd_w, by = c("distal_area","distal_town")) %>%
    mutate(omega = mean_ibd) %>%
    distinct(ref_town, ref_area, omega, distal_town_size,
             distal_town, distal_area, bin, w, d)
  
  return(omega_identity_by_descent)
}