pacman::p_load(dplyr, stringr, tidyr, dplyr, data.table)

source("~/Documents/Genizon/Genizon_Scripts/cleaner/get_loc_pairs.R")

omega_cross_coalescence <- function(loc_pairs_fn = "/Users/luke/Documents/Genizon/Data/loc_pairs_may2022.RDS", 
                                          wts_loc_fn = "/Users/luke/Documents/Genizon/Data/watershed_locations_feb2022.csv",
                                          pedigree_fn = "/Users/luke/Documents/Genizon/BALSAC/Balsac_aout_2021_v2/tout_balsac.csv",
                                          cross_coal_fn = "~/Documents/Genizon/cross_coalescence/cat_pairwise_bottleneck.csv",
                                          rerun = F) {
  if(rerun){
    bw <- 20  # bin width (km)
    dl <- 200 # dist limit (km)
    loc_pairs <- get_pairwise_locations(fread(wts_loc_fn), bin_width = bw, distance_limit = dl)
    saveRDS(loc_pairs, loc_pairs_fn)
  }
  
  N_probands <-
    fread(pedigree_fn) %>%
    filter(!ind %in% father,        # not fathers
           !ind %in% mother,        # or mothers
           !is.na(datem),           # who married
           datemp > 1900) %>%       # hard cap at 1900
    group_by(lieum) %>% 
    tally(name = "N") %>% 
    ungroup()
  
  # reference town
  ref_loc <- fread(wts_loc_fn) %>% 
    dplyr::select(lieu, a, t) %>% 
    dplyr::rename(lieum = lieu, 
                  ref_area = a, 
                  ref_town = t)
  
  ref_probands <-
    N_probands %>%  
    left_join(ref_loc, by = "lieum") %>%
    dplyr::rename(ref_town_size = N) %>%
    mutate(town_ref = lieum) %>%
    dplyr::select(town_ref, lieum, ref_town_size)
  
  # distal town
  distal_loc <- fread(wts_loc_fn) %>% 
    dplyr::select(lieu, a, t) %>% 
    dplyr::rename(lieump = lieu, 
                  distal_area = a, 
                  distal_town = t)
  
  distal_probands <-
    N_probands %>% 
    dplyr::rename(lieump = lieum) %>%
    left_join(distal_loc, by = "lieump") %>%
    dplyr::rename(distal_town_size = N) %>%
    mutate(town_distal = lieump) %>%
    dplyr::select(town_distal, lieump, distal_town_size)
  
  loc_pairs <- readRDS(loc_pairs_fn) %>%
    dplyr::rename(ref_town = t_2, ref_area = a_2,
                  distal_town = t_1, distal_area = a_1)
  
  one_sided_cross_coalescence <- 
    # this file contains ONE side of all pairwise comparisons
    fread(cross_coal_fn) %>%
    dplyr::rename(lieum = town_B, lieump = town_A) %>%
    rename_all(funs(
      stringr::str_replace_all(., '_A', '_distal') %>%
        stringr::str_replace_all(., '_B', '_ref')
    ) ) %>%
    left_join(ref_probands, by = "lieum") %>%
    left_join(distal_probands, by = "lieump") %>%
    left_join(ref_loc, by = "lieum") %>%
    left_join(distal_loc, by = "lieump") %>%
    left_join(loc_pairs, by = c("ref_town", "ref_area", "distal_town", "distal_area")) %>%
    dplyr::select(lieum, ref_town, ref_area, ref_town_size, kinship_ref, kinship_distal_ref,
           lieump, distal_town, distal_area, distal_town_size, kinship_distal, bin, w, d)
  
  pairwise_cross_coalescence <-
    # get pairwise by folding over and binding rows
    one_sided_cross_coalescence %>%
    dplyr::rename(
      lieum = lieump, 
      lieump = lieum,
      ref_town = distal_town,
      ref_area = distal_area,
      distal_town = ref_town,
      distal_area = distal_area,
      kinship_ref = kinship_distal,
      kinship_distal = kinship_ref,
      ref_town_size = distal_town_size,
      distal_town_size = ref_town_size
    ) %>%
    bind_rows(one_sided_cross_coalescence)
  
  omega_cross_coalescence <-  
    pairwise_cross_coalescence %>%
    mutate(distal_within_coalescence = 2 * kinship_distal / (distal_town_size * (distal_town_size - 1)),
           ref_within_coalescence = 2 * kinship_ref / (ref_town_size * (ref_town_size - 1)),
           cross_coalescence = kinship_distal_ref / (2 * ref_town_size * distal_town_size),
           percent_overlap = (2 * cross_coalescence) / (distal_within_coalescence + ref_within_coalescence) * 100,
           omega = percent_overlap
    ) %>%
    dplyr::select(omega, lieum, ref_town, ref_area, ref_town_size, lieump, distal_town, distal_area, distal_town_size, bin, w, d)
  
  return(omega_cross_coalescence)
}
