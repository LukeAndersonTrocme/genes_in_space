pacman::p_load(dplyr)
pacman::p_load(data.table)


get_ibd_in_quebec <- function(individual_ibd,
                              pedigree_file = "/Users/luke/Documents/Genizon/BALSAC/Balsac_aout_2021_v2/tout_balsac.csv",
                              cartagene_id_key = "/Users/luke/Documents/Genizon/BALSAC/ind_file111_IID.txt",
                              location_file = "/Users/luke/Documents/Genizon/Data/watershed_locations_feb2022.csv") {
  # Geographic locations
  geo_locs <- fread(location_file) %>% dplyr::rename(lieump = lieu)
  # ID conversion key for CAG
  cag_id <- fread(cartagene_id_key)
  # IDs from the Genizon cohort
  ibd_ids <- unlist(c(unique(individual_ibd$IID1), unique(individual_ibd$IID2)))

  # get all IDs and their locations
  iid_loc <-
    # load genealogy
    fread(pedigree_file) %>%
    # keep individuals with genetic data
    filter(ego_genizon %in% ibd_ids | ind %in% cag_id$ind) %>%
    # get conversion key
    left_join(. , cag_id) %>%
    # combine ID columns for IBD merge
    mutate(IID = case_when(!is.na(file111) ~ as.character(file111),
                           is.na(file111) ~ as.character(ego_genizon))) %>%
    inner_join(geo_locs, by = "lieump") %>%
    # keep columns of interest
    dplyr::select(IID, a, t) %>%

  return(iid_loc)
}
