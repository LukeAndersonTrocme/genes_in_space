pacman::p_load(dplyr)
pacman::p_load(tidyr)
pacman::p_load(geosphere)

get_pairwise_locations <- function(locations, bin_width = 25, distance_limit = 200){

loc <- locations %>% dplyr::select(Lon, Lat, t, a, WSCSSDA) # required columns

src <- loc                                                # set all towns as source towns `A`
names(src) <- paste0(names(src),"_1")                      # add suffix for permutation

snk <- loc                                                # set all towns as source towns `B`
names(snk) <- paste0(names(snk),"_2")                      # add suffix for permutation

bins <- c(0, 1, seq(bin_width, distance_limit, bin_width))# set distance bins

loc_pairs <-
  tidyr::crossing(snk, src) %>%                           # crossing all combinations of towns
  mutate(d = geosphere::distHaversine(cbind(Lon_1, Lat_1),           
                           cbind(Lon_2, Lat_2))/1e3,        # measure distance between locations
         bin = cut(d, bins, right = F),                   # clump pairs of towns into distance bins
         w = ifelse(WSCSSDA_1 == WSCSSDA_2,                 # determine if towns share the same watershed
                    "same_t", "diff_t")) %>%  
  dplyr::select(a_1, t_1, a_2, t_2, w, bin, d)   # subset columns of interest

return(loc_pairs)
}