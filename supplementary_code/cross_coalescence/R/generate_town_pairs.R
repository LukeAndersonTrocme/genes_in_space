# load libraries
library(dplyr)
library(geosphere) 
library(readr)

setwd("/home/luke1111/projects/ctb-sgravel/luke1111/cross_coalescence/R/")

# first we generate a list of towns to compute within-town coalesence
# define probands 
town_probands <- 
  data.table::fread("../data/tout_balsac.csv") %>% 
  filter(!ind %in% father,        # not fathers
         !ind %in% mother,        # or mothers
         !is.na(lieum),           # who married
         !is.na(datem),           # who married
         datemp > 1900)           # hard cap at 1900

# at least 100 people in a town
min_size <- town_probands %>% group_by(lieum) %>% tally() %>% filter(n>100) %>% arrange(-n)

write.table(min_size$lieum, file = "../list_of_towns.txt", col.names = FALSE, row.names = FALSE)

# load files
loc <- data.table::fread("../data/watershed_locations_feb2022.csv") %>% dplyr::rename(lieum = lieu)

src <- loc 
names(src) <- paste0(names(src),"A") 

snk <- loc 
names(snk) <- paste0(names(snk),"B")

loc_pairs <-
  tidyr::crossing(snk, src) %>%                              # crossing all combinations of towns
  mutate(d = geosphere::distHaversine(cbind(LonA, LatA),           
                                      cbind(LonB, LatB))/1e3,# distance between locations
         bin = cut(d, c(0, 1, seq(20, 200, 20)),             # clump into distance bins
                   right = F),
         w = ifelse(WSCSSDAA == WSCSSDAB,                    # determine if same watershed
                    "same", "diff"))  %>% 
  filter(!is.na(bin), d > 0) %>%
  dplyr::select(lieumA, aA, tA, lieumB, aB, tB, w, bin, d)

coalescence_town_pairs <- loc_pairs %>% 
  filter(lieumA %in% min_size$lieum, lieumB %in% min_size$lieum) %>% 
  dplyr::select(lieumA, lieumB) %>% 
  filter(lieumA < lieumB) %>% 
  arrange(lieumA, lieumB)

# write a file for each town with a list of all other town pairs
coalescence_town_pairs %>%
  group_by(lieumA) %>%
  group_walk(~ write_csv(.x, paste0("../town_pairs/", .y$lieumA, ".csv"), col_names = FALSE))