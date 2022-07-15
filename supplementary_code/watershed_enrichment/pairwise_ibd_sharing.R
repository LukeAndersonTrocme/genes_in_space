pacman::p_load(dplyr)
pacman::p_load(tidyr)

pairwise_ibd_sharing <- function(individual_ibd, individuals_in_quebec, loc_pairs){
  
  left_side <- setNames(individuals_in_quebec, paste0(names(individuals_in_quebec), '_1')) %>% dplyr::rename(IID1 = IID_1)
  right_side <- setNames(individuals_in_quebec, paste0(names(individuals_in_quebec), '_2')) %>% dplyr::rename(IID2 = IID_2)
  
  pairwise_individuals <-
    # cross all combinations of pairs of individual
    tidyr::crossing(left_side, right_side) %>%
    left_join(loc_pairs, by = c("a_1","t_1","a_2","t_2")) %>%
    filter(!is.na(bin), IID1 != IID2)
  
  # get all pairwise combinations of IDs in Quebec
  pairwise_ibd_individuals <-
    pairwise_individuals %>%
    left_join(individual_ibd, by = c("IID1","IID2")) %>%
    mutate(name_1 = paste(t_1, a_1),
           name_2 = paste(t_2, a_2))
  
  # phase pairwise
  # because some towns with have individuals on 'left' or 'right' sides
  # we have to phase the IBD, take the sum and then 'unfold' the ibd
  # this will leave us with a full pairwise comparison between towns
  
  a <-
    pairwise_ibd_individuals %>%
    filter(name_1 < name_2) %>%
    dplyr::select(IID1, IID2, small, medium, large, name_1, name_2)
  
  b <-
    pairwise_ibd_individuals %>%
    dplyr::rename(name_1 = name_2,
                  name_2 = name_1,
                  IID1 = IID2, 
                  IID2 = IID1) %>%
    filter(name_1 < name_2) %>%
    dplyr::select(IID1, IID2, small, medium, large, name_1, name_2)
  
  flat_ibd_loc <-
    full_join(a, b, by=c('IID1', 'IID2', 'name_2', 'name_1'),
              suffix=c('1', '2')) %>%
    # when no IBD segments, replace with zero
    replace_na(list(small1 = 0, medium1 = 0 , large1 = 0,
                    small2 = 0, medium2 = 0 , large2 = 0)) %>%
    mutate(s = small1 + small2,
           m = medium1 + medium2,
           l = large1 + large2, 
           ibd = s + m + l
    ) %>%
    dplyr::select(IID1, IID2, name_1, name_2, ibd)
  
  # because we've phased and 'folded' the pairwise combinations
  # we now need to add the inverse direction for completeness
  
  complete_ibd <-
    flat_ibd_loc %>%
    filter(IID1 != IID2) %>%
    dplyr::rename(IID1 = IID2,
                  IID2 = IID1,
                  name_1 = name_2,
                  name_2 = name_1) %>%
    bind_rows(flat_ibd_loc)
  
  return(complete_ibd)
}
