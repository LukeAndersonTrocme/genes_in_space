pacman::p_load(dplyr)
pacman::p_load(tidyr)

sum_over_sets_of_towns <- function(input_omega_metric) {
  
  within_watershed <- 
    input_omega_metric %>% 
    filter(w == "same_t") %>%
    group_by(ref_town, ref_area, bin) %>%
    dplyr::summarise(Omega_within = sum(omega, na.rm = T),
                     c_b_within = sum(weight, na.rm = T),
                     n_within = n()) %>%
    ungroup()
  
  between_watershed <- 
    input_omega_metric %>% 
    filter(w == "diff_t") %>%
    group_by(ref_town, ref_area, bin) %>%
    dplyr::summarise(Omega_between = sum(omega, na.rm = T),
                     c_b_between = sum(weight, na.rm = T),
                     n_between = n()) %>%
    ungroup()
  
  winthin_and_between <- 
    full_join(within_watershed, between_watershed, 
              by = c("ref_town", "ref_area", "bin")) %>%
    filter(!is.na(n_between), !is.na(n_within)) %>%
    replace_na(list(Omega_between = 0, Omega_within = 0)) %>%
    dplyr::mutate(eta = Omega_within / (Omega_within + Omega_between),
                  c_b = c_b_within / (c_b_within + c_b_between),
                  e_b = eta / c_b
    ) %>%
    filter(!is.na(e_b))
  
  return(winthin_and_between)
}
