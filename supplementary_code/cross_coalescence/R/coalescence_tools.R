# Get the deepest generation for each internal node
maximum_genealogical_depth <- function(pedigree, list_of_probands) {
  # probands are individuals at generation 0
  individuals <- data.frame(ind = list_of_probands, generation = 0)
  # initialize output to append to downstream
  gen_depth <- individuals
  # for loop to iterate through generations
  for (k in 1:18) {
    individuals_parents <-
      # start from the entire pedigree
      pedigree %>%
      # keep parents of individuals at generation 'k'
      filter(ind %in% individuals$ind) %>%
      # if parent is missing, flag individual as founder
      mutate(founder = ifelse(is.na(father) | is.na(mother), T, F))

    # stop climbing when there are only founders left
    if (nrow(filter(individuals_parents, founder == F)) == 0) {
      break
    }

    # to ascend genealogy, parents become reference individuals
    mothers <- unique(individuals_parents$mother)
    fathers <- unique(individuals_parents$father)
    # recursion step reassigns reference individuals
    individuals <-
      data.frame(ind = c(mothers, fathers), generation = k) %>%
      filter(!is.na(ind)) # removes missing IDs

    # output max depth per individual
    gen_depth <-
      # individuals already visited by climbing algorithm
      gen_depth %>%
      # remove duplicate individuals at 'k-1'
      filter(!ind %in% individuals$ind) %>%
      # append updated individuals 'k' generations deep
      bind_rows(individuals)
  }
  return(gen_depth)
}

# Sum of kinships of ancestors to a list of probands
estimated_genetic_contribution <- function(pedigree, list_of_probands) {
  # probands are individuals at generation 0
  individuals <- data.frame(
    ind = list_of_probands,
    number_of_probands = 1,
    expected_contribution = 1
  )
  # initialize output to append to downstream
  kinship_to_probands <- individuals
  # for loop to iterate through generations
  for (k in 1:18) {
    parents_of_individuals <-
      left_join(individuals, pedigree, by = "ind") %>%
      # if parent is missing, flag individual as founder
      mutate(
        founder = ifelse(is.na(father) | is.na(mother), T, F),
        mother = case_when(is.na(mother) ~ ind, TRUE ~ mother), # revert to selfing if founder
        father = case_when(is.na(father) ~ ind, TRUE ~ father)
      )

    # stop climbing when there are only founders left
    if (nrow(filter(parents_of_individuals, founder == F)) == 0) {
      break
    }

    probands_related_to_mother <-
      parents_of_individuals %>%
      group_by(mother) %>%
      summarise(number_of_probands = sum(number_of_probands)) %>%
      rename(ind = mother) %>%
      ungroup()

    probands_related_to_father <-
      parents_of_individuals %>%
      group_by(father) %>%
      summarise(number_of_probands = sum(number_of_probands)) %>%
      rename(ind = father) %>%
      ungroup()

    # recursion step reassigns reference individuals
    individuals <-
      bind_rows(
        probands_related_to_mother,
        probands_related_to_father
      ) %>%
      mutate(expected_contribution = 2^-(k) * number_of_probands) %>%
      filter(!is.na(ind)) # removes missing IDs

    kinship_to_probands <- bind_rows(kinship_to_probands, individuals)
  }
  estimated_genetic_contribution <-
    kinship_to_probands %>%
    group_by(ind) %>%
    # this sums across generations
    summarise(expected_contribution = sum(expected_contribution))

  return(estimated_genetic_contribution)
}

# Percolate coalescence statistics up the genealogy
coalescence_percolator <- function(pedigree, list_of_probands) {

  # initialize table that will be updated with each iteration
  lineage_table <-
    # maximum genealogical depth for ancestors related to a given list of probands
    maximum_genealogical_depth(pedigree, list_of_probands) %>%
    # append total expected contribution for ancestors related to a given list of probands
    left_join(estimated_genetic_contribution(pedigree, list_of_probands), by = "ind") %>%
    # append 'mother' and 'father' columns (note : left_join restricts to ancestors related to probands)
    left_join(pedigree, by = "ind") %>%
    # initialize coalescence statistics to 0 (exception : probands are ancestral to themselves)
    mutate(
      # p_individual_is_ancestral : p that an individual contains material ancestral to probands
      p_individual_is_ancestral = ifelse(generation == 0, 1, 0),
      # q_no_coalescence : expected number of lineages had there been no coalescence
      q_no_coalescence = 0,
      # coalescence_rate : approximate coalescence rate per individual
      coalescence_rate = 0,
      # estimate how much an individual in the genealogy contributes to kinship among probands
      realized_kinship = 0
    )

  # for loop to iterate through generations
  for (k in 1:max(lineage_table$generation)) {
    # list of ancestors at generation 'k'
    individuals_at_generation_k <- lineage_table %>%
      filter(generation == k) %>%
      pull(ind)

    # get offspring of mothers at generation k
    mothers_at_generation_k <- filter(lineage_table, mother %in% individuals_at_generation_k)

    mothers_realized_kinship <-
      # pairwise combination of all offspring sharing the same mother
      left_join(mothers_at_generation_k, mothers_at_generation_k, by = "mother") %>%
      # limit to distinct offspring
      filter(ind.x != ind.y) %>%
      rowwise() %>%
      # sort names to avoid reciprocal/duplicates
      mutate(name = toString(sort(c(ind.x, ind.y)))) %>%
      # limits to all unique pairs of offspring
      distinct(mother, name) %>%
      # separate into numeric IDs for merging lineage table
      separate(name, into = c("ind1", "ind2"), convert = T) %>%
      # add expected_contribution of offspring 1
      left_join(rename(lineage_table, ind1 = ind), by = c("mother", "ind1")) %>%
      # add expected_contribution of offspring 2
      left_join(rename(lineage_table, ind2 = ind), by = c("mother", "ind2"), suffix = c("1", "2")) %>%
      # coalescence probability is the product of both expected contributions
      mutate(coalescence_probability = expected_contribution1 * expected_contribution2) %>%
      # group by mothers to sum over all pairs of offspring
      group_by(mother) %>%
      # compute realized kinship for mothers at generation k
      summarise(mothers_realized_kinship = 0.25 * sum(coalescence_probability))

    mothers_coalescence_rate <-
      mothers_at_generation_k %>%
      # group offspring by their mothers
      group_by(mother) %>%
      # compute coalescence statistics for mothers at generation k
      summarise(
        p_mother_is_ancestral = 1 - prod((1 - p_individual_is_ancestral) / 2 + 0.5),
        q_mother_no_coalescence = sum(p_individual_is_ancestral / 2),
        mothers_coalescence_rate = q_mother_no_coalescence - p_mother_is_ancestral
      ) %>%
      ungroup() %>%
      left_join(mothers_realized_kinship, by = "mother") %>%
      rename(ind = mother)

    # get offspring of fathers at generation k
    fathers_at_generation_k <- filter(lineage_table, father %in% individuals_at_generation_k)


    fathers_realized_kinship <-
      # pairwise combination of all offspring sharing the same father
      left_join(fathers_at_generation_k, fathers_at_generation_k, by = "father") %>%
      # limit to distinct offspring
      filter(ind.x != ind.y) %>%
      rowwise() %>%
      # sort names to avoid reciprocal/duplicates
      mutate(name = toString(sort(c(ind.x, ind.y)))) %>%
      # limits to all unique pairs of offspring
      distinct(father, name) %>%
      # separate into numeric IDs for merging lineage table
      separate(name, into = c("ind1", "ind2"), convert = T) %>%
      # add expected_contribution of offspring 1
      left_join(rename(lineage_table, ind1 = ind), by = c("father", "ind1")) %>%
      # add expected_contribution of offspring 2
      left_join(rename(lineage_table, ind2 = ind), by = c("father", "ind2"), suffix = c("1", "2")) %>%
      # coalescence probability is the product of both expected contributions
      mutate(coalescence_probability = expected_contribution1 * expected_contribution2) %>%
      # group by fathers to sum over all pairs of offspring
      group_by(father) %>%
      # compute realized kinship for fathers at generation k
      summarise(fathers_realized_kinship = 0.25 * sum(coalescence_probability))

    fathers_coalescence_rate <-
      fathers_at_generation_k %>%
      # group offspring by their fathers
      group_by(father) %>%
      # compute coalescence statistics for fathers at generation k
      summarise(
        p_father_is_ancestral = 1 - prod((1 - p_individual_is_ancestral) / 2 + 0.5),
        q_father_no_coalescence = sum(p_individual_is_ancestral / 2),
        fathers_coalescence_rate = q_father_no_coalescence - p_father_is_ancestral
      ) %>%
      ungroup() %>%
      left_join(fathers_realized_kinship, by = "father") %>%
      rename(ind = father)

    # recursion step updates coalescence statistics in lineage table for individuals at generation k
    lineage_table <-
      lineage_table %>%
      # add coalescence statistics for mothers at generation k
      left_join(mothers_coalescence_rate, by = "ind") %>%
      # add coalescence statistics for fathers at generation k
      left_join(fathers_coalescence_rate, by = "ind") %>%
      # update coalescence statistics for individuals at generation k
      mutate(
        # if individual is a mother or father at generation k, update p_individual_is_ancestral
        p_individual_is_ancestral = case_when(
          !is.na(p_mother_is_ancestral) ~ p_mother_is_ancestral,
          !is.na(p_father_is_ancestral) ~ p_father_is_ancestral,
          TRUE ~ p_individual_is_ancestral
        ),
        # if individual is a mother or father at generation k, update q_no_coalescence
        q_no_coalescence = case_when(
          !is.na(q_mother_no_coalescence) ~ q_mother_no_coalescence,
          !is.na(q_father_no_coalescence) ~ q_father_no_coalescence,
          TRUE ~ q_no_coalescence
        ),
        # if individual is a mother or father at generation k, update coalescence_rate
        coalescence_rate = case_when(
          !is.na(mothers_coalescence_rate) ~ mothers_coalescence_rate,
          !is.na(fathers_coalescence_rate) ~ fathers_coalescence_rate,
          TRUE ~ coalescence_rate
        ),
        # if individual is a mother or father at generation k, update realized_kinship
        realized_kinship = case_when(
          !is.na(mothers_realized_kinship) ~ mothers_realized_kinship,
          !is.na(fathers_realized_kinship) ~ fathers_realized_kinship,
          TRUE ~ realized_kinship
        )
      ) %>%
      # keep original columns of lineage table (removes temp mother/father variables)
      dplyr::select(ind, mother, father, generation, p_individual_is_ancestral, q_no_coalescence, coalescence_rate, expected_contribution, realized_kinship)
  }
  return(lineage_table)
}

# Percolate coalescence statistics up the genealogy
cross_coalescence_percolator <- function(pedigree, lineage_table_A, lineage_table_B) {

  # initialize table that will be updated with each iteration
  lineage_table <-
    bind_rows(lineage_table_A, lineage_table_B) %>%
    # append 'mother' and 'father' columns (note : left_join restricts to ancestors related to probands)
    left_join(pedigree, by = "ind") %>%
    # initialize coalescence statistics to 0
    mutate(realized_kinship = -1)

  # for loop to iterate through generations
  for (k in 1:max(lineage_table$generation)) {
    # print(k)
    # list of ancestors at generation 'k'
    individuals_at_generation_k <- lineage_table %>%
      filter(generation == k) %>%
      pull(ind)

    # get offspring of mothers at generation k
    mothers_at_generation_k <- filter(lineage_table, mother %in% individuals_at_generation_k)

    mothers_realized_kinship <-
      # pairwise combination of all offspring sharing the same mother
      left_join(mothers_at_generation_k, mothers_at_generation_k, by = "mother") %>%
      # limit to distinct offspring
      filter(ind.x != ind.y, region.x != region.y, region.x > region.y) %>%
      # coalescence probability is the product of both expected contributions
      mutate(coalescence_probability = expected_contribution.x * expected_contribution.y) %>%
      # group by mothers to sum over all pairs of offspring
      group_by(mother) %>%
      # compute realized kinship for mothers at generation k
      summarise(mothers_realized_kinship = 0.25 * sum(coalescence_probability)) %>%
      rename(ind = mother)

    # get offspring of fathers at generation k
    fathers_at_generation_k <- filter(lineage_table, father %in% individuals_at_generation_k)


    fathers_realized_kinship <-
      # pairwise combination of all offspring sharing the same father
      left_join(fathers_at_generation_k, fathers_at_generation_k, by = "father") %>%
      # limit to distinct offspring
      filter(ind.x != ind.y, region.x != region.y, region.x > region.y) %>%
      # coalescence probability is the product of both expected contributions
      mutate(coalescence_probability = expected_contribution.x * expected_contribution.y) %>%
      # group by fathers to sum over all pairs of offspring
      group_by(father) %>%
      # compute realized kinship for fathers at generation k
      summarise(fathers_realized_kinship = 0.25 * sum(coalescence_probability)) %>%
      rename(ind = father)

    # recursion step updates coalescence statistics in lineage table for individuals at generation k
    lineage_table <-
      lineage_table %>%
      # add coalescence statistics for mothers at generation k
      left_join(mothers_realized_kinship, by = "ind") %>%
      # add coalescence statistics for fathers at generation k
      left_join(fathers_realized_kinship, by = "ind") %>%
      # update coalescence statistics for individuals at generation k
      mutate(
        # if individual is a mother or father at generation k, update realized_kinship
        realized_kinship = case_when(
          !is.na(mothers_realized_kinship) ~ mothers_realized_kinship,
          !is.na(fathers_realized_kinship) ~ fathers_realized_kinship,
          TRUE ~ realized_kinship
        )
      ) %>%
      # keep original columns of lineage table (removes temp mother/father variables)
      dplyr::select(ind, mother, father, region, generation, expected_contribution, realized_kinship)
  }
  # return output
  return(lineage_table)
}

