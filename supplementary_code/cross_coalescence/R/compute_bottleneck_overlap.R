#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=4) {
  stop("Four arguments must be supplied: town_A, town_B, input_directory, output_directory.n", call.=FALSE)
}

library(data.table)
library(tidyverse)

A <- args[1]
B <- args[2]
input_directory <- args[3]
output_directory <- args[4]

in_fileA <- paste0(input_directory, "lineage_table_", A, ".csv")
in_fileB <- paste0(input_directory, "lineage_table_", B, ".csv")
cross_A_B <- paste0(input_directory, "lineage_table_", A, "_", B, ".csv")

if (!file.exists(in_fileA)) {
  print(paste(in_fileA, "does not exist"))
  break
}

if (!file.exists(in_fileB)) {
  print(paste(in_fileB, "does not exist"))
  break
}

if (!file.exists(cross_A_B)) {
  print(paste(cross_A_B, "does not exist"))
  break
}

lineage_table_A <- fread(in_fileA) %>%
  mutate(location = "kinship_A") %>%
  filter(!is.na(ind), generation > 0, realized_kinship > 0) %>%
  distinct(ind, generation, realized_kinship, location)

lineage_table_B <- fread(in_fileB) %>%
  mutate(location = "kinship_B") %>%
  filter(!is.na(ind), generation > 0, realized_kinship > 0) %>%
  distinct(ind, generation, realized_kinship, location)

lineage_table_A_B <- fread(cross_A_B) %>%
  mutate(location = "kinship_A_B") %>%
  filter(!is.na(ind), generation > 0, realized_kinship > 0) %>%
  distinct(ind, generation, realized_kinship, location)

town_A_B_overlap <- 
  bind_rows(lineage_table_A, 
            lineage_table_B, 
            lineage_table_A_B) %>%
  group_by(location) %>%
  dplyr::summarise(total_realized_kinship = sum(realized_kinship, na.rm = T)) %>%
  pivot_wider(names_from = location, values_from = total_realized_kinship) %>%
  mutate(percent_overlap = kinship_A_B / (kinship_A + kinship_B) * 100,
         town_A = A, town_B = B)
  
fwrite(town_A_B_overlap, file = paste0(output_directory,"overlap_", A, "_", B, ".csv"))