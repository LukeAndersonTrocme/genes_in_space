#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 4) {
  stop("Four arguments must be supplied: pedigree, town_A, town_B and output_directory.n", call. = FALSE)
}

library(data.table)
library(tidyverse)
source("/home/luke1111/projects/ctb-sgravel/luke1111/cross_coalescence/R/coalescence_tools.R")

input_pedigree<-args[1]
A <- args[2]
B <- args[3]
output_directory <- args[4]

in_fileA <- paste0(output_directory, "lineage_table_", A, ".csv")
in_fileB <- paste0(output_directory, "lineage_table_", B, ".csv")

if (!file.exists(in_fileA)) {
  print(paste(in_fileA, "does not exist"))
  break
}

if (!file.exists(in_fileB)) {
  print(paste(in_fileB, "does not exist"))
  break
}

lineage_table_A <- fread(in_fileA) %>%
  mutate(region = A) %>%
  dplyr::select(ind, region, generation, expected_contribution)

lineage_table_B <- fread(in_fileB) %>%
  mutate(region = B) %>%
  dplyr::select(ind, region, generation, expected_contribution)

pedigree <- fread(input_pedigree)

start.time <- Sys.time()
lineage_table_A_B <- cross_coalescence_percolator(pedigree, lineage_table_A, lineage_table_B)
print(Sys.time() - start.time)

lineage_table_A_B$n_probands_A <- unique(lineage_table_A$n_probands)
lineage_table$n_probands_B <- unique(lineage_table_B$n_probands)

fwrite(lineage_table_A_B, file = paste0(output_directory, "lineage_table_", A, "_", B, ".csv"))
