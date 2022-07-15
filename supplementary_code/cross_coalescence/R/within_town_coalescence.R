#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("Three arguments must be supplied: input_pedigree, town_name, output_directory.n", call.=FALSE)
}

library(data.table)
library(tidyverse)
source("/home/luke1111/projects/ctb-sgravel/luke1111/cross_coalescence/R/coalescence_tools.R")

input_pedigree<-args[1]
town_name<-args[2]
output_directory<-args[3]

pedigree <- fread(input_pedigree)

town_probands <-
  pedigree %>%
  filter(!ind %in% father,        # not fathers
         !ind %in% mother,        # or mothers
         lieum == town_name,      # from focal town
         !is.na(datem),           # who married
         datemp > 1900) %>%       # hard cap at 1900
  pull(ind)

n_probands <- length(town_probands)

if (n_probands < 10) {
  stop("Town must have at least ten probands to climb.n", call.=FALSE)
}

start.time <- Sys.time()
lineage_table <- coalescence_percolator(pedigree, town_probands)
print(Sys.time() - start.time)

lineage_table$n_probands <- n_probands

fwrite(lineage_table, file = paste0(output_directory,"lineage_table_",town_name,".csv"))
