# Script to append superpopulation to sample list downloaded from 1kgp

library(tidyverse)

subject_file <- "data/1kgp/metadata/subjects_populations.tsv"
population_file <- "data/1kgp/metadata/20131219.populations.tsv"
out_file <- "data/1kgp/metadata/unrelated_samples.tsv"

subject_df <- read_tsv(subject_file, col_names = c("subject_name", "population"))
population_df <- read_tsv(population_file) %>%
  rename(population = `Population Code`,
         super_population = `Super Population`) %>%
  select(population, super_population)

subject_df <- subject_df %>%
  left_join(population_df, by = "population")

write_tsv(subject_df, out_file)
