# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
#source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
all_logfold <- read_tsv(file = "data/03_combined_log2fc.tsv")


# Wrangle data ------------------------------------------------------------
my_data_clean_aug %>% ...




# Visualise data ----------------------------------------------------------
my_data_clean_aug %>% ...


# Write data --------------------------------------------------------------
write_tsv(...)
ggsave(...)