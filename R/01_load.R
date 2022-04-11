# Info --------------------------------------------------------------------
#     Collapse data to a single file or convert .xlsx to .tsv
#     here we could imagine having an .xlsx-file with multiple sheets,
#     from which we create a single .tsv

# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
my_data_raw <- read_tsv(file = "data/_raw/my_raw_data.tsv")


# Wrangle data ------------------------------------------------------------
my_data <- my_data_raw # %>% ...


# Write data --------------------------------------------------------------
write_tsv(x = my_data,
          file = "data/01_my_data.tsv")
