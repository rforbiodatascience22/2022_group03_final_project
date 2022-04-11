# Info --------------------------------------------------------------------
#     Remove invalid data, e.g. if you have amino acid sequence data,
#     remove non-valid sequences containing X or other non-standard amino acid
#     characters or fix columns,
#     e.g. dates or when two labels are the same, but spelled differently

# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
my_data <- read_tsv(file = "data/01_my_data.tsv")


# Wrangle data ------------------------------------------------------------
my_data_clean <- my_data # %>% ...


# Write data --------------------------------------------------------------
write_tsv(x = my_data_clean,
          file = "data/02_my_data_clean.tsv")