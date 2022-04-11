# Info --------------------------------------------------------------------
#     Collapse data to a single file or convert .xlsx to .tsv
#     here we could imagine having an .xlsx-file with multiple sheets,
#     from which we create a single .tsv

# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
# my_data_raw <- read_tsv(file = "data/_raw/my_raw_data.tsv")
genecount_i_raw <- read_tsv(
  "data/_raw/GSE89408_GEO_count_matrix_rename.txt.gz"
)

genecount_treatment_i_raw <- read_tsv(
  "data/_raw/GSE97165_GEO_count_matrix.txt.gz"
)

metadata_treatment_i_raw <- read_csv(
  "data/_raw/metadata_treatment_i.txt"
)

metadata_i_raw <- read_csv(
  "data/_raw/metadata_i.txt"
)

# Wrangle data ------------------------------------------------------------
# my_data <- my_data_raw # %>% ...


# Write data --------------------------------------------------------------
# write_tsv(x = my_data,
#          file = "data/01_my_data.tsv")
