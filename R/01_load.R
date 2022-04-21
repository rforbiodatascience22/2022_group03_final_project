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
genecount_large_load <- read_tsv(
  "data/_raw/GSE89408_GEO_count_matrix_rename.txt.gz"
)

genecount_treatment_load <- read_tsv(
  "data/_raw/GSE97165_GEO_count_matrix.txt.gz"
)

metadata_large_load <- read_tsv(
  "data/_raw/meta_large_GSE89408.txt",
  skip = 29, # Skip 29 first rows in data.
  na = c("", "NA", "n/a")
)


metadata_treatment_load <- read_tsv(
  "data/_raw/meta_treatment_GSE97165.txt",
  skip = 28, # Skip 28 first rows in data.
  na = c("", "NA", "n/a")
  )

# Wrangle data ------------------------------------------------------------
# my_data <- my_data_raw # %>% ...


# Write data --------------------------------------------------------------
write_tsv(x = genecount_large_load,
          file = "data/01_large_count.tsv")
write_tsv(x = genecount_treatment_load,
          file = "data/01_treatment_count.tsv")
write_tsv(x = metadata_large_load,
          file = "data/01_large_meta.tsv")
write_tsv(x = metadata_treatment_load,
          file = "data/01_treatment_meta.tsv")
