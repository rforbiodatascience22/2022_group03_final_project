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
large_raw <- read_tsv(file = "./data/01_large_count.tsv")
treatment_raw <- read_tsv(file = "./data/01_treatment_count.tsv")
large_meta_raw <- read_tsv(file = "data/01_large_meta.tsv")
treatment_meta_raw <- read_tsv(file = "data/01_treatment_meta.tsv")

# Wrangle data ------------------------------------------------------------

# Large dataset

# Count data

large <- large_raw %>% 
  pivot_longer(cols = contains("tissue"), names_to = "id") %>% 
  pivot_wider(names_from = ...1, values_from = value)

# Meta data

large_meta <- large_meta_raw %>% 
  slice(1,9,10,11) %>% # Only take meaningful values
  select(-("!Sample_title")) %>% # Drop title column
  mutate(var_num = row_number()) %>% # Add index column instead
  pivot_longer(cols = contains("tissue"), names_to = "sample") %>% 
  pivot_wider(names_from = var_num, values_from = value) %>% 
  separate("2", sep = " ", into = c("o","sex")) %>% # Extract sex
  separate("3", sep = "age: ", into = c("o","age_str")) %>% # Extract age
  mutate(age = as.integer(age_str)) %>% # Convert age to integer
  separate("4", sep = "disease: ", into = c("o","disease")) %>% 
  rename(acc_num = "1") %>% # Add title to accession number column
  bind_cols(select(large,id)) %>% # Add correct id's (I have checked, they match)
  select(id, disease, sex, age, acc_num) # Extract final columns

# Join

large_w_meta <- right_join(large_meta, large, by = "id")

# Treatment dataset

# Count data

treatment <- treatment_raw %>% 
  pivot_longer(cols = starts_with("RA"), names_to = "id") %>% 
  pivot_wider(names_from = ...1, values_from = value)

# Meta data

treatment_meta <- treatment_meta_raw %>% 
  slice(1,9,10,12) %>% # Only take meaningful values
  select(-("!Sample_title")) %>% # Drop title column
  mutate(var_num = row_number()) %>% # Add index column instead
  pivot_longer(cols = starts_with("rheumatoid"), # Transpose data 1.
               names_to = "sample") %>% 
  pivot_wider(names_from = "var_num", values_from = "value") %>% # Transpose data 2.
  separate("2", sep = " ", into = c("o","sex")) %>% # Extract sex
  separate("3", sep = "age: ", into = c("o","age_str")) %>% # Extract age
  mutate(age = as.integer(age_str)) %>% # Convert age to integer
  separate("4", sep = "weeks: ", into = c("o","disease_duration_str")) %>% # Extract disease duration
  mutate(disease_duration = as.integer(disease_duration_str)) %>% # Convert duration to integer
  rename(acc_num = "1") %>% # Add title to accession number column
  mutate(id = str_extract(sample,"RA_[a-z]{3,4}_[[:digit:]]{1,2}")) %>% # Extract sample ID
  fill(age, .direction = "up") %>% # Fill age to paired sample
  fill(disease_duration, .direction = "up") %>% # Fill duration to paired sample
  mutate(treatment = str_detect(id, "post")) %>% # Add treatment column
  select(id,treatment,sex,age,disease_duration,acc_num) # Extract final columns


# Join

treatment_w_meta <- right_join(treatment_meta, treatment, by = "id")

# Write data --------------------------------------------------------------
write_tsv(x = large_w_meta,
          file = "data/02_large_w_meta_clean.tsv")

write_tsv(x = treatment_w_meta,
          file = "data/02_treatment_w_meta_clean.tsv")
