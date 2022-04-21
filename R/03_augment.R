# Info --------------------------------------------------------------------
#     Add new variables to your data

# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
treatment <- read_tsv(file = "data/06_treatment_w_meta_clean.tsv")


# Wrangle data ------------------------------------------------------------

# Adding log2 fold change

log_wide <- treatment %>% 
  select(-c(sex, acc_num, age, disease_duration, treatment)) %>% 
  pivot_longer(cols = -id, names_to = "gene") %>% 
  pivot_wider(names_from = id, values_from = value) %>%
  rowwise() %>% 
  mutate("RA_2fc_1" = log2(RA_post_1+1)-log2(RA_pre_1+1),
         "RA_2fc_2" = log2(RA_post_2+1)-log2(RA_pre_2+1),
         "RA_2fc_3" = log2(RA_post_3+1)-log2(RA_pre_3+1),
         "RA_2fc_4" = log2(RA_post_4+1)-log2(RA_pre_4+1),
         "RA_2fc_5" = log2(RA_post_5+1)-log2(RA_pre_5+1),
         "RA_2fc_6" = log2(RA_post_6+1)-log2(RA_pre_6+1),
         "RA_2fc_7" = log2(RA_post_7+1)-log2(RA_pre_7+1),
         "RA_2fc_8" = log2(RA_post_8+1)-log2(RA_pre_8+1),
         "RA_2fc_9" = log2(RA_post_9+1)-log2(RA_pre_9+1),
         "RA_2fc_10" = log2(RA_post_10+1)-log2(RA_pre_10+1),
         "RA_2fc_11" = log2(RA_post_11+1)-log2(RA_pre_11+1),
         "RA_2fc_12" = log2(RA_post_12+1)-log2(RA_pre_12+1),
         "RA_2fc_13" = log2(RA_post_13+1)-log2(RA_pre_13+1),
         "RA_2fc_14" = log2(RA_post_14+1)-log2(RA_pre_14+1),
         "RA_2fc_15" = log2(RA_post_15+1)-log2(RA_pre_15+1),
         "RA_2fc_16" = log2(RA_post_16+1)-log2(RA_pre_16+1),
         "RA_2fc_17" = log2(RA_post_17+1)-log2(RA_pre_17+1),
         "RA_2fc_18" = log2(RA_post_18+1)-log2(RA_pre_18+1),
         "RA_2fc_19" = log2(RA_post_19+1)-log2(RA_pre_19+1)
  )

log_long <- log_wide %>% select(-starts_with("RA_p")) %>%  
  pivot_longer(cols = starts_with("RA_2fc"), names_to = "id") %>% 
  pivot_wider(names_from = "gene", values_from = "value")

# Write data --------------------------------------------------------------

write_tsv(x = log_long, file = "data/03_treatment_log2fc.tsv")