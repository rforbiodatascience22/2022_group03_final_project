# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
large_w_meta <- read_tsv(file = "./data/02_large_w_meta_clean.tsv")


# Wrangle data ------------------------------------------------------------
large_mod <- large_w_meta %>% 
  select(!c(sex, age, acc_num)) %>% 
  filter(disease == c("Normal","Rheumatoid arthritis (early)")) %>% 
  mutate(outcome = case_when(disease == "Normal" ~ 0,
                             disease == "Rheumatoid arthritis (early)" ~ 1)) %>% 
  relocate(outcome)

# Model data
earlyRA.model <- large_mod %>%
  glm(outcome ~ LXN + IL8,
      data = .,
      family = binomial(link = "logit"))


# Visualise data ----------------------------------------------------------
my_data_clean_aug %>% ...


# Write data --------------------------------------------------------------
write_tsv(...)
ggsave(...)