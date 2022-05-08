# Load libraries ----------------------------------------------------------
library("tidyverse")
library("ggridges")
library("viridisLite")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
long <- read_tsv(file = "data/02_large_w_meta_clean.tsv")
drug <- read_tsv(file = "data/02_treatment_w_meta_clean.tsv")


# Wrangle data ------------------------------------------------------------
all_metadata <- 
  drug %>%
  select(id, treatment, sex, age) %>% 
  mutate_if(is.logical, as.character) %>% 
  filter(treatment == "TRUE") %>%
  dplyr::rename(disease = treatment) %>% 
  mutate(disease = case_when(disease == 'TRUE' ~ 'Treatment cohort')) %>% 
  bind_rows(long %>% 
              select(id, disease, sex, age) 
  ) 

# Visualise data ----------------------------------------------------------

age_distrib <-
  all_metadata %>% 
  drop_na(age) %>% 
  ggplot(mapping = aes(x = age, fill = disease)) +
  geom_density() +
  facet_wrap(vars(disease), ncol = 1) +
  theme_minimal() +
  scale_fill_viridis_d(alpha = 0.5) +
  theme(legend.position="none") +
  labs(title = "Age distribution per condition")

sex_repres <-
  all_metadata %>% 
  ggplot(mapping = aes(x = sex, fill = sex)) +
  geom_bar(position="dodge") +
  facet_wrap(vars(disease)) +
  theme_minimal() + 
  scale_fill_viridis_d(alpha = 0.8) +
  labs(title = "Sex participation per condition") +
  xlab("") 

# ages of both sex for all conditions
age_n_sex <-
  all_metadata %>% 
  drop_na(age) %>% 
  ggplot(mapping = aes(x = sex, y = age, fill = sex)) + 
  geom_violin(scale = "count", alpha = 0.8, color = NA) +
  geom_jitter(width = 0.1, size = 0.3 ) +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE) +
  facet_wrap(vars(disease)) 

# Write data --------------------------------------------------------------
ggsave("07_age-distribution.png",
       plot = age_distrib,
       path = "/cloud/project/results")
ggsave("07_sex_repres.png",
       plot = sex_repres,
       path = "/cloud/project/results")
ggsave("07_age_n_sex-distribution.png",
       plot = age_n_sex,
       path = "/cloud/project/results")