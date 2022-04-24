# Load libraries ----------------------------------------------------------
library("tidyverse")
library("broom")
library("compositions")

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
large <- read_tsv(file = "data/02_large_w_meta_clean.tsv")


# Wrangle data ------------------------------------------------------------



# -------------------------------------------------------------------------
###### Diseases:

# Model data
pca_fit <- large %>%
  select(-c(id,disease,sex,acc_num,age)) %>%
  clr() %>% 
  prcomp(scale = TRUE)


# Visualise data
pca_fit %>%
  augment(large) %>% 
  ggplot(mapping = aes(x = .fittedPC1, 
                       y =.fittedPC2, 
                       color = factor(disease))) +
  geom_point() +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  labs(color = 'disease')



# -------------------------------------------------------------------------
###### Male vs. Female (diseased):

# Model data

large_no_normal <- large %>%
  filter(disease != "Normal")

pca_fit <- large_no_normal %>%
  select(-c(id,disease,sex,acc_num,age)) %>%
  clr() %>% 
  prcomp(scale = TRUE)


# Visualise data
pca_fit %>% 
  augment(large_no_normal) %>% 
  ggplot(mapping = aes(x = .fittedPC1, 
                       y =.fittedPC2, 
                       color = factor(sex))) +
  geom_point() +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  labs(color = 'sex')



# -------------------------------------------------------------------------
###### Age group (diseased):

# Model data

large_no_normal <- large %>%
  filter(disease != "Normal")

pca_fit <- large_no_normal %>%
  select(-c(id,disease,sex,acc_num,age)) %>%
  clr() %>% 
  prcomp(scale = TRUE)


# Visualise data
pca_fit %>% 
  augment(large_no_normal) %>% 
  ggplot(mapping = aes(x = .fittedPC1, 
                       y =.fittedPC2, 
                       color = age)) +
  geom_point() +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  labs(color = 'age') +
  scale_colour_gradient(low = "white", high = "black", na.value = NA)

# Write data --------------------------------------------------------------
#write_tsv(...)
#ggsave(...)