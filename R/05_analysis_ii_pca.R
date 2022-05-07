library(tidyverse)
library(compositions)
library(cowplot)
library(ggplot2)
library(scales)
library(tidyr)
library(broom)

# Reading data
combined <- read_tsv(file = "./data/02_combined_vst.tsv.gz")
combined_meta <- read_tsv(file = "./data/02_combined_meta.tsv")
# Joining data with metadata
combined_w_meta <- combined %>% 
  left_join(combined_meta)

# PCA - differentiated by sex only for healthy and RA
male_combined_w_meta <- combined_w_meta %>% 
  filter(sex == 'M' & (disease == 'Normal' | disease == 'Rheumatoid arthritis (early)' | disease == 'Rheumatoid arthritis (established)'))
female_combined_w_meta <- combined_w_meta %>% 
  filter(sex == 'F' & (disease == 'Normal' | disease == 'Rheumatoid arthritis (early)' | disease == 'Rheumatoid arthritis (established)'))

# Extracting only numeric data for analysis
male_combined_w_meta_numeric <- male_combined_w_meta %>%
  select(!(age)) %>% 
  select(where(is.numeric)) 

# Extracting non-zero data and performing PCA
male_pca_fit <- male_combined_w_meta_numeric[, which(colSums(abs(male_combined_w_meta_numeric)) != 0)]  %>%
  prcomp() # do PCA on scaled data

# Extracting only numeric data for analysis
female_combined_w_meta_numeric <- female_combined_w_meta %>%
  select(!(age)) %>% 
  select(where(is.numeric))    # retain only numeric columns

# Extracting non-zero data and performing PCA
female_pca_fit <- female_combined_w_meta_numeric[, which(colSums(abs(female_combined_w_meta_numeric)) != 0)] %>%
  prcomp() # do PCA on scaled data

# Ploting PCAs
male_pca_plot <- male_pca_fit %>%
  augment(male_combined_w_meta) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, 
             .fittedPC2, 
             color = disease)) + 
  geom_point(size = 1.5) +
  labs(x = 'PC1',
       y = 'PC2',
       title = 'PCA plot of diseases in males') +
  theme_minimal()

female_pca_plot <- female_pca_fit %>%
  augment(female_combined_w_meta) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, 
             .fittedPC2, 
             color = disease)) + 
  geom_point(size = 1.5) +
  labs(x = 'PC1',
       y = 'PC2',
       title = 'PCA plot of diseases in females') +
  theme_minimal()

# Extracting a legend
legend <- get_legend(male_pca_plot + theme(legend.position = 'bottom'))

prow <- plot_grid(male_pca_plot + theme(legend.position = 'null') + theme_cowplot(), 
                  female_pca_plot + theme(legend.position = 'null'), 
                  rel_widths = c(2, 2))
two_pcas <- plot_grid(prow, legend, ncol = 1, rel_heights = c(1, .1))


## KMEANS analysis - females
# Tidying PC analysis data
tidy_female_pca <- female_pca_fit %>% 
  tidy() 

# Extracting first two components
female_points <- tidy_female_pca %>% 
  filter(PC==1 | PC==2) %>% 
  pivot_wider(names_from = 'PC', values_from = 'value') %>% 
  select(-row)

# Renaming columns
female_points <- female_points %>% 
  rename(PC1 = "1") %>% 
  rename(PC2 = "2")

# Performing kmeans clustering for 2 and 3 clusters
kclusts_female <- 
  tibble(k = 2:3) %>%
  mutate(
    kclust_female = map(k, ~kmeans(female_points, .x)),
    tidied_female = map(kclust_female, tidy),
    glanced_female = map(kclust_female, glance),
    augmented_female = map(kclust_female, augment, female_points)
  )

clusters_female <- kclusts_female %>%
  unnest(cols = c(tidied_female))

assignments_female <- kclusts_female %>% 
  unnest(cols = c(augmented_female))

clusterings_female <- kclusts_female %>%
  unnest(cols = c(glanced_female))

# Ploting clustering's results
kmeans_plot_female <- assignments_female %>% 
  ggplot(aes(x = PC1, 
             y = PC2)) +
  geom_point(aes(color = .cluster), 
             alpha = 0.8) +
  labs(title = 'Plot of clustered PCA components',
       color = 'cluster') +
  facet_wrap(~k) +
  theme_minimal()

## KMEANS analysis - females
# Tidying PC analysis data
tidy_male_pca <- male_pca_fit %>% 
  tidy() 

# Extracting first two components
male_points <- tidy_male_pca %>% 
  filter(PC==1 | PC==2) %>% 
  pivot_wider(names_from = 'PC', values_from = 'value') %>% 
  select(-row)

# Renaming columns
male_points <- male_points %>% 
  rename(PC1 = "1") %>% 
  rename(PC2 = "2")

# Performing kmeans clustering for 2 and 3 clusters
kclusts_male <- 
  tibble(k = 2:3) %>%
  mutate(
    kclust_male = map(k, ~kmeans(male_points, .x)),
    tidied_male = map(kclust_male, tidy),
    glanced_male = map(kclust_male, glance),
    augmented_male = map(kclust_male, augment, male_points)
  )

clusters_male <- kclusts_male %>%
  unnest(cols = c(tidied_male))

assignments_male <- kclusts_male %>% 
  unnest(cols = c(augmented_male))

clusterings <- kclusts_male %>%
  unnest(cols = c(glanced_male))

# Ploting clustering's results
kmeans_plot_male <- assignments_male %>% 
  ggplot(aes(x = PC1, 
             y = PC2)) +
  geom_point(aes(color = .cluster), 
             alpha = 0.8) +
  labs(title = 'Plot of clustered PCA components',
       color = 'cluster') +
  facet_wrap(~k) +
  theme_minimal()


# Plotting combined plots
prow_clus_female <- plot_grid(female_pca_plot + 
                                theme(legend.position = 'null') + 
                                guides(colour = guide_legend(nrow = 3)), 
                              kmeans_plot_female + 
                                theme(legend.position = 'null') + 
                                guides(colour = guide_legend(nrow = 3)), 
                              rel_widths = c(2, 4))

prow_clus_male <- plot_grid(male_pca_plot + 
                              theme(legend.position = 'bottom') + 
                              guides(colour = guide_legend(nrow = 3)), 
                            kmeans_plot_male + 
                              theme(legend.position = 'bottom') + 
                              guides(colour = guide_legend(nrow = 3)), 
                            rel_widths = c(2, 4))

pcol_clus <- plot_grid(prow_clus_female, 
                       prow_clus_male, 
                       nrow = 2,
                       rel_heights = c(1, 1.3))

pcol_clus

ggsave('results/pcas.png', pcol_clus, width = 15, height = 8, bg = "transparent")
