# Load libraries ----------------------------------------------------------
library("tidyverse")
library("broom")
library("compositions")
library("preprocessCore")

# Define functions --------------------------------------------------------
#source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
large <- read_tsv(file = "data/02_large_w_meta_clean.tsv")
treatment <- read_tsv(file = "data/02_treatment_w_meta_clean.tsv")

# Wrangle data ------------------------------------------------------------

# Join the two data sets:

large_no_meta_wide <- large %>% 
  select(-c(disease,sex,acc_num,age)) %>% 
  pivot_longer(cols = -id, names_to = "gene") %>% 
  pivot_wider(names_from = id, values_from = value)

treatment_no_meta_wide <- treatment %>%
  select(-c(treatment,sex,age,disease_duration, acc_num)) %>% 
  pivot_longer(cols = -id, names_to = "gene") %>% 
  pivot_wider(names_from = id, values_from = value)

combined <- treatment_no_meta_wide %>% 
  inner_join(large_no_meta_wide, by = "gene") %>% 
  pivot_longer(cols = -gene, names_to = "id") %>% 
  pivot_wider(names_from = gene, values_from = value)

combined_RA_normal_filtered <- treatment_no_meta_wide %>% 
  inner_join(large_no_meta_wide, by = "gene") %>% 
  select(gene, starts_with("RA_p"), starts_with("normal")) %>%
  rowwise() %>% 
  filter(sum(c_across(-gene)) > 100) %>%  
  pivot_longer(cols = -gene, names_to = "id") %>% 
  pivot_wider(names_from = gene, values_from = value)

combined_RA_normal_filtered_wide <- treatment_no_meta_wide %>% 
  inner_join(large_no_meta_wide, by = "gene") %>% 
  select(gene, starts_with("RA_p"), starts_with("normal")) %>%
  rowwise() %>% 
  filter(sum(c_across(-gene)) > 100)

combined_RA_normal_filtered_wide_integer <- 
  combined_RA_normal_filtered_wide %>% 
  mutate(across(where(is.numeric), ~ as.integer(.x)))

vst_counts <- combined_RA_normal_filtered_wide_integer %>% 
  select(-gene) %>% 
  as.matrix() %>% 
  unname() %>%
  DESeq2::vst()

colnames(vst_counts) <- pull(combined_RA_normal_filtered, id)
rownames(vst_counts) <- pull(combined_RA_normal_filtered_wide_integer, gene)

vst_complete <- vst_counts %>% 
  as_tibble(rownames = "gene") %>% 
  pivot_longer(cols = -gene, names_to = "id") %>% 
  pivot_wider(names_from = gene, values_from = value)


rlog_counts <- combined_RA_normal_filtered_wide_integer %>% 
  select(-gene) %>% 
  as.matrix() %>% 
  unname() %>%
  DESeq2::rlog()

colnames(rlog_counts) <- pull(combined_RA_normal_filtered, id)
rownames(rlog_counts) <- pull(combined_RA_normal_filtered_wide_integer, gene)

rlog_complete <- rlog_counts %>% 
  as_tibble(rownames = "gene") %>% 
  pivot_longer(cols = -gene, names_to = "id") %>% 
  pivot_wider(names_from = gene, values_from = value)

write_tsv(rlog_complete, file = "./experiments/rlog_treatment_vs_normal.tsv")

# -------------------------------------------------------------------------
###### Make linear model. Predict early RA tissue from normal:

rlog_complete <- read_tsv(file = "./experiments/rlog_treatment_vs_normal.tsv")

early_vs_normal <- rlog_complete %>% 
  filter(str_detect(id, "RA_post", negate = TRUE)) %>% 
  mutate(response = case_when(str_detect(id,"normal") ~ 0, 
                              TRUE ~ 1)) %>% 
  select(-id)

long <- early_vs_normal %>% 
  pivot_longer(cols = -response, 
               names_to = "gene", 
               values_to = "log_expr")

nested <- long %>% 
  group_by(gene) %>% 
  nest() %>% 
  ungroup()

mdl_data <- nested %>% mutate(mdl = map(data, 
                                        ~glm(response ~ log_expr, 
                                             data = ., 
                                             family = binomial(link = "logit"))))

mdl_tidy <- mdl_data %>% 
  mutate(tidy = map(mdl, tidy, conf.int = TRUE)) %>% 
  unnest(tidy)

mdl_coef <- mdl_tidy %>% filter(term == "log_expr")

mdl_coef %>% 
  select(-c(data, mdl)) %>% 
  write_tsv(file = "./experiments/mdl_coef.tsv")


mdl_sign <- mdl_coef %>% filter(p.value < (0.05/200)) %>% 
  mutate(significant = if_else(p.value < (0.05/300), 1, 0))

mdl_plot <- mdl_sign %>% mutate(neg_log10 = -log10(p.value))


# Plot
ggplot(data = mdl_plot, 
       mapping = aes(x = fct_rev(fct_reorder(gene,neg_log10)),
                     y = neg_log10,
                     color = factor(significant))) +
  geom_point() +
  geom_abline(intercept = -log10(0.05/300), slope = 0)

ggplot(data = mdl_plot, mapping = aes(x = estimate, 
                                      y = fct_rev(fct_reorder(gene, estimate)),
                                      color = factor(significant))) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 0, linetype = 2)





# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
###### Batch effect. Combined dataset colored by sample.

combined <- read_tsv(file = "data/02_combined_vst.tsv.gz") %>% 
  filter(str_detect(id, "OA|AG|undiff", negate = TRUE))

# Model data
pca_fit <- combined %>%
  select(-c(id)) %>%
  prcomp(scale = TRUE)

# Augment

aug <- combined %>% mutate(sample = str_detect(id, "RA_p"))

# Visualise data

pca_fit %>%
  augment(aug) %>% 
  ggplot(mapping = aes(x = .fittedPC1, 
                       y =.fittedPC3,
                       color = sample)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  labs(color = "Sample")

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
###### Fix batch effect.

combined <- read_tsv(file = "data/02_combined_vst.tsv.gz") %>% 
  filter(str_detect(id, "OA|AG|undiff", negate = TRUE))

# Model data
pca_fit <- combined %>%
  select(-c(id)) %>%
  prcomp()

# Augment

aug <- combined %>% mutate(sample = str_detect(id, "RA_p"))

# Visualise data

pca_fit %>%
  augment(aug) %>% 
  ggplot(mapping = aes(x = .fittedPC1, 
                       y =.fittedPC3,
                       color = sample)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  labs(color = "Sample")



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
###### Make linear model. Predict RA_tissue from normal (Large set):

# Filter heavily, and subsample RA to 28
combined <- read_tsv(file = "./data/02_combined.tsv") %>% 
  filter(str_detect(id, "RA_tissue|normal_tissue"))

combined_filtered_wide <- combined %>% 
  pivot_longer(cols = -id, names_to = "gene") %>% 
  pivot_wider(names_from = id, values_from = value) %>% 
  rowwise() %>% 
  filter(sum(c_across(-gene)) > 200) %>% 
  ungroup()

combined_filtered_wide_vst <- combined_filtered_wide %>% 
  mutate(across(where(is.numeric), ~ as.integer(.x))) %>% 
  select(-gene) %>% 
  as.matrix() %>% 
  unname() %>%
  DESeq2::vst()

colnames(combined_filtered_wide_vst) <- pull(combined, id)
rownames(combined_filtered_wide_vst) <- pull(combined_filtered_wide, gene)


combined_vst_long <- combined_filtered_wide_vst %>% 
  as_tibble(rownames = "gene") %>% 
  pivot_longer(cols = -gene, names_to = "id") %>% 
  pivot_wider(names_from = gene, values_from = value)

View(combined_vst_long)

# Subsample RA_tissue
set.seed(18)
RA_tissue <- combined_vst_long %>% 
  filter(str_detect(id, "RA_tissue")) %>%
  sample_n(28)


# Combine with healthy tissue:
RA_vs_normal <- combined_vst_long %>% 
  filter(str_detect(id, "normal")) %>% 
  bind_rows(RA_tissue) %>% 
  mutate(response = if_else(str_detect(id, "RA"), 1, 0))


# Create long form for modelling:
long <- RA_vs_normal %>% 
  select(-id) %>% 
  pivot_longer(cols = -response, 
               names_to = "gene", 
               values_to = "log_expr")

# Nest data:
nested <- long %>% 
  group_by(gene) %>% 
  nest() %>% 
  ungroup()

# Model:
mdl_data <- nested %>% mutate(mdl = map(data, 
                                        ~glm(response ~ log_expr, 
                                             data = ., 
                                             family = binomial(link = "logit"))))
# Tidy:
mdl_tidy <- mdl_data %>% 
  mutate(tidy = map(mdl, tidy, conf.int = TRUE)) %>% 
  unnest(tidy)



# Extract only coefficient:
mdl_coef <- mdl_tidy %>% 
  filter(term == "log_expr") %>% 
  select(-c(data,mdl))

write_tsv(mdl_coef, file = "experiments/model_coeff.tsv")

# Plot
View(mdl_tidy)
# Gene on x, coeff on y:




# -------------------------------------------------------------------------
###### Make linear model. Predict RA_tissue from normal (Large set):

vst_complete <- read_tsv(file = "./data/02_combined_vst.tsv.gz")

RA_vs_normal <- vst_complete %>% 
  filter(str_detect(id, "RA_tissue|normal_tissue")) %>%  
  mutate(response = if_else(str_detect(id, "RA"), 1, 0)) %>% 
  select(-id)

long <- RA_vs_normal %>% 
  pivot_longer(cols = -response, 
               names_to = "gene", 
               values_to = "log_expr")

nested <- long %>% 
  group_by(gene) %>% 
  nest() %>% 
  ungroup()

mdl_data <- nested %>% mutate(mdl = map(data, 
                                        ~glm(response ~ log_expr, 
                                             data = ., 
                                             family = binomial(link = "logit"))))

mdl_tidy <- mdl_data %>% 
  mutate(tidy = map(mdl, tidy, conf.int = TRUE)) %>% 
  unnest(tidy)

mdl_coef <- mdl_tidy %>% filter(term == "log_expr")

mdl_coef %>% 
  select(-c(data, mdl)) %>% 
  write_tsv(file = "./experiments/mdl_coef.tsv")


mdl_sign <- mdl_coef %>% filter(p.value < (0.05/200)) %>% 
  mutate(significant = if_else(p.value < (0.05/300), 1, 0))

mdl_plot <- mdl_sign %>% mutate(neg_log10 = -log10(p.value))


# Plot
ggplot(data = mdl_plot, 
       mapping = aes(x = fct_rev(fct_reorder(gene,neg_log10)),
                     y = neg_log10,
                     color = factor(significant))) +
  geom_point() +
  geom_abline(intercept = -log10(0.05/300), slope = 0)

ggplot(data = mdl_plot, mapping = aes(x = estimate, 
                                      y = fct_rev(fct_reorder(gene, estimate)),
                                      color = factor(significant))) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 0, linetype = 2)




# -------------------------------------------------------------------------
###### Diseases (clr()):

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
###### Diseases (normalize.quantiles()):

View(large)

# Model data
pca_fit <- large %>%
  filter(disease != "Osteoarthritis") %>% 
  select(-c(id,disease,sex,acc_num,age)) %>%
  as.matrix() %>% 
  normalize.quantiles() %>% 
  prcomp(scale = TRUE)

# Change augment data:

large_aug <- large %>% 
  filter(disease != "Osteoarthritis")

# Visualise data
pca_fit %>%
  augment(large_aug) %>% 
  ggplot(mapping = aes(x = .fittedPC1, 
                       y =.fittedPC2, 
                       color = factor(disease))) +
  geom_point() +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  labs(color = 'disease')

# -------------------------------------------------------------------------
###### Diseases vs. normal (combined) (normalize.quantiles()):


# Model data
pca_fit <- combined_RA_normal %>%
  filter(str_detect(id, "^RA_p") | str_detect(id, "normal")) %>% 
  select(-c(id)) %>%
  as.matrix() %>% 
  limma::normalizeQuantiles() %>% 
  prcomp(scale = TRUE)

# Change augment data:

combined_aug <- combined_RA_normal %>% 
  filter(str_detect(id, "^RA_p") | str_detect(id, "normal")) %>% 
  mutate(normal = str_detect(id, "normal"))

# Visualise data
pca_fit %>%
  augment(combined_aug) %>% 
  ggplot(mapping = aes(x = .fittedPC1, 
                       y =.fittedPC2,
                       color = normal)) +
  geom_point() +
  geom_text(aes(label = id)) +
  theme_minimal() +
  theme(legend.position = 'bottom')

# -------------------------------------------------------------------------
###### Diseases vs. normal (combined and filtered) (vst()):


# Model data
pca_fit <- vst_complete %>%
  filter(str_detect(id, "^RA_p") | str_detect(id, "normal")) %>% 
  select(-c(id)) %>%
  as.matrix() %>%
  prcomp(scale = TRUE)

# Change augment data:

vst_aug <- vst_complete %>% 
  filter(str_detect(id, "^RA_p") | str_detect(id, "normal")) %>% 
  mutate(normal = str_detect(id, "normal"))

# Visualise data
pca_fit %>%
  augment(combined_aug) %>% 
  ggplot(mapping = aes(x = .fittedPC1, 
                       y =.fittedPC2,
                       color = normal)) +
  geom_point() +
  geom_text(aes(label = id)) +
  theme_minimal() +
  theme(legend.position = 'bottom')

# -------------------------------------------------------------------------
###### Diseases vs. normal (combined and filtered) (rlog()):


# Model data
pca_fit <- rlog_complete %>%
  filter(str_detect(id, "^RA_p") | str_detect(id, "normal")) %>% 
  select(-c(id)) %>%
  as.matrix() %>%
  prcomp(scale = TRUE)

# Change augment data:

rlog_aug <- rlog_complete %>% 
  filter(str_detect(id, "^RA_p") | str_detect(id, "normal")) %>% 
  mutate(normal = str_detect(id, "normal"))

# Visualise data
pca_fit %>%
  augment(rlog_aug) %>% 
  ggplot(mapping = aes(x = .fittedPC1, 
                       y =.fittedPC2,
                       color = normal)) +
  geom_point() +
  theme_minimal() +
  theme(legend.position = 'bottom')

# -------------------------------------------------------------------------
###### Normal vs. RA large vst():

combined_meta <- read_tsv(file = "./data/02_combined_meta.tsv")

vst_complete <- read_tsv(file = "./data/02_combined_vst.tsv")

normal_vs_RA_vst <- vst_complete %>% 
  filter(str_detect(id, "normal|RA_tissue")) %>%
  mutate(normal = str_detect(id,"normal", negate = TRUE))

# Model data
pca_fit <- normal_vs_RA_vst %>%
  select(-c(id, normal)) %>%
  as.matrix() %>%
  prcomp()

# Augment metadata:

augment_data <- normal_vs_RA_vst %>% left_join(combined_meta, by = "id")

# Visualise data
pca_fit %>%
  augment(augment_data) %>% 
  ggplot(mapping = aes(x = .fittedPC1, 
                       y =.fittedPC2,
                       color = sex)) +
  geom_point() +
  geom_text(aes(label = sex)) +
  theme_minimal() +
  theme(legend.position = 'bottom')

# -------------------------------------------------------------------------
###### RA large. Male vs. female vst():

combined_meta <- read_tsv(file = "./data/02_combined_meta.tsv")

vst_complete <- read_tsv(file = "./data/02_combined_vst.tsv")

RA_vst <- vst_complete %>% 
  filter(str_detect(id, "RA_tissue"))

# Model data
pca_fit <- RA_vst %>%
  select(-c(id)) %>%
  as.matrix() %>%
  prcomp()

# Augment metadata:

augment_data <- RA_vst %>% left_join(combined_meta, by = "id")

# Visualise data
pca_fit %>%
  augment(augment_data) %>% 
  ggplot(mapping = aes(x = .fittedPC1, 
                       y =.fittedPC2,
                       color = sex)) +
  geom_point() +
  geom_text(aes(label = sex)) +
  theme_minimal() +
  theme(legend.position = 'bottom')

# -------------------------------------------------------------------------
###### RA_post. Male vs. female vst():

combined_meta <- read_tsv(file = "./data/02_combined_meta.tsv")

vst_complete <- read_tsv(file = "./data/02_combined_vst.tsv")

RA_post_vst <- vst_complete %>% 
  filter(str_detect(id, "RA_post"))

# Model data
pca_fit <- RA_post_vst %>%
  select(-c(id)) %>%
  as.matrix() %>%
  prcomp()

# Augment metadata:

augment_data <- RA_post_vst %>% left_join(combined_meta, by = "id")

# Visualise data
pca_fit %>%
  augment(augment_data) %>% 
  ggplot(mapping = aes(x = .fittedPC1, 
                       y =.fittedPC2,
                       color = sex)) +
  geom_point() +
  geom_text(aes(label = sex)) +
  theme_minimal() +
  theme(legend.position = 'bottom')


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