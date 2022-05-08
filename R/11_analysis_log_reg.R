# Load libraries ----------------------------------------------------------
library("tidyverse")
library("broom")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------

# Increase connection size:
Sys.setenv("VROOM_CONNECTION_SIZE" = 2*131072)

# Load data:
vst_complete <- read_tsv(file = "./data/02_combined_vst.tsv.gz")

# Wrangle data ------------------------------------------------------------

# Extract samples, create response variable and convert to long form:

normal_pre_long <- vst_complete %>% 
  filter(str_detect(id, "RA_pre|normal_tissue")) %>%  
  mutate(response = if_else(str_detect(id, "RA"), 1, 0)) %>% 
  select(-id) %>% 
  pivot_longer(cols = -response, 
               names_to = "gene", 
               values_to = "log_expr")


# Nest data:

nested <- normal_pre_long %>% 
  group_by(gene) %>% 
  nest() %>% 
  ungroup()


# Model data --------------------------------------------------------------

# Model logistic regression for all genes:

mdl_data <- nested %>% mutate(mdl = map(data, 
                                        ~glm(response ~ log_expr, 
                                             data = ., 
                                             family = binomial(link = "logit"))))

# Extract model parameters:

mdl_tidy <- mdl_data %>% 
  mutate(tidy = map(mdl, tidy, conf.int = TRUE)) %>% 
  unnest(tidy)

# Filter only estimates for log coefficients:

mdl_coef <- mdl_tidy %>% filter(term == "log_expr")

mdl_coef <- read_tsv(file = "data/11_log_reg_coef.tsv")

# Add variable for significance:

mdl_sign <- mdl_coef %>% 
  mutate(significant = if_else(p.value < (0.05), 1, 0))

mdl_plot <- mdl_sign %>% mutate(neg_log10 = -log10(p.value))

# Visualise data ----------------------------------------------------------
log_reg <- mdl_plot %>% 
  ggplot(mapping = aes(x =gene,
                       y = estimate,
                       color = factor(significant))) + 
  geom_point(alpha = 0.4) +
  scale_color_discrete(na.translate = F, label = c("False","True")) +
  theme_minimal(base_family = "Avenir",
                base_size = 10) +
  labs(x = "Gene", 
       y = "Regression coefficient",
       color = "Adjusted p < 0.05") +
  theme(legend.position = "bottom",
        axis.text.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())


# Write data --------------------------------------------------------------

# Write model coefficients:
mdl_coef %>% 
  select(-c(data, mdl)) %>% 
  write_tsv(file = "./data/11_log_reg_coef.tsv")

ggsave(plot = log_reg, 
       filename = "results/11_log_reg.png", 
       width = 10, 
       height = 6, 
       units = "in")
