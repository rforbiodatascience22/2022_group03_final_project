# Load libraries ----------------------------------------------------------
library("tidyverse")
library("ggplot2")

# Define functions --------------------------------------------------------
#source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
all <- read_tsv(file = "./data/02_combined_vst.tsv.gz")


# Wrangle data ------------------------------------------------------------
ra_pre <- all %>%  
  filter(str_detect(id, "RA_pre_*")) %>% 
  mutate(condition = "Early RA")
ra_i <- ra_pre %>% 
  slice_sample(n = nrow(ra_pre)/2, replace = FALSE ) %>% 
  mutate(outcome = 1)
ra_ii <- ra_pre %>% 
  anti_join(ra_i, by = "id")

# early_ra_long <- all %>% 
#   filter(str_detect(id, "RA_tissue_*")) %>% 
#   slice_head(n = 57) %>%  # id < RA_tissue_58 are early RA patients
#   mutate(condition = "Early RA")
# ra_iii <- early_ra_long %>% 
#   slice_sample(n = nrow(early_ra_long)/2, replace = FALSE ) %>% 
#   mutate(outcome = 1)
# ra_iv <-  early_ra_long %>% 
#   anti_join(ra_iii, by = "id")

normal <- all %>% 
  filter(str_detect(id, "normal_tissue_*"))%>% 
  mutate(condition = "Normal")
# split normal for the training and testing
normal_i <- normal %>% 
  slice_sample(n = nrow(normal)/2, replace = FALSE ) %>% 
  mutate(outcome = 0)
#selects all rows from normal that are not present in normal_i
normal_ii <- normal %>% 
  anti_join(normal_i, by = "id")

train_data <-
  normal_i %>% 
  bind_rows(ra_i) %>% 
  # bind_rows(ra_iii) %>% 
  relocate(outcome)


test_data <-
  normal_ii %>% 
  bind_rows(ra_ii)
  # bind_rows(ra_iv)


# Model data
earlyRA.model_i <- train_data %>%
  glm(outcome ~ LXN ,
      data = .,
      family = binomial(link = "logit"))

# earlyRA.model_ii <- train_data %>%
#   glm(outcome ~ LXN + IL8 ,
#       data = .,
#       family = binomial(link = "logit"))

model_sum <- summary(earlyRA.model_i)

### test the model
test_data$outcome = predict(earlyRA.model_i, test_data, type="response")
test_data_outcome <- test_data %>% 
  relocate(outcome)



# Visualise data ----------------------------------------------------------
LXN_curve <- train_data %>% 
  ggplot(mapping =  aes(x=LXN, y=outcome)) + 
  geom_point(alpha=.5) +
  stat_smooth(method="glm",
              se=FALSE,
              method.args = list(family=binomial),
              #col="purple",
              lty=2
  ) +
  # geom_point(data = test_data_outcome,
  #            mapping = aes(x=LXN, 
  #                          y=outcome,
  #                          colour = condition,
  #                          shape = condition),
  #            size = 4,
  #            fill = "white") +
  geom_point(data = test_data_outcome,
             mapping = aes(
               x=LXN,
               y=outcome,
               colour = condition
             ), 
             size = 4) +
  geom_point(data = test_data_outcome,
             mapping = aes(
               x=LXN,
               y=outcome
               ),
             colour = "white", size = 1.5) +
  theme_minimal() +
  labs(
    title = "Simple model to distinguish early RA from normal condition",
    subtitle = "Based on RNA-seq data for the gene LXN (CXCL8)",
    caption = "Input data were normalized reads."
  )


# Write data --------------------------------------------------------------
ggsave("final_model_LXN_curve.png",
       LXN_curve,
       path = "/cloud/project/results")
