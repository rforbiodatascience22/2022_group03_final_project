# Load libraries ----------------------------------------------------------
library("tidyverse")
library("ggplot2")

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
large_w_meta <- read_tsv(file = "./data/02_large_w_meta_clean.tsv")
treatment_w_meta <- read_tsv(file = "./data/02_treatment_w_meta_clean.tsv")


# Wrangle data ------------------------------------------------------------
long_model <- large_w_meta %>% 
  select(!c(sex, age, acc_num)) %>% 
  filter(disease == c("Normal","Rheumatoid arthritis (early)")) %>% 
  mutate(outcome = case_when(disease == "Normal" ~ 0,
                             disease == "Rheumatoid arthritis (early)" ~ 1)) %>% 
  relocate(outcome)

test_data<- treatment_w_meta %>% 
  filter(treatment == "FALSE") %>% 
  select(!c(treatment, sex, age, disease_duration, acc_num))
  
  

# Model data
earlyRA.model <- long_model %>%
  glm(outcome ~ LXN + IL8,
      data = .,
      family = binomial(link = "logit"))


  

# summary(earlyRA.model)

### test the model
test_data$y_pred = predict(earlyRA.model, test_data, type="response")


# Visualise data ----------------------------------------------------------
LXN_curve <- long_model %>% 
  ggplot(mapping =  aes(x=LXN, y=outcome)) + 
  geom_point(alpha=.5) +
  stat_smooth(method="glm",
              se=FALSE,
              method.args = list(family=binomial)
  )

IL8_curve <- long_model %>% 
  ggplot(mapping =  aes(x=IL8, y=outcome)) + 
  geom_point(alpha=.5) +
  stat_smooth(method="glm",
              se=FALSE,
              method.args = list(family=binomial)
  )


LXN_curve
IL8_curve



# Write data --------------------------------------------------------------
ggsave("model1_LXN_curve.png",
       LXN_curve,
       path = "/cloud/project/results")

ggsave("model1_IL8_curve.png",
       IL8_curve,
       path = "/cloud/project/results")

