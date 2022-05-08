# Loading datasets
data8 <- read_tsv("/cloud/project/data/02_large_w_meta_clean.tsv")
drug <- read_tsv("/cloud/project/data/02_treatment_w_meta_clean.tsv")

# Data wrangling
data8 %>% 
  drop_na(age) %>% 
  filter(disease=="Normal") %>% 
  select(age) %>% 
  mean("age") 

druggood <- drug %>%
  as_tibble() %>% 
  select(1:4) %>% 
  mutate_if(is.logical, as.character) %>% 
  filter(treatment == "TRUE") %>%
  dplyr::rename(disease = treatment) %>% 
  mutate(disease = case_when(disease == 'TRUE' ~ 'Treatment cohort')) %>% 
  mutate(sex = case_when(sex == 'F' ~ 'Female')) %>%
  mutate(sex = case_when(sex == 'M' ~ 'Male')) %>%
  bind_rows(data8 %>% 
              select(1:4) 
  )

# New tibble with age means
means <- druggood %>% 
  drop_na(age) %>% 
  group_by(disease) %>% 
  summarize(Mean = round(mean(age, na.rm=TRUE))) %>% 
  as_tibble()

# Plots
density <- druggood %>% 
  drop_na(age) %>% 
  ggplot(aes(age, fill=disease)) +
  scale_fill_viridis_d(alpha = 0.5) +
  geom_density(alpha=0.4) 

histogram <- druggood %>% 
  drop_na(sex) %>% 
  ggplot(aes(disease, fill=sex)) +
  scale_fill_viridis_d(alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x.top = element_blank()) +
  geom_bar(alpha=0.4)

sex_distribution_plot <- histogram +
  facet_wrap(vars(sex)) +
  scale_x_discrete(guide = guide_axis(angle = 90))+
  xlab("Patient group") +
  ylab("Number of patients")+
  scale_fill_viridis_d(alpha = 0.5) +
  theme_minimal()

age_distribution_plot <- density + 
  geom_vline(data = means, aes(xintercept=Mean), size=0.5, color="red") +
  geom_label(data=means, size=2.5, 
             aes(x=Mean, y=0.015, label=paste0(Mean)),
             show.legend = FALSE) +
  ylab("Density")+
  xlab("Age of patients")+
  guides(fill=guide_legend(title="Condition"))+
  facet_wrap(vars(disease), ncol=1) +
  scale_fill_viridis_d(alpha = 0.5) +
  theme_minimal()

age_distribution_plot
sex_distribution_plot

# Saving plots
ggsave('/cloud/project/results/04_age_distribution.png', 
       age_distribution_plot, 
       width = 8, 
       height = 8, 
       bg = "transparent")

ggsave('/cloud/project/results/04_sex_distribution.png', 
       sex_distribution_plot, 
       width = 8, 
       height = 8, 
       bg = "transparent")

