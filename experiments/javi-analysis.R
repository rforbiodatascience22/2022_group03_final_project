# Info --------------------------------------------------------------------
# heat plot 1: with genes significantly differential expressed for RA established
#             samples, heatmap for normal, baseline and post treatment samples. 


# Load libraries ----------------------------------------------------------
library("tidyverse")

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data ---------------------------------------------------------------

logfold_treatment <- read_tsv(file = "data/03_treatment_log2fc.tsv")
combined_meta <- read_tsv(file = "data/02_combined_meta.tsv")
combined_vst <- read_tsv(file = "data/02_combined_vst.tsv.gz")
combined_dataset <- read_tsv(file = "data/02_combined.tsv")

# Wrangle data ------------------------------------------------------------
normalized_dataset <- combined_vst %>% 
  filter(grepl("normal", id) | grepl("RA_", id)) %>% 
  pivot_longer(cols = -c(id),
               names_to = "Genes" , 
               values_to = "Reads") %>% 
  pivot_wider(names_from = id, values_from = Reads)

normal_mean <- normalized_dataset %>% 
  select(Genes, starts_with("normal")) %>% 
  mutate(mean_reads = rowMeans(across(where(is.numeric)))) %>% 
  select(Genes, mean_reads)

normalized_dataset <- normalized_dataset %>% 
  left_join(normal_mean, by = "Genes")

# log fold change (FC = sample/Normal_mean)
only_logFCs <- normalized_dataset %>% 
  transmute(
    across(
      !c(Genes, mean_reads),
      ~ log2(. + 1) - log2(mean_reads + 1)
    )
  )

combined_logfc <- normalized_dataset %>% 
  select(Genes) %>% 
  bind_cols(only_logFCs) %>% 
  pivot_longer(cols = -c(Genes), names_to = "id", values_to = "logfc") %>% 
  pivot_wider(names_from = "Genes", values_from = "logfc") %>% 
  left_join(combined_meta, by = "id")

# adding meta data for analysis
combined_logfc_meta <- combined_logfc %>% 
  filter(is.na(disease) 
         | disease == "Normal" 
         | disease == "Rheumatoid arthritis (established)") %>% 
  mutate(tag = case_when(
    treatment == TRUE ~ "RA post tDMARD",
    treatment == FALSE ~ "RA baseline",
    disease == "Normal" ~ "Normal",
    disease == "Rheumatoid arthritis (established)" ~ "RA established"
  )) %>% 
  select(-treatment, -disease_duration, -disease, -acc_num) %>% 
  relocate(tag, .after = id) %>% 
  relocate(sex, .after = tag) %>% 
  relocate(age, .after = sex)

# selecting genes significantly differential expressed from RA established 
selected_genes <- combined_logfc_meta %>% 
  filter(tag == "RA established") %>% 
  select(-age, -sex, -tag) %>% 
  pivot_longer(cols = -c(id),
               names_to = "Genes" , 
               values_to = "logfc") %>% 
  pivot_wider(names_from = id, values_from = logfc) %>% 
  mutate(RA_mean = rowMeans(across(where(is.numeric)))) %>% 
  relocate(RA_mean, .after = Genes) %>% 
  filter(RA_mean > 0.55 | RA_mean < -0.55) %>% 
  pull(Genes)


combined_logfc_selected <- combined_logfc_meta %>% 
  select(id, tag, any_of(selected_genes))


# heatmmap selected genes (normal, baseline, post treatment)
heatmap1 <- combined_logfc_selected %>% 
  filter(tag == "Normal" | tag == "RA baseline" | tag == "RA post tDMARD") %>% 
  pivot_longer(cols = -c(id, tag),
               names_to = "Genes" , 
               values_to = "logfc") %>% 
  ggplot(mapping = aes(x = id,
                       y = Genes,
                       fill = logfc )) +
  geom_tile(alpha = 0.5) +
  facet_wrap(~ tag, scales = 'free_x', ncol =) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))

#heatmap selected genes (normal, RA established)
heatmap2 <- combined_logfc_selected %>% 
  filter(tag == "Normal" | tag == "RA established") %>% 
  pivot_longer(cols = -c(id, tag),
               names_to = "Genes" , 
               values_to = "logfc") %>% 
  ggplot(mapping = aes(x = id,
                       y = Genes,
                       fill = logfc )) +
  geom_tile(alpha = 0.5) +
  facet_wrap(~ tag, scales = 'free_x', ncol =) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))

# Write data --------------------------------------------------------------
ggsave("heatmap_fig2B.png",
       heatmap1,
       path = "/cloud/project/results")

ggsave("heatmap2.png",
       heatmap2,
       path = "/cloud/project/results")