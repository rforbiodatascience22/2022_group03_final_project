# Info --------------------------------------------------------------------
#     add this to augment later. Here i'm creating a table for the heatmaps

# Load libraries ----------------------------------------------------------
library("tidyverse")
library("forcats")
#library("DESeq2)

# Load data ---------------------------------------------------------------

combined_meta <- read_tsv(file = "data/02_combined_meta.tsv")
combined_vst <- read_tsv(file = "data/02_combined_vst.tsv.gz")
combined_dataset <- read_tsv(file = "data/02_combined.tsv")

# -------------------------------------------------------------------------
# select the data needed
# -------------------------------------------------------------------------

normalized_dataset <- combined_vst %>% 
  filter(grepl("normal", id) | grepl("RA_", id)) %>% 
  pivot_longer(cols = -c(id),
               names_to = "Genes" , 
               values_to = "Reads") %>% 
  pivot_wider(names_from = id, values_from = Reads)

# -------------------------------------------------------------------------
# Mean of normal samples expression
# -------------------------------------------------------------------------

normal_mean <- normalized_dataset %>% 
  select(Genes, starts_with("normal")) %>% 
  mutate(mean_reads = rowMeans(across(where(is.numeric)))) %>% 
  select(Genes, mean_reads)

normalized_dataset <- normalized_dataset %>% 
  left_join(normal_mean, by = "Genes")

# -------------------------------------------------------------------------
# log fold change (FC = sample/Normal_mean)
# -------------------------------------------------------------------------

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

# -------------------------------------------------------------------------
# Wrangle data
# -------------------------------------------------------------------------

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


# Write data --------------------------------------------------------------

write_tsv(x = combined_logfc_meta, file = "data/03_heatmap_log2fc.tsv")
