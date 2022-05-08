# Info --------------------------------------------------------------------
#     Calculates log2 fold changes between conditions using DESeq2.

# Load libraries ----------------------------------------------------------
library("tidyverse")
library("forcats")
#library("DESeq2)

# Define functions --------------------------------------------------------

# Only used here:

calc_logfold <- function(data) {
  
  # Returns log2 fold changes from data frame
  
  # Filter and convert counts to integer
  count <- data %>%
    select(-condition) %>% 
    pivot_longer(cols = -id, names_to = "gene") %>% 
    pivot_wider(names_from = id, values_from = value) %>% 
    rowwise() %>% 
    filter(sum(c_across(-gene)) > 200) %>% 
    ungroup() %>% 
    mutate(across(where(is.numeric), ~ as.integer(.x))) %>% 
    as.data.frame()
  
  # Extract meta data for DESeq
  meta <- data %>% 
    select(c(id,condition)) %>% 
    as.data.frame()
  
  # Convert to DESeq dataset
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count,
                                        colData = meta,
                                        design = ~condition,
                                        tidy = TRUE)
  # Run DESeq
  dds <- DESeq2::DESeq(dds)
  
  # Extract and return results as tibble:
  
  DESeq2::results(dds) %>% as_tibble(rownames = "gene")
  }


# Load data ---------------------------------------------------------------
combined_dataset <- read_tsv(file = "data/02_combined.tsv")
combined_meta <- read_tsv(file = "data/02_combined_meta.tsv")
combined_vst <- read_tsv(file = "data/02_combined_vst.tsv.gz")


# Wrangle data ------------------------------------------------------------

# Add conditions column:

combined <- combined_dataset %>% 
  mutate("condition" = case_when(str_detect(id,"pre") ~ "pre",
                                 str_detect(id,"post") ~ "post",
                                 str_detect(id,"normal") ~ "normal"))


# -------------------------------------------------------------------------
# Log2 Fold Change: Normal vs. RA_pre
# -------------------------------------------------------------------------

normal_vs_pre <- combined %>% 
  filter(str_detect(id, "RA_pre|normal_tissue")) %>% 
  arrange(id)

results_normal_vs_pre <- normal_vs_pre %>% 
  calc_logfold() %>% 
  filter(!is.na(padj)) %>% 
  mutate(significant = if_else(padj < 0.05, TRUE, FALSE),
         "differentiation" = "Healthy vs. Pre treatment")

# -------------------------------------------------------------------------
# Log2 Fold Change: RA_pre vs. RA_post
# -------------------------------------------------------------------------

pre_vs_post <- combined %>% 
  filter(str_detect(id, "RA_p")) %>% 
  arrange(desc(id))

results_pre_vs_post <- pre_vs_post %>% 
  calc_logfold() %>% 
  filter(!is.na(padj)) %>% 
  mutate(significant = if_else(padj < 0.05, TRUE, FALSE),
         "differentiation" = "Pre vs. Post Treatment")

# -------------------------------------------------------------------------
# Log2 Fold Change: Normal vs. RA_post
# -------------------------------------------------------------------------

normal_vs_post <- combined %>% 
  filter(str_detect(id, "RA_post|normal")) %>% 
  arrange(desc(id))

results_normal_vs_post <- normal_vs_post %>% 
  calc_logfold() %>% 
  filter(!is.na(padj)) %>% 
  mutate(significant = if_else(padj < 0.05, TRUE, FALSE),
         "differentiation" = "Healthy vs. Post Treatment")

# -------------------------------------------------------------------------
# Combine results to one data frame
# -------------------------------------------------------------------------

all_results <- results_normal_vs_pre %>% 
  bind_rows(results_pre_vs_post) %>%
  bind_rows(results_normal_vs_post)



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

write_tsv(x = all_results, file = "data/03_combined_log2fc.tsv")
write_tsv(x = combined_logfc_meta, file = "data/03_heatmap_log2fc.tsv")
