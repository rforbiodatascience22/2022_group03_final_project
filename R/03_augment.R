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
combined <- read_tsv(file = "./data/02_combined.tsv")

# Wrangle data ------------------------------------------------------------

# Add conditions column:

combined <- combined %>% 
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


# Write data --------------------------------------------------------------

write_tsv(x = all_results, file = "data/03_combined_log2fc.tsv")
