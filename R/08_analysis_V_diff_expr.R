# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
#source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
all_logfold <- read_tsv(file = "data/03_combined_log2fc.tsv")


# Wrangle data ------------------------------------------------------------

# Find top 100 over and under expressed genes between healthy and pre
# treatment samples:

# Extract lists of genes:

over_expr <- all_logfold %>% 
  filter(differentiation == "Healthy vs. Pre treatment",
         significant) %>% 
  arrange(desc(log2FoldChange))
  select(gene) %>% 
  slice_head(n = 100) %>% 
  mutate("diff_expr" = "over")

under_expr <- all_logfold %>% 
  filter(differentiation == "Healthy vs. Pre treatment",
         significant) %>%
  arrange(desc(log2FoldChange)) %>%
  select(gene) %>% 
  slice_tail(n = 100) %>% 
  mutate("diff_expr" = "under")

# Combine lists:

diff_expr_genes <- over_expr %>% bind_rows(under_expr)

# Visualise data ----------------------------------------------------------
my_data_clean_aug %>% ...


# Write data --------------------------------------------------------------

# List of differentially expressed genes:
write_tsv(x = diff_expr_genes,
          file = "data/08_diff_expr_genes.tsv")
ggsave(...)