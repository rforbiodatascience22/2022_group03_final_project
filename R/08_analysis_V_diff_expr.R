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
  arrange(desc(log2FoldChange)) %>% 
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

all_logfold %>% 
  left_join(diff_expr_genes, by = "gene") %>%
  mutate(differentiation = fct_relevel(differentiation, 
                                       c("Healthy vs. Pre treatment",
                                         "Healthy vs. Post Treatment",
                                         "Pre vs. Post Treatment"
                                       ))) %>% 
  ggplot(mapping = aes(x = gene, 
                       y = log2FoldChange,
                       color = factor(significant),
                       fill = factor(diff_expr))) + 
  geom_point(alpha = 0.4, shape = 21) +
  scale_color_manual(values = c("black", "dodgerblue3"), na.translate = F) +
  scale_fill_manual(values = c("under" = "brown1", "over" = "chartreuse"), na.translate = F) +
  theme_minimal(base_family = "Avenir",
                base_size = 12) +
  labs(x = "Gene", 
       y = "log2 Fold Change",
       fill = "Differential expression",
       color = "Adjusted p < 0.05") +
  lims(y = c(-15,15)) +
  theme(legend.position = "bottom",
        axis.text.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  facet_grid(rows = vars(differentiation))


# Write data --------------------------------------------------------------

# List of differentially expressed genes:
write_tsv(x = diff_expr_genes,
          file = "data/08_diff_expr_genes.tsv")
ggsave(...)