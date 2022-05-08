# Info --------------------------------------------------------------------
# heat plot 1: with genes significantly differential expressed for RA established
#             samples, heatmap for normal, baseline and post treatment samples. 


# Load libraries ----------------------------------------------------------
library("tidyverse")

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data ---------------------------------------------------------------

heatmap_table <- read_tsv(file = "data/03_heatmap_log2fc.tsv")
heatmap_genes <- read_tsv(file = "data/08_diff_expr_genes.tsv")

# getting selected genes 

selected_genes_under <- heatmap_genes %>% 
  filter(diff_expr == "under") %>% 
  pull(gene)

selected_genes_over <- heatmap_genes %>% 
  filter(diff_expr == "over") %>% 
  pull(gene)


# heatmmap selected genes (normal, baseline, post treatment)
heatmap_under <- heatmap_table %>%
  select(id, tag, any_of(selected_genes_under)) %>% 
  filter(tag == "Normal" | tag == "RA baseline" | tag == "RA post tDMARD") %>% 
  pivot_longer(cols = -c(id, tag),
               names_to = "Genes" , 
               values_to = "logfc") %>% 
  ggplot(mapping = aes(x = id,
                       y = Genes,
                       fill = logfc )) +
  geom_tile(alpha = 0.5) +
  facet_wrap(~ tag, scales = 'free_x', strip.position = "bottom") +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 16)) +
  labs(title = "Genes significantly under regulated in Rheumatoid Arthritis")

heatmap_over <- heatmap_table %>%
  select(id, tag, any_of(selected_genes_over)) %>% 
  filter(tag == "Normal" | tag == "RA baseline" | tag == "RA post tDMARD") %>% 
  pivot_longer(cols = -c(id, tag),
               names_to = "Genes" , 
               values_to = "logfc") %>% 
  ggplot(mapping = aes(x = id,
                       y = Genes,
                       fill = logfc )) +
  geom_tile(alpha = 0.5) +
  facet_wrap(~ tag, scales = 'free_x', strip.position = "bottom") +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 16)) +
  labs(title = "Genes significantly over regulated in Rheumatoid Arthritis")


heatmap_over
heatmap_under
# Write data --------------------------------------------------------------
ggsave("heatmap_under.png",
       heatmap_under,
       dpi = 300,
       width = 8, 
       height = 10,
       units = "in",
       path = "/cloud/project/results")

ggsave("heatmap_over.png",
       heatmap_over,
       dpi = 300,
       width = 8, 
       height = 10,
       units = "in",
       path = "/cloud/project/results")
