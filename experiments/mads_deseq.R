library("tidyverse")

# DESeq2
# -------------------------------------------------------------------------
# Normal vs. RA_pre
# -------------------------------------------------------------------------

combined <- read_tsv(file = "./data/02_combined.tsv") %>% 
  filter(str_detect(id, "RA_pre|normal_tissue")) %>% 
  arrange(id)

count <- combined %>% 
  pivot_longer(cols = -id, names_to = "gene") %>% 
  pivot_wider(names_from = id, values_from = value) %>% 
  rowwise() %>% 
  filter(sum(c_across(-gene)) > 200) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric), ~ as.integer(.x))) %>% 
  as.data.frame()

meta <- combined %>% 
  select(id) %>% 
  mutate(condition = if_else(str_detect(id,"normal"), 
                             "normal",
                             "pre")) %>% 
  as.data.frame()

dds <- DESeq2::DESeqDataSetFromMatrix(countData = count,
                               colData = meta,
                               design = ~condition,
                               tidy = TRUE)
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds)

results_pre_normal <- res %>% 
  as_tibble(rownames = "gene") %>% 
  filter(!is.na(padj)) %>% 
  mutate(significant = if_else(padj < 0.05, TRUE, FALSE),
         "differentiation" = "Normal vs. Early"
         )

results_pre_normal %>% 
  ggplot(mapping = aes(x = gene, 
                       y = log2FoldChange,
                       color = significant)) + 
  geom_point(alpha = 0.5) +
  theme_minimal()

# -------------------------------------------------------------------------
# Normal vs. RA_tissue
# -------------------------------------------------------------------------


combined <- read_tsv(file = "./data/02_combined.tsv") %>% 
  filter(str_detect(id, "RA_tissue|normal_tissue"))

# Subsample RA_tissue
set.seed(18)
RA_tissue <- combined %>% 
  filter(str_detect(id, "RA_tissue")) %>%
  sample_n(28)

# Combine with healthy tissue:
combined <- combined %>% 
  filter(str_detect(id, "normal")) %>% 
  bind_rows(RA_tissue)

count <- combined %>% 
  pivot_longer(cols = -id, names_to = "gene") %>% 
  pivot_wider(names_from = id, values_from = value) %>% 
  rowwise() %>% 
  filter(sum(c_across(-gene)) > 200) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric), ~ as.integer(.x))) %>% 
  as.data.frame()

meta <- combined %>% 
  select(id) %>% 
  mutate(condition = if_else(str_detect(id,"normal"), 
                             "normal",
                             "RA")) %>% 
  as.data.frame()

dds <- DESeq2::DESeqDataSetFromMatrix(countData = count,
                                      colData = meta,
                                      design = ~condition,
                                      tidy = TRUE)
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds)

results_RA_normal <- res %>% 
  as_tibble(rownames = "gene") %>% 
  filter(!is.na(padj)) %>% 
  mutate(significant = if_else(padj < 0.05, TRUE, FALSE),
         "differentiation" = "Normal vs. Late")


results_RA_normal %>% 
  ggplot(mapping = aes(x = gene, 
                       y = log2FoldChange,
                       color = significant)) + 
  geom_point(alpha = 0.5) +
  theme_minimal()

# -------------------------------------------------------------------------
# RA_pre vs. RA_post
# -------------------------------------------------------------------------


combined <- read_tsv(file = "./data/02_combined.tsv") %>% 
  filter(str_detect(id, "RA_p")) %>% 
  arrange(desc(id))


count <- combined %>% 
  pivot_longer(cols = -id, names_to = "gene") %>% 
  pivot_wider(names_from = id, values_from = value) %>% 
  rowwise() %>% 
  filter(sum(c_across(-gene)) > 200) %>% 
  ungroup() %>% 
  mutate(across(where(is.numeric), ~ as.integer(.x))) %>% 
  as.data.frame()

meta <- combined %>% 
  select(id) %>% 
  mutate(condition = if_else(str_detect(id,"pre"), 
                             "pre",
                             "post")) %>% 
  as.data.frame()

dds <- DESeq2::DESeqDataSetFromMatrix(countData = count,
                                      colData = meta,
                                      design = ~condition,
                                      tidy = TRUE)
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds)
DESeq2::summary(res)

results_pre_post <- res %>% 
  as_tibble(rownames = "gene") %>% 
  mutate(significant = if_else(padj < 0.05, TRUE, FALSE),
         "differentiation" = "Early vs. Post Treatment")


results_pre_post %>% 
  ggplot(mapping = aes(x = gene, 
                       y = log2FoldChange,
                       color = significant)) + 
  geom_point(alpha = 0.5) +
  theme_minimal()


# All results:

all_results <- results_pre_normal %>% 
  bind_rows(results_RA_normal) %>%
  bind_rows(results_pre_post)

write_tsv(x = all_results, file = "./experiments/all_deseq2_results.tsv")



# Top 100 of head and tail from:

# Pre vs. normal:
pre_normal_head <- 
  results_pre_normal %>% 
  filter(significant) %>% 
  arrange(padj) %>%
  select(gene) %>% 
  slice_head(n = 200)

pre_normal_tail <- 
  results_pre_normal %>% 
  filter(significant) %>% 
  arrange(padj) %>%
  select(gene) %>% 
  slice_tail(n = 750)

# RA vs. normal:
RA_normal_head <- 
  results_RA_normal %>% 
  filter(significant) %>% 
  arrange(padj) %>%
  select(gene) %>% 
  slice_head(n = 200)

RA_normal_tail <- 
  results_RA_normal %>% 
  filter(significant) %>% 
  arrange(padj) %>%
  select(gene) %>% 
  slice_tail(n = 750)

# Intersect between heads and tails:

head <- pre_normal_head %>% intersect(RA_normal_head)
tail <- pre_normal_tail %>% intersect(RA_normal_tail)


# Intersect with post:

results_pre_post %>% 
  filter(significant) %>% 
  select(gene) %>% 
  intersect(tail)

head

results_pre_post %>%
  filter(significant) %>% 
  semi_join(tail)

results_pre_post


# Pre vs. Post:
post_normal_head <- 
  results_pre_post %>% 
  filter(significant) %>% 
  arrange(padj) %>% 
  slice_head(n = 100)

post_normal_tail <- 
  results_pre_normal %>% 
  filter(significant) %>% 
  arrange(padj) %>% 
  slice_tail(n = 100)




all_results %>% 
  ggplot(mapping = aes(x = gene, 
                       y = log2FoldChange,
                       color = factor(significant))) + 
  geom_point(alpha = 0.3) +
  scale_color_manual(values = c("black","dodgerblue3")) +
  theme_minimal(base_family = "Avenir",
                base_size = 12) +
  labs(x = "Gene", 
       y = "log2 Fold Change",
       color = "Adjusted p < 0.05") +
  theme(legend.position = "bottom",
        axis.text.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  facet_grid(rows = vars(differentiation),
             scales = "free")
