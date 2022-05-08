# Load libraries ----------------------------------------------------------
library("tidyverse")
library("ggplot2")

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
vst_data <- read_tsv(file = "data/02_combined_vst.tsv.gz")


# Wrangle data ------------------------------------------------------------
vst_data_long <-
  vst_data %>% 
  select(id,
         CD3D,
         MS4A1,
         CTLA4,
         CD19,
         IL10,
         CLEC12A,
         AURKA,
         abParts,
         CLEC2B,
         CD58) %>%
  separate(col = id, 
           into = c("a", "b"),
           sep = "_",
           remove = FALSE) %>% 
  unite("condition", a:b, sep= " ") %>% 
  mutate( condition = case_when(
    condition == "RA pre" ~ "RA baseline",
    condition == "RA post" ~ "RA post-tDMARD",
    condition == "normal tissue" ~ "Normal",
    condition == "OA tissue" ~ "Osteoarthritis",
    condition == "AG tissue" ~ "Arthralgia",
    condition == "undiff tissue" ~ "Undifferentiated",
    condition == "RA tissue" ~ "RA early & late",
  )
  ) %>% 
  pivot_longer(cols = -c(id, condition),
               names_to = "gene" , 
               values_to = "expression_level")


# Visualise data ----------------------------------------------------------

# Boxplot with all the conditions (both long and treatment data)
allcond_sele_genes <-
  vst_data_long %>% 
  ggplot(mapping = aes(
    x = condition,
    y = expression_level
  )
  ) +
  geom_boxplot(outlier.size = .5) +
  geom_jitter(aes(colour = condition),
              size = 0.7
  )+
  facet_wrap(~factor(gene, levels =c("CD3D",
                                     "MS4A1",
                                     "CTLA4",
                                     "CD19",
                                     "abParts",
                                     "CLEC12A",
                                     "AURKA",
                                     "IL10",
                                     "CLEC2B",
                                     "CD58")
                     ),
             scales = 'free_y', nrow = 2) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank()
        ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )+
  ylab("Expression Level") +
  labs(col="")

# Boxplot of only treatment and normal data
drug_sele_genes <-
  vst_data_long %>% 
  filter(condition == "RA baseline" 
         | condition == "RA post-tDMARD"
         | condition == "Normal"
  ) %>% 
  ggplot(mapping = aes(
    x = condition,
    y = expression_level
  )
  ) +
  geom_boxplot(outlier.size = .5) +
  geom_jitter(aes(fill = condition),
              colour = "black",
              pch = 21, 
              size = 2
  )+
  facet_wrap(~factor(gene, levels =c("CD3D",
                                     "MS4A1",
                                     "CTLA4",
                                     "CD19",
                                     "abParts",
                                     "CLEC12A",
                                     "AURKA",
                                     "IL10",
                                     "CLEC2B",
                                     "CD58")
                     ), 
                     scales = 'free_y',
                     nrow = 2) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank()
        ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
        )+
  ylab("Expression Level") +
  scale_fill_viridis_d(alpha = 0.5) +
  labs(col="")
  
# Write data --------------------------------------------------------------
ggsave("08_boxp_sele_genes_allcond.png",
       plot = allcond_sele_genes,
       bg = "transparent",
       dpi = 300,
       width = 12, 
       height = 10,
       units = "in",
       path = "/cloud/project/results")

ggsave("08_bxpl_sele_genes_drug.png",
       plot = drug_sele_genes,
       bg = "transparent",
       dpi = 300,
       width = 12, 
       height = 10,
       units = "in",
       path = "/cloud/project/results")
