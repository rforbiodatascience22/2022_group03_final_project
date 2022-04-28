# Load libraries ----------------------------------------------------------
library("tidyverse")

# Load dataset
pca_fit_male <- read_tsv(file = "data/04_pca_fit_male.tsv")

kclust <- kmeans(points, centers = 3)
kclust