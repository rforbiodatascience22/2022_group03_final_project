# Run all scripts ---------------------------------------------------------
source(file = "R/01_load.R")
source(file = "R/02_clean.R")
source(file = "R/03_augment.R")
source(file = "R/04_analysis_exploration.R")
source(file = "R/05_analysis_pca.R")
source(file = "R/06_analysis_diff_expr.R")
source(file = "R/07_analysis_heatmap.R")
source(file = "R/08_analysis_boxpl_nrml.R")
source(file = "R/09_analysis_model.R")
#source(file = "R/11_analysis_log_reg.R") # Optional and VERY slow.
source(file = "R/14_analysis_eda.R")
rmarkdown::render("./doc/slides.Rmd", "ioslides_presentation")
