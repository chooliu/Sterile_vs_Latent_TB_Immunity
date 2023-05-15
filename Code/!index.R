# ==============================================================================
# !index.R
# master script list / table of contents
# ==============================================================================




# version control / export date to avoid overwriting
current_date_label <- "20220925"
dir.create(paste0("./FiguresTables/", current_date_label))

# basic R environment setup
source("Code/00_load_fxns_and_libraries.R")

# load metadata, activation markers, EDA
source("Code/01_load_metadata_and_markers.R")
source("Code/02_examine_distributions.R")
source("Code/03_multivariate_viz.R")
source("Code/04_table1.R")

# statistical testing - compare lineage/functional %s by group 
source("Code/05_unstimulated_comparisons.R")
source("Code/06_stimulated_comparisons.R") 
source("Code/07_test_phenograph_clusters.R") 
source("Code/08_phenograph_manual_bcg_cluster.R")

# misc. multivariate 
source("Code/09_heatmaps.R") # visualize between-marker correlation, frequencies
source("Code/10_SPICE_fdr_correct.R") # FDR correction of SPICE values

# figures for paper
source("Code/11_OR_forest_plots.R") # plots of ln(OR) & 95% confints
source("Code/12_targeted_marker_viz.R") # scatter plots for paper

# compile final results in Excel file
source("Code/13_export_results.R")


