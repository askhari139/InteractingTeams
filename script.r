library(funcsKishore)
pwd <- "/Users/kishorehari/Desktop/PostDoc/MultiAxial/Networks/MiDAS2/InteractingTeams"
setwd(pwd)
figDir <- "figures"
moduleTSdata <- read_csv("Data_export_for_sharing/panel1_data/panel1_data_tidy.csv")

##outputs

pca <- read_csv("Data_export_for_sharing/pca_metrics/pca_metrics_by_folder.csv")
