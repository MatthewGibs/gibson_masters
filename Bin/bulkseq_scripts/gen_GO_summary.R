##### Create and write to file an essential summary of all CP results #####
summarise_GO <- function(mainDir) {
  
  ##### Setup #####
  setwd(mainDir)
  
  clusterProfiler_plots_packages <- c("data.table", "clusterProfiler", "ggplot2", "tidyverse", "dplyr", "org.Hs.eg.db", "enrichplot", "forcats")
  
  suppressMessages(lapply(clusterProfiler_plots_packages, require, character.only = TRUE))
  
  # Create the output directory
  output_dir <- "clusterProfiler/GO_enrichment"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Creating 'clusterProfiler_summary' directory.")
  } else {}
  
  ##### Generate useful list summarising all GO results #####
  # Read in all GO enrichment results
  GMCSF_both <- read_csv("clusterProfiler/GO_enrichment/GMCSF/GO_BOTH_significant.csv")
  GMCSF_down <- read_csv("clusterProfiler/GO_enrichment/GMCSF/GO_DOWN_significant.csv")
  GMCSF_up <- read_csv("clusterProfiler/GO_enrichment/GMCSF/GO_UP_significant.csv")
  
  MCSF_both <- read_csv("clusterProfiler/GO_enrichment/MCSF/GO_BOTH_significant.csv")
  MCSF_down <- read_csv("clusterProfiler/GO_enrichment/MCSF/GO_DOWN_significant.csv")
  MCSF_up <- read_csv("clusterProfiler/GO_enrichment/MCSF/GO_UP_significant.csv")
  
  Common_both <- read_csv("clusterProfiler/GO_enrichment/GMCSF_MCSF_common/GO_BOTH_significant.csv")
  Common_up <- read_csv("clusterProfiler/GO_enrichment/GMCSF_MCSF_common/GO_UP_significant.csv")
  Common_down <- read_csv("clusterProfiler/GO_enrichment/GMCSF_MCSF_common/GO_DOWN_significant.csv")
  
  GMCSF_unique_both <- read_csv("clusterProfiler/GO_enrichment/GMCSF_unique/GO_BOTH_significant.csv")
  GMCSF_unique_up <- read_csv("clusterProfiler/GO_enrichment/GMCSF_unique/GO_UP_significant.csv")
  GMCSF_unique_down <- read_csv("clusterProfiler/GO_enrichment/GMCSF_unique/GO_DOWN_significant.csv")
  
  MCSF_unique_both <- read_csv("clusterProfiler/GO_enrichment/MCSF_unique/GO_BOTH_significant.csv")
  MCSF_unique_up <- read_csv("clusterProfiler/GO_enrichment/MCSF_unique/GO_UP_significant.csv")
  MCSF_unique_down <- read_csv("clusterProfiler/GO_enrichment/MCSF_unique/GO_DOWN_significant.csv")
  
  GMCSF <- list(GMCSF_both, GMCSF_up, GMCSF_down)
  MCSF <- list(MCSF_both, MCSF_up, MCSF_down)
  Common <- list(Common_both, Common_up, Common_down)
  GMCSF_unique <- list(GMCSF_unique_both, GMCSF_unique_up, GMCSF_unique_down)
  MCSF_unique <- list(MCSF_unique_both, MCSF_unique_up, MCSF_unique_down)
  names(GMCSF) <- c("GMCSF_both", "GMCSF_up", "GMCSF_down")
  names(MCSF) <- c("MCSF_both", "MCSF_up", "MCSF_down")
  names(Common) <- c("Common_both", "Common_up", "Common_down")
  names(GMCSF_unique) <- c("GMCSF_unique_both", "GMCSF_unique_up", "GMCSF_unique_down")
  names(MCSF_unique) <- c("MCSF_unique_both", "MCSF_unique_up", "MCSF_unique_down")
  
  # finalize the list
  GO_results <- list(GMCSF, MCSF, Common, GMCSF_unique, MCSF_unique)
  names(GO_results) <- c("GMCSF", "MCSF", "Common", "GMCSF_unique", "MCSF_unique")
  
  rm(GMCSF, MCSF, Common, GMCSF_unique, MCSF_unique, 
     GMCSF_both, GMCSF_up, GMCSF_down,
     MCSF_both, MCSF_up, MCSF_down,
     Common_both, Common_up, Common_down,
     GMCSF_unique_both, GMCSF_unique_up, GMCSF_unique_down,
     MCSF_unique_both, MCSF_unique_up, MCSF_unique_down)
  
  # Save this list to file for later use
  saveRDS(GO_results, file = "clusterProfiler/GO_enrichment/GO_summary.rds")
}
