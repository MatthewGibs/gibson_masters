##### Identify protein-protein interactions of biologically-relevant transcription factors. #####

collate_tf_DEGs <- function(mainDir # Path to the 'Work directory
                            ) {
  setwd(mainDir)
  string_packages <- c("clusterProfiler", 
                     "AnnotationDbi", 
                     "org.Hs.eg.db",
                     "R.utils", 
                     "biomaRt",
                     "dplyr",
                     "ggplot2",
                     "data.table",
                     "tidyverse")
  
  suppressMessages(lapply(string_packages, require, character.only = TRUE))
  
  # Read in and annotate the two lists of DEGs
  DEGs <- readRDS("DEG_lists/GMCSF_and_MCSF_comparison.rds")
  
  # Convert ENSG ID to gene name
  source(sprintf("%s/biomaRt_annotate.R", bin_path))
  generate_ensembl(mainDir)
  
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "external_gene_name", # The ID type you want in your final object.
             object = DEGs$GMCSF$Gene_name # character vector containing the 'initial ID' type that you want to convert
  )
  DEGs$GMCSF$Gene_name <- object_final$external_gene_name
  GMCSF_DEGs <- DEGs$GMCSF
  
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "external_gene_name", # The ID type you want in your final object.
             object = DEGs$MCSF$Gene_name # character vector containing the 'initial ID' type that you want to convert
  )
  DEGs$MCSF$Gene_name <- object_final$external_gene_name
  MCSF_DEGs <- DEGs$MCSF
  
  TFlist <- list()
  TFs_to_DEGs <- function(TF_names, GO_term, folder_name) {
    for (i in TF_names) {
      
      print(sprintf("Identifying DEGs associated with %s", i))
      
      interactions <- read_tsv(sprintf("String/%s/%s/%s_string_interactions.tsv", GO_term, folder_name, i))
      interactions <- interactions[interactions$`#node1` == sprintf("%s", i), c(1,2,9)]
      colnames(interactions) <- c("TF", "Gene_name", "Coexpression")
      interactions[nrow(interactions)+1, ] <- list(sprintf("%s", i), sprintf("%s", i), 1)
      
      GMCSF <- left_join(interactions, GMCSF_DEGs, by = "Gene_name")
      MCSF <- left_join(interactions, MCSF_DEGs, by = "Gene_name")
      
      combined <- data.frame(TF = GMCSF$TF, 
                             Gene_name = GMCSF$Gene_name, 
                             Coexpression = GMCSF$Coexpression,
                             GMCSF_L2FC = GMCSF$log2FoldChange, 
                             MCSF_L2FC = MCSF$log2FoldChange)
      
      # Remove genes which are not DE for both treatments
      combined <- combined[!with(combined,is.na(GMCSF_L2FC) & is.na(MCSF_L2FC)),]
      
      GMCSF_combined <- data.frame(Gene_name = combined$Gene_name,
                                   Treatment = c(rep("GMCSF", nrow(combined))),
                                   L2FC = combined$GMCSF_L2FC)
      
      MCSF_combined <- data.frame(Gene_name = combined$Gene_name,
                                  Treatment = c(rep("MCSF", nrow(combined))),
                                  L2FC = combined$MCSF_L2FC)
      
      combined2 <- rbind(GMCSF_combined, MCSF_combined)
      
      TFlist[[folder_name]][[i]] <- combined
      TFlist[["Network"]] <- rbind(TFlist[["Network"]], combined2)
      TFlist[["Full"]] <- rbind(TFlist[["Full"]], combined)
    }
    
    assign("TFlist", TFlist, envir = .GlobalEnv)
  }
  
  # Define the set of genes for examining protein-protein interaction - regulation of inflammatory response.
  MCSF_TFs <- c("CLOCK", "RB1")
  GMCSF_TFs <- c("NFKB1", "PPARD", "NR1D2")
  
  TFs_to_DEGs(TF_names = MCSF_TFs, 
              GO_term = "Regulation_of_inflammatory_response",
              folder_name = "MCSF")
  
  TFs_to_DEGs(TF_names = GMCSF_TFs, 
              GO_term = "Regulation_of_inflammatory_response",
              folder_name = "GMCSF")
  str(TFlist)
  # Produce a file for cytoscape network construction
  cytoscape_data <- TFlist[["Network"]]
  
  cytoscape_data <- unique(cytoscape_data)
  str(cytoscape_data)
  Full_TFs <- c(MCSF_TFs,
                GMCSF_TFs)
  Full_TFs <- data.frame(Gene_name = Full_TFs, Status = c(rep(1, length(Full_TFs))))
  str(Full_TFs)
  cytoscape_data <- left_join(cytoscape_data, Full_TFs, by = "Gene_name")
  
  write.csv(cytoscape_data, "String/regualtion_of_inflammatory_response_TF_data.csv")
  
  # Define the set of genes for examining protein-protein interaction - myeloid cell differentiation.
  MCSF_TFs <- c("SPI1", "RUNX1")
  GMCSF_TFs <- c("MYC", "MAFB")
  
  TFlist <- list()
  
  TFs_to_DEGs(TF_names = MCSF_TFs, 
              GO_term = "Myeloid_cell_differentiation",
              folder_name = "MCSF")
  
  TFs_to_DEGs(TF_names = GMCSF_TFs, 
              GO_term = "Myeloid_cell_differentiation",
              folder_name = "GMCSF")
  
  # Produce a file for cytoscape network construction
  cytoscape_data <- TFlist[["Network"]]
  
  cytoscape_data <- unique(cytoscape_data)
  
  Full_TFs <- c(MCSF_TFs, GMCSF_TFs)
  Full_TFs <- data.frame(Gene_name = Full_TFs, Status = c(rep(1, length(Full_TFs))))
  
  cytoscape_data <- left_join(cytoscape_data, Full_TFs, by = "Gene_name")
  
  write.csv(cytoscape_data, "String/myeloid_cell_differentiation_TF_data.csv")
}
