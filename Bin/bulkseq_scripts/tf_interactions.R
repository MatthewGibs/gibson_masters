##### Identify protein-protein interactions of biologically-relevant transcription factors. #####
collate_tf_DEGs <- function(mainDir, # Path to the 'Work' directory
                            TF_names, # List of transcription factors associated with the GO term of interest
                            GO_term # Name of the GO term of interest (for example: regulation_of_inflammatory_response)
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
  
  TFlist <- list() # Initialise the list object which will store TF to DEG matches
  
  for (i in TF_names) {
    print(sprintf("Identifying DEGs associated with %s", i))
    
    interactions <- read_tsv(sprintf("String/%s/%s_string_interactions.tsv", GO_term, i)) # Read in the .tsv file for the transcription factor
    interactions <- interactions[interactions$`#node1` == sprintf("%s", i), c(1,2)] # Subset to retain only gene names
    colnames(interactions) <- c("TF", "Gene_name") # Set the column names
    interactions[nrow(interactions)+1, ] <- list(sprintf("%s", i), sprintf("%s", i)) # Append the transcription factor itself to the end of the dataframe
    
    # Perform left joins to apply info from the DEG lists
    GMCSF <- left_join(interactions, GMCSF_DEGs, by = "Gene_name")
    MCSF <- left_join(interactions, MCSF_DEGs, by = "Gene_name")
    
    combined <- data.frame(TF = GMCSF$TF, 
                           Gene_name = GMCSF$Gene_name,
                           GMCSF_L2FC = GMCSF$log2FoldChange, 
                           MCSF_L2FC = MCSF$log2FoldChange)
    
    # Remove genes which are not DE for both treatments
    combined <- combined[!with(combined,is.na(GMCSF_L2FC) & is.na(MCSF_L2FC)),]
    
    # Generate combined objects which will be added to the list object
    GMCSF_combined <- data.frame(Gene_name = combined$Gene_name,
                                 Treatment = c(rep("GMCSF", nrow(combined))),
                                 L2FC = combined$GMCSF_L2FC)
    
    MCSF_combined <- data.frame(Gene_name = combined$Gene_name,
                                Treatment = c(rep("MCSF", nrow(combined))),
                                L2FC = combined$MCSF_L2FC)
    
    combined2 <- rbind(GMCSF_combined, MCSF_combined)
    
    TFlist[["Network"]] <- rbind(TFlist[["Network"]], combined2)
  }
  
  # Produce a file for cytoscape network construction
  cytoscape_data <- TFlist[["Network"]]
  
  cytoscape_data <- unique(cytoscape_data)
  
  Full_TFs <- data.frame(Gene_name = TF_names, Status = c(rep(1, length(TF_names))))
  
  cytoscape_data <- left_join(cytoscape_data, Full_TFs, by = "Gene_name")
  
  write.csv(cytoscape_data, sprintf("String/%s_TF_data.csv", GO_term))
  
}
