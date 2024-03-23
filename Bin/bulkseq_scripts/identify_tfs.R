##### Identify transcription factors from supplied lists of DEGs #####
# This script contains a function to generate lists of unique and common DEGs based on two input lists of DEGs. 
# It then saves these to file.
identify_tfs <- function(mainDir, # path to the working directory
                         DEG_list, # path to a .csv containing DEGs
                         output_name # The prefix used to name the file, typically a treatment/condition type that dictated the DEGs.
) {
  # Move to the working directory
  setwd(mainDir)
  
  # Read in the biomart script
  source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
  generate_ensembl(mainDir)
  
  df <- read.csv(sprintf("%s", DEG_list), row.names = 1) #generate a df from a saved .csv file
  colnames(df) <- c("Gene_name", "log2FoldChange", "padj") # set the column names of the data frame
  rownames(df) <- df$Gene_name
  
  # Gather necessary attributes to identify transcription factors
  gather_attributes(initial_ID = "ensembl_gene_id", # This is the starting ID type contained in the 'object'.
                    target_attributes = c("name_1006"), # A vector listing all the target attributes you wish to gather information for.
                    object = df$Gene_name # character vector containing the 'initial ID' type that you want to convert
  )
  
  # Filter for genes with known DNA-binding domains
  transcription_factors <- object_final[grepl("sequence-specific DNA binding", object_final$name_1006), ]
  
  # Further filter this list to retain unique rows and essential information
  transcription_factors <- data.frame(transcription_factors[1])
  transcription_factors <- unique(transcription_factors)
  colnames(transcription_factors)[1] <- "Gene_name"
  
  # Regain information on log2 fold change.
  final_df <- left_join(transcription_factors, df, by = "Gene_name")
  
  # Write the list of DE transcription factors to file
  write.csv(final_df, file = sprintf("%s/DEG_lists/%s_tfs.csv", mainDir, output_name))
  print("Transcription factors have been identified for the supplied set of differentially expressed genes.")
  print(sprintf("The resultant file can be found at '%s/DEG_lists/%s_tfs.csv'", mainDir, output_name))
  print("------------------------------")
}