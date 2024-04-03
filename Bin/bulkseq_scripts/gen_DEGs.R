##### Generate extra DEG lists #####
# This script contains a function to generate lists of unique and common DEGs based on two input lists of DEGs. 
# It then saves these to file.
compare_DEGs <- function(mainDir, # path to the working directory
                         DEG_list1, # path to the first .csv containing DEGs
                         DEG_list2, # path to the second .csv containing DEGs
                         orig_names # a character vector indicating how you would like to name the 2 original lists of DEGs
                         ) 
  {
  
  setwd(mainDir)
  DEGs_list <- list() # create an empty list
  res <- c(DEG_list1, DEG_list2)
  for (i in res) {
    df <- read.csv(sprintf("%s", i), row.names = 1) #generate a df from a saved .csv file
    colnames(df) <- c("Gene_name", "log2FoldChange", "padj") # set the column names of the data frame
    rownames(df) <- df$Gene_name
    DEGs_list[[length(DEGs_list)+1]] <- df # add the df to the list in a new position
  }
  
  # Add three new dataframes to DEG_list, each of which are detailed below:
  
  # Determines which genes are common between the two DEG sets and saves it to slot 3 in the list.
  DEGs_list[[length(DEGs_list)+1]] <- DEGs_list[[1]][DEGs_list[[1]][["Gene_name"]] %in% DEGs_list[[2]][["Gene_name"]], ]
  
  # Determines which genes are unique to the first DEG set and saves it to slot 4 in the list
  DEGs_list[[length(DEGs_list)+1]] <- DEGs_list[[1]][!DEGs_list[[1]][["Gene_name"]] %in% DEGs_list[[2]][["Gene_name"]], ]
  
  # Determines which genes are unique to the second DEG set and saves it to slot 5 in the list
  DEGs_list[[length(DEGs_list)+1]] <- DEGs_list[[2]][!DEGs_list[[2]][["Gene_name"]] %in% DEGs_list[[1]][["Gene_name"]], ]
  
  names(DEGs_list) <- c(orig_names[1], 
                        orig_names[2],
                        sprintf("%s_%s_common", orig_names[1], orig_names[2]),
                        sprintf("%s_unique", orig_names[1]),
                        sprintf("%s_unique", orig_names[2])) # name each data frame within the list appropriately
  
  # Save this compound list containing 5 sets of DEGs to file
  saveRDS(DEGs_list, file = sprintf("%s/DEG_lists/%s_and_%s_comparison.rds", mainDir, orig_names[1], orig_names[2]))
  
  # Save individual .csv's to file
  
  # Common DEGs
  write.csv(DEGs_list[[3]], file = sprintf("%s/DEG_lists/%s_%s_common.csv", mainDir, orig_names[1], orig_names[2]))
  
  # DEGs unique to the first set of DEGs
  write.csv(DEGs_list[[4]], file = sprintf("%s/DEG_lists/%s_unique.csv", mainDir, orig_names[1]))
  
  # DEGs unique to the second set of DEGs
  write.csv(DEGs_list[[5]], file = sprintf("%s/DEG_lists/%s_unique.csv", mainDir, orig_names[2]))
  
  print(sprintf("Common and unique DEGs have been identified for the %s and %s lists of DEGs.", orig_names[1], orig_names[2]))
  print("------------------------------")
}