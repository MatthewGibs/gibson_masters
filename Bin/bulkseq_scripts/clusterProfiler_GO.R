################### clusterProfiler Script ###################
# This script contains two functions detailed below:
# - cluster_Profiler_GO: This function performs individual GO enrichment (enrichGO) on each set of genes given in a list format.
# - compare_cluster_GO: This function performs collective GO enrichment (compareCluster) on a list containing two or more gene sets.
# Whilst they are similar, these functions produce different results for the purposes of visualization. Both should typically be run.

########## Individual GO enrichment function ##########
# This function utilizes clusterProfiler's 'enrichGO' function to perform GO enrichment on individual sets of genes.
cluster_Profiler_GO <- function(mainDir, # The path to the working directory
                                DESeq2_output, # The path to the .rds files output from the DESeq2 script
                                DEG_list, # An object of the 'list' type containing dataframes with a column named 'Gene_name'. The Gene_name column should contain gene names in the ENSG ID format.
                                lfc_type = "ALL" # choose between "UP", "DOWN", or "BOTH". Selecting "ALL" will perform all three, and is the default option.
                                )
{
  ##### Setup #####
  # required packages
  clusterProfiler_packages <- c("readr",
                                "clusterProfiler",
                                "org.Hs.eg.db",
                                "DOSE",
                                "DESeq2",
                                "data.table", 
                                "ggplot2", 
                                "tidyverse", 
                                "dplyr", 
                                "enrichplot", 
                                "forcats")
  
  suppressMessages(lapply(clusterProfiler_packages, require, character.only = TRUE))
  
  setwd(mainDir)
  
  # Read in necessary functions
  source(sprintf("%s/biomaRt_annotate.R", bin_path)) # performs annotation
  source(sprintf("%s/GO_bulk.R", bin_path)) # performs clusterProfiler's enrichGO function
  
  generate_ensembl(mainDir)
  
  # Read in the set of background genes
  background_genes <- readRDS(sprintf("%s/universe.rds", DESeq2_output)) # Pull the list of universal genes from file.
  
  # Match the ENSG IDs to gene symbols
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = background_genes # character vector containing the 'initial ID' type that you want to convert
  )
  # Finalize the conversion to gene symbols
  background_genes <- object_final$external_gene_name
  
  ##### Create all the necessary directories for outputs #####
  # Creates the primary output directory for this script
  output_dir <- "clusterProfiler"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'clusterProfiler' directory.")
  }
  
  output_dir <- "clusterProfiler/GO_enrichment"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'GO_enrichment' directory. Output files will be saved there.")
  } else {
    print("This script outputs files into the 'clusterProfiler/GO_enrichment' directory.")
  }
  
  # Create the 'plots' directory if necessary.
  output_dir <- "plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created 'plots' directory.")
  }
  
  # Create the 'clusterProfiler_plots' directory if necessary.
  output_dir <- "plots/clusterProfiler_plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'clusterProfiler_plots' directory. Figures generated from clusterProfiler can be found there.")
  } else {
    print("Figures generated from clusterProfiler can be found in the 'plots/clusterProfiler_plots' directory.")
  }
  
  # Create the 'plotdata' directory if necessary.
  output_dir <- "plots/plotdata"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'plotdata' directory. Data used to generate figures can be found there.")
  } else {
    print("Data used to generate figures can be found in the 'plots/plotdata' directory.")
  }
  
  rm(output_dir)
  
  print("All the necessary directories have been created or already exist, proceeding with clusterProfiler.")
  print("------------------------------")
  
  ##### Perform GO enrichment #####
  
  # This for loop extracts the ENSG IDs from each df in 'DEGs_list'. It then converts the IDs to symbols, and performs GO enrichment.
  DEG_groups <- names(DEG_list)
  for (i in DEG_groups) {
    
    DEGs <- DEG_list[[i]]
    
    convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "external_gene_name", # The ID type you want in your final object.
             object = DEGs$Gene_name # character vector containing the 'initial ID' type that you want to convert
    )
    
    # Replace the gene IDs with the gene symbols
    DEGs[1] <- object_final[,2]
    
    # Creates a directory in which the results of the GO enrichment will be deposited.
    output_dir <- sprintf("clusterProfiler/GO_enrichment/%s", i)
    if (!dir.exists(output_dir)){
      dir.create(output_dir)
      print(sprintf("Created the 'clusterProfiler/GO_enrichment/%s' directory.", i))
    }
    
    print(sprintf("Performing GO enrichment on the '%s' set of differentially expressed genes.", i))
    # This if loop checks the specified log2fold change type. If 'ALL' was indicated when calling the cluster_Profiler funtion
    # then it runs the GO_enrich function for all 3 permutations of log2fold change type.
    if (lfc_type == "ALL"){
      Regs <- c("UP", "DOWN", "BOTH")
      for (k in Regs) {
        Regulation <- k
        # Perform the enrichment
        GO_enrich(DEGs = DEGs, # List of differentially expressed genes.
                  background_genes =  background_genes,
                  Regulation = Regulation, # choose between "UP", "DOWN", or "BOTH".
                  GO_Term = "BP", # choose between "BP" (biological process), "CC" (cellular component) or "MF" (molecular function).
                  output_dir = output_dir) # This directory is named based on the 'results' vector
      }
    } else {
      Regulation <- lfc_type
      # Perform the enrichment
      GO_enrich(DEGs = DEGs, # List of differentially expressed genes.
                background_genes =  background_genes,
                Regulation = Regulation, # choose between "UP", "DOWN", or "BOTH".
                GO_Term = "BP", # choose between "BP" (biological process), "CC" (cellular component) or "MF" (molecular function).
                output_dir = output_dir) # This directory is named based on the 'results' vector
    }
    
    print(sprintf("GO enrichment complete for the '%s' genes. Results can be found in the 'clusterProfiler/%s_enrichment' directory.", i, i))
    print("------------------------------")
  }
  print("GO enrichment complete for all gene lists.")
  print("------------------------------")
}

########## Compare Cluster function ##########
compare_cluster_GO <- function(mainDir, # The path to the working directory
                               DESeq2_output, # The path to the .rds files output from the DESeq2 script
                               DEG_list, # An object of the 'list' type containing dataframes with a column named 'Gene_name'. The Gene_name column should contain gene names in the ENSG ID format.
                               comparison_name # A character indicating what the results of this comparison should be named.
                               )
  {
  ##### Setup #####
  # required packages
  clusterProfiler_packages <- c("readr",
                                "clusterProfiler",
                                "org.Hs.eg.db",
                                "DOSE",
                                "DESeq2",
                                "ggplot2",
                                "enrichplot")
  
  suppressMessages(lapply(clusterProfiler_packages, require, character.only = TRUE))
  
  setwd(mainDir)
  
  # Read in necessary functions
  source(sprintf("%s/biomaRt_annotate.R", bin_path)) # performs annotation
  
  generate_ensembl(mainDir)
  
  # Read in the set of background genes
  background_genes <- readRDS(sprintf("%s/universe.rds", DESeq2_output)) # Pull the list of universal genes from file.
  
  # Match the ENSG IDs to gene symbols
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = background_genes # character vector containing the 'initial ID' type that you want to convert
  )
  
  # Finalize the conversion to gene symbols
  background_genes <- object_final$external_gene_name
  
  ##### Create all the necessary directories for outputs #####
  # Creates the primary output directory for this script
  output_dir <- "clusterProfiler"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'clusterProfiler' directory.")
  }
  
  output_dir <- "clusterProfiler/GO_enrichment"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'GO_enrichment' directory. Output files will be saved there.")
  } else {
    print("This script outputs files into the 'clusterProfiler/GO_enrichment' directory.")
  }
  
  # Create the 'plots' directory if necessary.
  output_dir <- "plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created 'plots' directory.")
  }
  
  # Create the 'clusterProfiler_plots' directory if necessary.
  output_dir <- "plots/clusterProfiler_plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'clusterProfiler_plots' directory. Figures generated from clusterProfiler can be found there.")
  } else {
    print("Figures generated from clusterProfiler can be found in the 'plots/clusterProfiler_plots' directory.")
  }
  
  # Create the 'plotdata' directory if necessary.
  output_dir <- "plots/plotdata"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'plotdata' directory. Data used to generate figures can be found there.")
  } else {
    print("Data used to generate figures can be found in the 'plots/plotdata' directory.")
  }
  
  rm(output_dir)
  
  print("All the necessary directories have been created or already exist, proceeding with clusterProfiler.")
  print("------------------------------")
  
  ##### Perform GO enrichment through cluster Comparisons #####
  
  # This for loop replaces the ENSG IDs from each df in 'DEGs_list' with the corresponding gene symbol.
  DEG_groups <- names(DEG_list)
  
  # Create an empty list that will be filled with annotated gene names.
  sublist <- list()
  
  for (i in DEG_groups) {
    
    DEGs <- DEG_list[[i]]
    
    convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "external_gene_name", # The ID type you want in your final object.
             object = DEGs$Gene_name # character vector containing the 'initial ID' type that you want to convert
    )
    
    # Write the gene symbols into the new list
    sublist[[i]] <- object_final$external_gene_name
  }
  
  print("Running compareCluster, this may take some time.")
  # Run the compareCluster funtion.
  cclust <- compareCluster(geneCluster = sublist, 
                           fun = enrichGO,
                           universe = background_genes,
                           OrgDb = org.Hs.eg.db,
                           keyType = "SYMBOL",
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.01,
                           qvalueCutoff = 0.05,	
                           readable = TRUE)
  print("compareCluster is complete, simplifying GO enrichments and performing visualisations.")
  
  # Use 'simplify' to reduce redundancy of GO terms.
  cclust <- clusterProfiler::simplify(x = cclust, cutoff = 0.7)
  
  # Generate a dotplot of GO terms across the different lists.
  cclust_dotplot <- dotplot(cclust, showCategory = 10)
  
  ggsave(sprintf("%s/plots/clusterProfiler_plots/%s_comparecluster_dotplot.png", mainDir, comparison_name), 
         bg = "white",
         dpi = 300,
         width = 10,
         height = 18,
         units = "in")
  
  # Save the 'dotplot' object to file.
  saveRDS(cclust_dotplot, file = sprintf("%s/plots/plotdata/%s_comparecluster_dotplot.rds", mainDir, comparison_name))
  print(sprintf("The dotplot has been saved. It can be found at: %s/plots/clusterProfiler_plots/%s_comparecluster_dotplot.png", mainDir, comparison_name))
  
  # Calculate pairwise termism between GO terms.
  cclust <- pairwise_termsim(cclust)
  
  cclust_emapplot <- emapplot(cclust, pie.params = list(legend_n = 2))
  
  ggsave(sprintf("%s/plots/clusterProfiler_plots/%s_comparecluster_emapplot.png", mainDir, comparison_name), 
         bg = "white",
         dpi = 300,
         width = 10,
         height = 10,
         units = "in")
  
  # Save the 'emapplot' object to file.
  saveRDS(cclust_emapplot, file = sprintf("%s/plots/plotdata/%s_comparecluster_emapplot.rds", mainDir, comparison_name))
  print(sprintf("The emapplot has been saved. It can be found at: %s/plots/clusterProfiler_plots/%s_comparecluster_emapplot.png", mainDir, comparison_name))
  
  # Save the 'cclust' object to file.
  saveRDS(cclust, file = sprintf("%s/clusterProfiler/GO_enrichment/%s_comparecluster.rds", mainDir, comparison_name))
  
  print("All plots and objects have been saved to file.")
  print("------------------------------")
}
