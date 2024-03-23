################# ClusterProfiler KEGG pathway script #####################
# This function utilizes clusterProfiler's 'enrichKEGG' function to perform GO enrichment on individual sets of genes.
cluster_Profiler_KEGG <- function(mainDir, # The path to the working directory
                                DESeq2_output, # The path to the .rds files output from the DESeq2 script
                                DEG_list # An object of the 'list' type containing dataframes with a column named 'Gene_name'. The Gene_name column should contain gene names in the ENSG ID format.
)
{
  ##### Setup #####
  KEGG_packages <- c("clusterProfiler", 
                     "org.Hs.eg.db",
                     "R.utils", 
                     "biomaRt",
                     "dplyr")
  
  suppressMessages(lapply(KEGG_packages, require, character.only = TRUE))
  
  setwd(mainDir)
  
  # Read in necessary functions
  source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R") # performs annotation
  
  generate_ensembl(mainDir)
  
  # Read in the set of background genes
  background_genes <- readRDS(sprintf("%s/universe.rds", DESeq2_output)) # Pull the list of universal genes from file.
  
  # Match the ENSG IDs to Entrez gene IDs
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "entrezgene_id", # The ID type you want in your final object.
             object = background_genes # character vector containing the 'initial ID' type that you want to convert
  )
  
  # Finalize the conversion to entrez ID
  background_genes <- object_final$entrezgene_id
  
  # Ensure that clusterProfiler can access the KEGG database
  R.utils::setOption("clusterProfiler.download.method","auto")
  
  ##### Create all the necessary directories for outputs #####
  # Creates the primary output directory for this script
  output_dir <- "clusterProfiler"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'clusterProfiler' directory.")
  }
  
  output_dir <- "clusterProfiler/KEGG_enrichment"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'KEGG_enrichment' directory. Output files will be saved there.")
  } else {
    print("This script outputs files into the 'clusterProfiler/KEGG_enrichment' directory.")
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
  
  print("All the necessary directories have been created or already exist, proceeding with KEGG enrichment.")
  print("------------------------------")
  
  ##### Perform KEGG enrichment #####
  # This for loop replaces the ENSG IDs from each df in 'DEGs_list' with the corresponding Entrez ID.
  DEG_groups <- names(DEG_list)
  
  for (i in DEG_groups) {
    
    DEGs <- DEG_list[[i]]
    
    # Match the ENSG IDs to Entrez gene IDs
    convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
               final_ID = "entrezgene_id", # The ID type you want in your final object.
               object = DEGs$Gene_name # character vector containing the 'initial ID' type that you want to convert
    )
    
    # Finalize the conversion to entrez ID
    DEGs <- object_final$entrezgene_id
    
    # KEGG pathway over-representation analysis
    kk <- enrichKEGG(gene = DEGs,
                     organism = "hsa",
                     pvalueCutoff = 0.05,
                     universe = background_genes,
                     qvalueCutoff = 0.05)
    
    results <- kk@result
    results <- subset(results, `p.adjust` < 0.05)
    
    write.csv(results, file = sprintf("clusterProfiler/KEGG_enrichment/%s_sig_KEGG.csv", i))
  }
  print("KEGG enrichment complete for all supplied lists of genes.")
  print("------------------------------")
}

########## Function for plotting a specifed KEGG pathway ##########
plot_pathway <- function(mainDir, # The path to the working directory.
                         DEG_output, # An dataframe containing ENSG ID's in a Gene_name column, as well as matching log2fold changes in a log2FoldChange column.
                         pathway_IDs, # A character vector containing eight-digit characters specifying the pathway(s) to visualize. (example: "hsa04657")
                         test_name # A character indicating the name of the DEGs provided, usually treatment type.
                         )
  {
  ##### Setup #####
  KEGG_packages <- c("clusterProfiler", 
                     "AnnotationDbi", 
                     "org.Hs.eg.db",
                     "R.utils", 
                     "biomaRt",
                     "dplyr",
                     "pathview")
  
  suppressMessages(lapply(KEGG_packages, require, character.only = TRUE))
  
  setwd(mainDir)
  
  ##### Create all the necessary directories for outputs #####
  # Create the 'plots' directory if necessary.
  output_dir <- "plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created 'plots' directory.")
  }
  
  # Create the 'KEGG_plots' directory if necessary.
  output_dir <- "plots/KEGG_plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'KEGG_plots' directory. Figures generated from KEGG enrichment can be found there.")
  } else {
    print("Figures generated from KEGG enrichment can be found in the 'plots/KEGG_plots' directory.")
  }
  
  # Create the sub-directory for the test name if necessary.
  output_dir <- sprintf("plots/KEGG_plots/%s", test_name)
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print(sprintf("Created the 'KEGG_plots/%s' directory. Figures generated from this pathview run can be found there.", test_name))
  } else {
    print(sprintf("Figures generated this pathview can be found in the 'plots/KEGG_plots/%s' directory.", test_name))
  }
  
  rm(output_dir)
  
  print("All the necessary directories have been created or already exist, proceeding with KEGG pathway visualisation.")
  print("------------------------------")
  
  ##### Manipulate the input data #####
  genes <- DEG_output$log2FoldChange
  
  # Read in necessary functions
  source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R") # performs annotation
  
  generate_ensembl(mainDir)
  
  # Match the ENSG IDs to Entrez gene IDs
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "entrezgene_id", # The ID type you want in your final object.
             object = DEG_output$Gene_name # character vector containing the 'initial ID' type that you want to convert
  )
  
  # Finalize the conversion to entrez ID
  entrez_IDs <- object_final$entrezgene_id
  
  names(genes) <- entrez_IDs
  
  ##### Run Pathview for all IDs in the pathway_IDs vector #####
  for (i in pathway_IDs) 
    {
    specific_ID <- i
    
    # Create the sub-directory for the specified pathway if necessary.
    output_dir <- sprintf("plots/KEGG_plots/%s/%s", test_name, specific_ID)
    if (!dir.exists(output_dir)){
      dir.create(output_dir)
    }
    
    # As the pathview function creates the image in the current directory, we must switch to where it must be saved.
    setwd(sprintf("%s/plots/KEGG_plots/%s/%s", mainDir, test_name, specific_ID))
    
    # Certain pathways cannot be visualised. The block of code below allows these errors to be reported and the loop continues.
    tryCatch(
      {
      hsa <- pathview(
        gene.data  = genes,
        pathway.id = specific_ID,
        species    = "hsa",
        limit      = list(gene = 10, cpd = 1), 
        low        = "steelblue4",
        mid        = "grey90",
        high       = "coral"
      )
      
      saveRDS(hsa, file = sprintf("%s/plots/KEGG_plots/%s/%s/%s_%s_pathview.rds", mainDir, test_name, specific_ID, test_name, specific_ID))
      
      print(sprintf("Pathway %s visualised successfully.", specific_ID))
      
      
    }, error = function(e) {
      print(sprintf("There was an error visualising pathway %s: %s", specific_ID, conditionMessage(e)))
    })
    setwd(mainDir)
    print("------------------------------")
  }
  
}
