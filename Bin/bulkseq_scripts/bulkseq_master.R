##### Master bulk RNA-seq script #####
# Written by Matthew Gibson
# This script runs an analysis pipeline comprising tximport, DESeq2 and clusterProfiler for the purpose of DEA.
##### Pipeline Structure #####
# This RNA-seq pipeline contains the following modules:
# 1: Run tximport on the output of the salmon quasi-mapping package.
# 2: Generate a metadata file for the default dataset used by this pipeline.
# 3: Run DESeq2 on the count matrix produced by Tximport.
# 4: Run clusterProfiler for GO enrichment.
# 5: Run clusterProfiler for KEGG enrichment.
# 6: Produce auxiliary figures based on all of the above.
##### tximport #####
# Generate a count matrix from quant.sf files using Tximport.
source("~/Matthew_Masters/Bin/bulkseq_scripts/tximport_script.R")

# set the path to the working directory
mainDir <- "~/Matthew_Masters/Work"

# Create a character vector containing the names of each sample
quant_directories <- c("SRR16816545", 
                       "SRR16816546", 
                       "SRR16816547",
                       "SRR16816548", 
                       "SRR16816549", 
                       "SRR16816550", 
                       "SRR16816551", 
                       "SRR16816552", 
                       "SRR16816553")

# CHARACTER vector naming each sample. This must be in the same order as the directories given in 'quant_directories'
sample_names <- c("donor1_untreated",
                  "donor2_untreated",
                  "donor3_untreated",
                  "donor1_MCSF",
                  "donor2_MCSF",
                  "donor3_MCSF",
                  "donor1_GMCSF",
                  "donor2_GMCSF",
                  "donor3_GMCSF")

tx2gene_path <- c("~/Matthew_Masters/Raw/Data/bulkseq_data/txp2gene.tsv")

# Run the tximport funtion
tximport_salmon(mainDir, # path to the 'work' directory
                quant_directories, # CHARACTER vector indication the location of all quant.sf files to import in the format "folder.quant/folder"
                tx2gene_path, # provide the path to the file that matches transcript and gene IDs
                sample_names, # CHARACTER vector naming each sample. This must be in the same order as the directories given in 'quant_directories'
                annotation_type = "ensembl_gene_id" # Optional variable for setting the type of gene ID or name as row names in the count matrix. The default is "external_gene_name".
                )
if(exists("object_final")) {rm(object_final)}

##### Metadata #####
# Read in the script for generating metadata.
source("~/Matthew_Masters/Bin/bulkseq_scripts/gen_metadata.R")

# set the path to the working directory
mainDir <- "~/Matthew_Masters/Work" 

# Call the funtion to generate metadata
generate_metadata(mainDir = mainDir)

##### DESeq2 #####
# Perform Differential Expression Analysis using DESeq2.
source("~/Matthew_Masters/Bin/bulkseq_scripts/DESeq2_script.R")

mainDir <- "~/Matthew_Masters/Work" # set the path to the working directory

# path to a tximport output file (saved in .rds format)
txi_mtx <- "~/Matthew_Masters/Work/tximport/txi.count.matrix.rds"

# path to the metadata file DESeq will use to inform design
meta_file <- "~/Matthew_Masters/Work/Annotation/sample_metadata.rds"

# the column of the sample_metadata file which will be used for the design of the DEA
design <- ~ treatment

DESeq2_DEA(mainDir = mainDir, # path to the 'work' directory
           txi_mtx = txi_mtx, # name of the object containing the metadata for the txi_mtx
           meta_file = meta_file, # path to the file containing the metadata information for the txi_mtx
           design = design,
           lfcThreshold = 1)

##### Generate extra lists of DEGs #####
source("~/Matthew_Masters/Bin/bulkseq_scripts/gen_DEGs.R")

mainDir <- "~/Matthew_Masters/Work" # set the path to the working directory

DEG_list1 <- sprintf("%s/DEG_lists/treatment_GMCSF_vs_none_DEGs.csv", mainDir)

DEG_list2 <- sprintf("%s/DEG_lists/treatment_MCSF_vs_none_DEGs.csv", mainDir)

orig_names <- c("GMCSF", "MCSF")

compare_DEGs(mainDir, # path to the working directory
             DEG_list1, # path to the first .csv containing DEGs
             DEG_list2, # path to the second .csv containing DEGs
             orig_names # a character vector indicating how you would like to name the 2 original lists of DEGs
             )

##### Identify transcription factors #####
source("~/Matthew_Masters/Bin/bulkseq_scripts/identify_tfs.R")
source("~/Matthew_Masters/Bin/bulkseq_scripts/gen_DEGs.R")

mainDir <- "~/Matthew_Masters/Work" # set the path to the working directory

# Identify the transcription factors within the list of M-CSF differentially expressed genes.
identify_tfs(mainDir, # path to the working directory
             DEG_list = sprintf("%s/DEG_lists/treatment_GMCSF_vs_none_DEGs.csv", mainDir), # path to a .csv containing DEGs
             output_name = "GMCSF" # The prefix used to name the file, typically a treatment/condition type that dictated the DEGs.
)

# Identify the transcription factors within the list of M-CSF differentially expressed genes.
identify_tfs(mainDir, # path to the working directory
             DEG_list = sprintf("%s/DEG_lists/treatment_MCSF_vs_none_DEGs.csv", mainDir), # path to a .csv containing DEGs
             output_name = "MCSF" # The prefix used to name the file, typically a treatment/condition type that dictated the DEGs.
)

# Colate the above lists and identify common and unique transcription factors.
DEG_list1 <- sprintf("%s/DEG_lists/GMCSF_tfs.csv", mainDir)

DEG_list2 <- sprintf("%s/DEG_lists/MCSF_tfs.csv", mainDir)

orig_names <- c("GMCSF_tfs", "MCSF_tfs")

compare_DEGs(mainDir, # path to the working directory
             DEG_list1, # path to the first .csv containing DEGs
             DEG_list2, # path to the second .csv containing DEGs
             orig_names # a character vector indicating how you would like to name the 2 original lists of DEGs
)

##### clusterProfiler GO Enrichment #####
# Performs GO enrichment through clusterProfiler.
source("~/Matthew_Masters/Bin/bulkseq_scripts/clusterProfiler_GO.R")

mainDir <- "~/Matthew_Masters/Work"

DESeq2_output <- "~/Matthew_Masters/Work/DESeq2"

# Read in the list containing all sets of DEGs
DEG_list <- readRDS(sprintf("%s/DEG_lists/GMCSF_and_MCSF_comparison.rds", mainDir)) # Pull the list of universal genes from file.

cluster_Profiler_GO(mainDir, # The path to the working directory
                    DESeq2_output, # The path to the .rds files output from the DESeq2 script
                    DEG_list, # An object of the 'list' type containing dataframes with a column named 'Gene_name'. The Gene_name column should contain gene names in the ENSG ID format.
                    )

# Run compareCluster on the full list of 5 DEG sets. The DEG lists are reordered for clarity.

# Read in the list containing all sets of DEGs
DEG_list <- readRDS(sprintf("%s/DEG_lists/GMCSF_and_MCSF_comparison.rds", mainDir)) # Pull the list of universal genes from file.
DEG_list <- list(DEG_list[[1]], DEG_list[[4]], DEG_list[[3]], DEG_list[[5]], DEG_list[[2]])
names(DEG_list) <- c("GM-CSF", "GM-CSF_unique", "Common", "M-CSF_unique", "M-CSF")

# Run the function
compare_cluster_GO(mainDir, # The path to the working directory
                   DESeq2_output, # The path to the .rds files output from the DESeq2 script
                   DEG_list, # An object of the 'list' type containing dataframes with a column named 'Gene_name'. The Gene_name column should contain gene names in the ENSG ID format.
                   comparison_name = "All" # A character indicating what the results of this comparison should be named.
)

# Run compareCluster on the DEGs for GMCSF and MCSF. This requires subsetting of the DEG_list object.

# Read in the list containing all sets of DEGs
DEG_list <- readRDS(sprintf("%s/DEG_lists/GMCSF_and_MCSF_comparison.rds", mainDir)) # Pull the list of universal genes from file.
DEG_list <- list(DEG_list[[1]], DEG_list[[2]])
names(DEG_list) <- c("GM-CSF", "M-CSF")

# Run the function
compare_cluster_GO(mainDir, # The path to the working directory
                   DESeq2_output, # The path to the .rds files output from the DESeq2 script
                   DEG_list, # An object of the 'list' type containing dataframes with a column named 'Gene_name'. The Gene_name column should contain gene names in the ENSG ID format.
                   comparison_name = "GMCSF_MCSF" # A character indicating what the results of this comparison should be named.
)

# Create a summary list of all generated results
source("~/Matthew_Masters/Bin/bulkseq_scripts/gen_GO_summary.R")
summarise_GO(mainDir = "~/Matthew_Masters/Work")

##### clusterProfiler KEGG enrichment #####
# Performs KEGG enrichment through clusterProfiler.
source("~/Matthew_Masters/Bin/bulkseq_scripts/clusterProfiler_KEGG.R")

mainDir <- "~/Matthew_Masters/Work"

DESeq2_output <- "~/Matthew_Masters/Work/DESeq2"

# Read in the list containing all sets of DEGs
DEG_list <- readRDS(sprintf("%s/DEG_lists/GMCSF_and_MCSF_comparison.rds", mainDir)) # Pull the list of universal genes from file.
names(DEG_list) <- c("GM-CSF", "M-CSF", "Common", "GM-CSF_unique", "M-CSF_unique")

# This function utilizes clusterProfiler's 'enrichKEGG' function to perform GO enrichment on individual sets of genes.
cluster_Profiler_KEGG(mainDir, # The path to the working directory
                      DESeq2_output, # The path to the .rds files output from the DESeq2 script
                      DEG_list # An object of the 'list' type containing dataframes with a column named 'Gene_name'. The Gene_name column should contain gene names in the ENSG ID format.
)

##### Visualise and save KEGG pathways #####
# Performs KEGG enrichment through clusterProfiler.
source("~/Matthew_Masters/Bin/bulkseq_scripts/clusterProfiler_KEGG.R")

mainDir <- "~/Matthew_Masters/Work"

DEG_output <- read.csv("~/Matthew_Masters/Work/DEG_lists/treatment_GMCSF_vs_none_DEGs.csv", row.names = 1)

pathway_IDs <- read.csv("~/Matthew_Masters/Work/clusterProfiler/KEGG_enrichment/GM-CSF_sig_KEGG.csv")
pathway_IDs <- pathway_IDs$ID

plot_pathway(mainDir = mainDir, # The path to the working directory.
             DEG_output = DEG_output, # An dataframe containing ENSG ID's in a Gene_name column, as well as matching log2fold changes in a log2FoldChange column.
             pathway_IDs = pathway_IDs, # A character vector containing eight-digit characters specifying the pathway(s) to visualize. (example: "hsa04657")
             test_name = "GMCSF" # A character indicating the name of the DEGs provided, usually treatment type.
)

DEG_output <- read.csv("~/Matthew_Masters/Work/DEG_lists/treatment_MCSF_vs_none_DEGs.csv", row.names = 1)

pathway_IDs <- read.csv("~/Matthew_Masters/Work/clusterProfiler/KEGG_enrichment/M-CSF_sig_KEGG.csv")
pathway_IDs <- pathway_IDs$ID

plot_pathway(mainDir = mainDir, # The path to the working directory.
             DEG_output = DEG_output, # An dataframe containing ENSG ID's in a Gene_name column, as well as matching log2fold changes in a log2FoldChange column.
             pathway_IDs = pathway_IDs, # A character vector containing eight-digit characters specifying the pathway(s) to visualize. (example: "hsa04657")
             test_name = "MCSF" # A character indicating the name of the DEGs provided, usually treatment type.
)

##### Generate Figures and tables #####

source("~/Matthew_Masters/Bin/bulkseq_scripts/dissertation_plots.R")

mainDir <- "~/Matthew_Masters/Work"

gen_all_plots(mainDir)

circular_GO_heatmap(mainDir, # path to the 'Work' directory
                    GO_ID = "GO:0050727", 
                    plot_name = "inflamm_heatmap")

circular_GO_heatmap(mainDir, # path to the 'Work' directory
                    GO_ID = "GO:0030099", 
                    plot_name = "myeloid_diff_heatmap")

# L2FC = 1 pathways
KEGG_compare_pathways(mainDir, specific_ID = "hsa05202")
KEGG_compare_pathways(mainDir, specific_ID = "hsa05203")
KEGG_compare_pathways(mainDir, specific_ID = "hsa04613")
KEGG_compare_pathways(mainDir, specific_ID = "hsa05322")

# L2FC = 2 pathways
#KEGG_compare_pathways(mainDir, specific_ID = "hsa04640")
#KEGG_compare_pathways(mainDir, specific_ID = "hsa04657")
#KEGG_compare_pathways(mainDir, specific_ID = "hsa04668")
#KEGG_compare_pathways(mainDir, specific_ID = "hsa03320")
#KEGG_compare_pathways(mainDir, specific_ID = "hsa04110")
