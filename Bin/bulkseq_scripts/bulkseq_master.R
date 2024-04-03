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

##### Key parameters ####
# This pipelines operates based on a specific directory structure using a number of scripts and files.
# The following paths are used but can be modified to suit your needs:

bin_path <- "~/git/gibson_masters/Bin/bulkseq_scripts" # Set path to the 'Bin' directory which contains essential scripts
mainDir <- "~/git/gibson_masters/Work" # Set the path to the 'Work' directory, which stores most outputs


##### tximport #####
# Call script to generate a count matrix from quant.sf files using Tximport.
source(sprintf("%s/tximport_script.R", bin_path))

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

# Run the tximport function
tximport_salmon(mainDir, # path to the 'work' directory
                quant_directories, # CHARACTER vector indication the location of all quant.sf files to import in the format "folder.quant/folder"
                tx2gene_path, # provide the path to the file that matches transcript and gene IDs
                sample_names, # CHARACTER vector naming each sample. This must be in the same order as the directories given in 'quant_directories'
                annotation_type = "ensembl_gene_id" # Optional variable for setting the type of gene ID or name as row names in the count matrix. The default is "external_gene_name".
                )
if(exists("object_final")) {rm(object_final)}

##### Metadata #####
# Read in the script for generating metadata.
source(sprintf("%s/gen_metadata.R", bin_path))

# Call the funtion to generate metadata
generate_metadata(mainDir = mainDir)

##### DESeq2 #####
# Perform Differential Expression Analysis using DESeq2.
source(sprintf("%s/DESeq2_script.R", bin_path))

# path to a tximport output file (saved in .rds format)
txi_mtx <- sprintf("%s/tximport/txi.count.matrix.rds", mainDir)

# path to the metadata file DESeq will use to inform design
meta_file <- sprintf("%s/Annotation/sample_metadata.rds", mainDir)

# the column of the sample_metadata file which will be used for the design of the DEA
design <- ~ treatment

DESeq2_DEA(mainDir = mainDir, # path to the 'work' directory
           txi_mtx = txi_mtx, # name of the object containing the metadata for the txi_mtx
           meta_file = meta_file, # path to the file containing the metadata information for the txi_mtx
           design = design,
           lfcThreshold = 1)

##### Generate lists of common and unique DEGs #####
source(sprintf("%s/gen_DEGs.R", bin_path))

DEG_list1 <- sprintf("%s/DEG_lists/treatment_GMCSF_vs_none_DEGs.csv", mainDir)

DEG_list2 <- sprintf("%s/DEG_lists/treatment_MCSF_vs_none_DEGs.csv", mainDir)

orig_names <- c("GMCSF", "MCSF")

compare_DEGs(mainDir, # path to the working directory
             DEG_list1, # path to the first .csv containing DEGs
             DEG_list2, # path to the second .csv containing DEGs
             orig_names # a character vector indicating how you would like to name the 2 original lists of DEGs
             )

##### Identify transcription factors #####
source(sprintf("%s/identify_tfs.R", bin_path))
source(sprintf("%s/gen_DEGs.R", bind_path))

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
source(sprintf("%s/clusterProfiler_GO.R", bin_path))

DESeq2_output <- sprintf("%s/DESeq2", mainDir)

# Read in the list containing all sets of DEGs
DEG_list <- readRDS(sprintf("%s/DEG_lists/GMCSF_and_MCSF_comparison.rds", mainDir)) # Pull the list of universal genes from file.

cluster_Profiler_GO(mainDir, # The path to the working directory
                    DESeq2_output, # The path to the .rds files output from the DESeq2 script
                    DEG_list, # An object of the 'list' type containing dataframes with a column named 'Gene_name'. The Gene_name column should contain gene names in the ENSG ID format.
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
source(sprintf("%s/gen_GO_summary.R", bin_path))
summarise_GO(mainDir = mainDir)

##### Identify key transcription factors and their interacting partners #####
source(sprintf("%s/tf_interactions.R", bin_path))

collate_tf_DEGs(mainDir)

##### Generate Figures and tables #####
source(sprintf("%s/dissertation_plots.R", bin_path))

gen_all_plots(mainDir)
