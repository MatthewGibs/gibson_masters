# This script contains the function to run tximport on salmon output.

tximport_salmon <- function(mainDir, # path to the 'work' directory
                            quant_directories, # CHARACTER vector indicating the location of all quant.sf files to import in the format "folder.quant/folder"
                            tx2gene_path, # provide the path to the file that matches transcript and gene IDs
                            sample_names, # CHARACTER vector naming each sample. This must be in the same order as the directories given in 'quant_directories'
                            annotation_type = "external_gene_name") # Optional variable for setting the type of gene ID or name as row names in the count matrix. The default is "external_gene_name".
{
  ##### Setup #####
  # Required packages
  import_packages <- c("tximport",
                       "biomaRt",
                       "dplyr",
                       "GenomicFeatures")
  
  suppressMessages(lapply(import_packages, require, character.only = TRUE))

  setwd(mainDir)
  
  # Set path to the quant.sf files of all listed samples
  files <- file.path(mainDir, "salmon", quant_directories, "quant.sf")
  
  ##### Acquire tx2gene #####
  # Check if there is an existing tx2gene file. If not, create one and save it to file for future use.
  if (file.exists(sprintf("%s", tx2gene_path))) {
    
    print("The 'tx2gene' file already exists, importing file.")
    
    # Read in the transcripts to gene file.
    tx2gene <- read.table(tx2gene_path, header = FALSE, col.names = c("TXNAME", "GENEID")) 
    
  } else {
    
    print("The 'tx2gene' file does not already exist. The file will now be produced, this may take some time.")
    
    # Create the tx2gene object matching transcript and gene IDs.
    txdb <- makeTxDbFromGFF("~/Matthew_Masters/Raw/Data/bulkseq_data/gencode.v44.genome.gtf.gz")
    k <- keys(txdb, keytype = "TXNAME")
    tx2gene <- select(txdb, k, "GENEID", "TXNAME")
    
    print("writing 'tx2gene.tsv' to file")
    
    # Save tx2gene to file.
    write.table(tx2gene, file = sprintf("%s", tx2gene_path), row.names = FALSE, col.names = FALSE, sep="\t")
  }
  
  ##### Run tximport #####
  # Run tximport specifying 'salmon' as the file source.
  txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)
  
  # Set the column names of the count matrix as specified by 'sample_names'.
  colnames(txi.salmon$counts) <- sample_names 
  
  # Create 'cts', an object used as input for the annotation script.
  cts <- txi.salmon$counts
  
  ##### Annotation #####
  # If "ensembl_gene_id" is specified as the annotation type, this means that the 'annotate' function is unnecessary.
  if (annotation_type == "ensembl_gene_id") {
    # Cut off the version numbers present in the rownames of the 'cts' object.
    rownames(cts) <- substr(rownames(cts), 1, 15)
    print("Ensembl IDs are being retained, skipping annotation.")
  } else {
    # Read in the function for acquiring annotation info
    source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
    generate_ensembl(mainDir)
    
    genes <- substr(rownames(cts), 1, 15)
    
    convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
               final_ID = annotation_type, # The ID type you want in your final object.
               object = genes # character vector containing the 'initial ID' type that you want to convert
    )
    
    print(sprintf("Changing the row names of the count matrix from ensembl_gene_id to %s", annotation_type))
    # Finally to change the rownames of the count matrix to the external gene names
    rownames(cts) <- object_final[,2]
  }
  
  txi.salmon$counts <- cts
  
  ##### Write output to file #####
  # Create the 'tximport' directory if it doesnt already exist.
  output_dir <- "tximport"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'tximport' directory. The count matrix will be saved there.")
  }
  
  # Write the output to file
  saveRDS(object = txi.salmon, 
          file = sprintf("%s/tximport/txi.count.matrix.rds", mainDir))
  
  print(sprintf("The output count matrix has been saved to the '%s/tximport' directory.", mainDir))
  print("tximport complete")
  print("------------------------------")
}