# This script contains the function to perform differential expression analysis using DESeq2.

# Let's create the function
DESeq2_DEA <- function(mainDir, # path to the 'work' directory
                       txi_mtx, # path to a tximport output file (saved in .rds format)
                       meta_file, # path to the file containing the metadata information for the txi_mtx
                       design, # the column of the sample_metadata file which will be used for the design of the DEA
                       padjThreshold = 0.05, # Set the threshold for padj when running DESeq2 and producing a filtered count matrix. Default is set to 0.05.
                       lfcThreshold = 2 # Set the threshold for log2FoldChange when running DESeq2 and producing a filtered count matrix. Default is set to 1.
                       )
{
  ##### Setup #####
  # Required packages
  DESeq_package <- c("DESeq2", "dplyr", "ggplot2", "RColorBrewer")
  
  suppressMessages(lapply(DESeq_package, require, character.only = TRUE))
  
  # Move to the specified main directory.
  setwd(mainDir)
  
  # Read in the count matrix.
  cts <- readRDS(txi_mtx)
  
  # Read in the metadata file
  meta <- readRDS(meta_file)
  
  ##### Create all the necessary directories for outputs #####
  # Create the 'DESeq2' directory if necessary.
  output_dir <- "DESeq2"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'DESeq2' directory. Files for downstream use are stored there.")
  } else {
    print("Files for downstream use are stored in the 'DESeq2' directory.")
  }
  
  # Create the 'plots' directory if necessary.
  output_dir <- "plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Creating the 'plots' directory.")
  }
  
  setwd(output_dir)
  # Create the 'DESeq2_plots' directory if necessary.
  output_dir <- "DESeq2_plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'DESeq2_plots' directory. Figures generated from DESeq2 can be found there.")
  } else {
    print("Figures generated from DESeq2 can be found in the 'plots/DESeq2_plots' directory.")
  }
  # Create the 'plotdata' directory if necessary.
  output_dir <- "plotdata"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'plotdata' directory. Data used to generate figures can be found there.")
  } else {
    print("Data used to generate figures can be found in the 'plots/plotdata' directory.")
  }
  
  setwd(mainDir)
  
  # Create the 'DEG_lists' directory if necessary.
  output_dir <- "DEG_lists"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'DEG_lists' directory. Lists of differentially expressed genes can be found there.")
  } else {
    print("Lists of differentially expressed genes can be found in the 'DEG_lists' directory.")
  }
  
  rm(output_dir)
  
  print("All the necessary directories have been created or already exist, proceeding with DESeq2.")
  print("------------------------------")
  
  ##### Filtering #####
  # The count matrix may contain an inflated number of genes. These genes will likely only have a small number of counts.
  # We shall filter the matrix based on rowsums greater than 5 (any gene with 5 or fewer counts is removed).
  cts_abu <- cts$abundance
  cts_mtx <- cts$counts
  cts_len <- cts$length
  
  # for each gene, compute total count and compare to threshold
  # keeping outcome in vector of 'logicals' (ie TRUE or FALSE, or NA)
  keep <- rowSums(cts_mtx) > 10
  # summary of test outcome: number of genes in each class:
  table(keep, useNA="always")
  
  # subset genes where test was TRUE.
  cts_abu <- cts_abu[keep,]
  cts_mtx <- cts_mtx[keep,]
  cts_len <- cts_len[keep,]
  
  # Apply the changes to the cts object.
  cts$abundance <- cts_abu
  cts$counts <- cts_mtx
  cts$length <- cts_len
  
  universe <- rownames(cts_mtx)
  # Save the universe of genes to file for later use.
  saveRDS(object = universe, 
          file = sprintf("%s/DESeq2/universe.rds", mainDir))
  
  print("The count matrix has been filtered to remove low counts.")
  ##### Running DESeq2 #####
  # Construct a DESeqDataSet from the cts object and sample information in samples.
  print("Constructing a DESeq dataset from tximport results")
  dds <- DESeqDataSetFromTximport(cts,
                                  colData = meta,
                                  design = design)
  print("Performing differential gene expression analysis using DESeq2")
  dds <- DESeq(dds)
  res <- resultsNames(dds)
  
  # Save the dds to file for later use in plotting counts.
  saveRDS(object = dds, 
          file = sprintf("%s/DESeq2/dds.rds", mainDir))
  
  # This loop outputs filtered results for each conditional comparison based on the experimental design.
  for (i in res) {
    if (i == "Intercept"){} 
    else {
      print(sprintf("Acquiring the list of DEGs for %s", i))
      
      
      # Build a results table based on a target padj value and log fold change threshold
      res_i <- results(dds, 
                       name = sprintf("%s", i),
                       alpha = padjThreshold,
                       lfcThreshold = lfcThreshold)
      
      # Save the object to file and name it according to its comparison details.
      saveRDS(object = res_i, 
              file = sprintf("%s/DESeq2/%s.rds", mainDir, i))
      print(sprintf("Saved %s.rds to file", i))
      
      # Shrink LFC in order to visualize data and save MA plots.
      print("Performing log fold change shrinkage using apeglm. This is for certain plots.")
      suppressMessages(resLFC <- lfcShrink(dds, 
                                           coef = sprintf("%s", i), 
                                           type = "apeglm", 
                                           lfcThreshold = lfcThreshold)
                       )
      
      
      # Create the MA plot and save it to file
      png(file = sprintf("%s/plots/DESeq2_plots/%s_MA.png", mainDir, i),
          res = 300,
          width = 5.25,
          height = 7,
          units = 'in')
      plotMA(resLFC, ylim=c(-7.5, 7.5)) + 
        abline(h=c(-lfcThreshold, lfcThreshold), col="dodgerblue", lwd=2)
      dev.off()
      
      # Save the data used to generate each MA plot to file
      saveRDS(resLFC, file = sprintf("%s/plots/plotdata/%s_MAdata.rds", mainDir, i))
      
      # Further save a specific data frame containing a list of all significantly DEGs
      gene_names <- res_i@rownames
      res_i <- as.data.frame(res_i@listData) # Create a count matrix as a data frame
      rownames(res_i) <- gene_names
      res_i <- filter(res_i, padj < padjThreshold) # filter the count matrix for DEGs
      res_i <- filter(res_i, abs(log2FoldChange) > lfcThreshold)
      DEGs <- rownames(res_i)
      DEGs <- data.frame(DEGs, res_i$log2FoldChange, res_i$padj)
      colnames(DEGs) <- c("Gene_name", "log2FoldChange", "padj")
      # print(nrow(DEGs))
      
      # Save the object to file and name it according to its comparison details.
      print(sprintf("Saving %s_DEGs.csv to file", i))
      write.csv(DEGs, 
              file = sprintf("%s/DEG_lists/%s_DEGs.csv", mainDir, i))
      
    }
    
  }
  
  ##### Finally, perform visualzations for useful plots #####
  
  # Plot the total counts of normalised and unnormalised samples
  total_counts_unnorm <- data.frame(meta$donor_number, 
                                    meta$condition, 
                                    meta$treatment, 
                                    colSums(cts$counts))
  colnames(total_counts_unnorm) <- c("donor_number", "condition", "treatment", "total_counts")
  
  total_counts_norm <- data.frame(meta$donor_number, 
                                  meta$condition, 
                                  meta$treatment, 
                                  colSums(counts(dds, normalized=TRUE)))
  colnames(total_counts_norm) <- c("donor_number", "condition", "treatment", "total_counts")
  
  # Create plots
  p1 <- ggplot(total_counts_unnorm, aes(y = total_counts, x = donor_number, fill = treatment)) + 
    xlab("Donor number") + 
    ylab("Total counts") +
    geom_bar(stat = "identity") + 
    facet_grid(~ treatment) + 
    ggtitle("Unnormalised") + 
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("#46bac2", "#371ea3", "#b3eebe"))
  
  p2 <- ggplot(total_counts_norm, aes(y = total_counts, x = donor_number, fill = treatment)) + 
    xlab("Donor number") + 
    ylab("Total counts") +
    geom_bar(stat = "identity") + 
    facet_grid(~ treatment) + 
    ggtitle("Normalised") + 
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("#46bac2", "#371ea3", "#b3eebe"))
  
  cowplot::plot_grid(p1, p2, ncol = 2)
  
  # Save the plot to file
  ggsave(sprintf("%s/plots/DESeq2_plots/total_counts_barplot.png", mainDir), 
         dpi = 300, 
         width = 10, 
         height = 6, 
         units = "in")
  
  # Save the plot's data to file
  barplot_data <- list(total_counts_unnorm, total_counts_norm)
  saveRDS(barplot_data, file = sprintf("%s/plots/plotdata/total_counts_barplot.rds", mainDir))
  
  # Perform Variance Stabilizing Transformation (vst) for visualization of transformed data
  vsd <- vst(dds, blind=FALSE)
  
  pcaData <- plotPCA(vsd, intgroup=c("condition", "treatment"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  PCA_plot <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=treatment)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  PCA_plot
  
  # Save plot to file
  ggsave(filename = sprintf("%s/plots/DESeq2_plots/PCA.png", mainDir),  # name of the image file
         dpi = 300
  )
  
  # Save the plot's data to file
  saveRDS(PCA_plot, file = sprintf("%s/plots/plotdata/PCA.rds", mainDir))
  
  print("DESeq2 is complete.")
  print("------------------------------")
}
