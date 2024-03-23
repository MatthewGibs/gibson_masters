##### Plots Script #####
# This script contains functions to produce plots for the following sections of the data pipeline:
# - FastQC/MultiQC (Summary heatmap)
# - Salmon (Barplot of mapping rates)
# - DESeq2 DEGs (Venn Diagram)
# - DESeq2 DEGs which are common yet change L2FC sign.
# - Numerous clusterProfiler plots including: Enrichment maps (emapplot), Gene-concept networks (cnetplot), lollipop plots

##### MultiQC Heatmap #####
# mainDir <- "~/Matthew_Masters/Work" # set the path to the output files of running MultiQC
multiqc_heatmap <- function(mainDir # Path to the Work directory
                            )
{
  # Read in necessary packages
  multiqc_packages <- c("ggplot2", "dplyr", "tidyverse", "tidyr")
  suppressMessages(lapply(multiqc_packages, require, character.only = TRUE))
  
  # Read in TSV with MultiQC metrics
  Multiqc_data <- read_tsv(sprintf("%s/FastQC/MultiQC/multiqc_data/multiqc_fastqc.txt", mainDir))
  
  # Subset to keep metrics for heatmap and file names
  Multiqc_data <- Multiqc_data[,c(1, 12:21)]
  
  # Indicate the order for metrics to appear on the heatmap's x-axis
  metric_order <- c("Per Base Sequence Quality", 
                    "Per Tile Sequence Quality",
                    "Per Sequence Quality Scores",
                    "Per Base Sequence Content",
                    "Per Sequence GC Content",
                    "Per Base N Content",
                    "Sequence Length Distribution",
                    "Sequence Duplication Levels",
                    "Overrepresented Sequences",
                    "Adapter Content")
  
  # Transform the data
  colnames(Multiqc_data) <- c("Filename", metric_order)
  Multiqc_data$Filename <- factor(Multiqc_data$Filename)
  Multiqc_data <- gather(Multiqc_data, Metric, Value, 2:11)
  
  # Covert 'pass', 'warn', and 'fail' to representative numeric values
  Multiqc_data$Value <- as.numeric(gsub("pass", 1, gsub("warn", 0.5, gsub("fail", 0.25, Multiqc_data$Value))))
  Multiqc_data$Metric <- factor(Multiqc_data$Metric, levels = metric_order)
  
  # Plot the heatmap
  multiqc_heatmap_plot <- ggplot(Multiqc_data, aes(x = Metric, y = Filename, fill = Value)) + 
    geom_tile() + 
    theme(legend.position="none") + 
    scale_x_discrete("Quality Metric", labels = ~ str_wrap(as.character(metric_order), 15))
  multiqc_heatmap_plot
  
  # Write the plot to file
   ggsave(sprintf("%s/plots/multiqc_heatmap.png", mainDir), 
          bg = "white", 
          dpi = 300, 
          width = 12, 
          height = 8, 
          units = "in")
   
   saveRDS(multiqc_heatmap_plot, file = sprintf("%s/plots/plotdata/multiqc_heatmap.rds", mainDir))
}

##### Salmon Mapping Rates #####
# mainDir <- "~/Matthew_Masters/Work" # set the path to the output files of running Salmon
salmon_bar <- function(mainDir # Path to the 'Work' directory
) {
  # Load required packages.
  salmon_packages <- c("ggplot2", "dplyr", "tidyverse", "rjson")
  suppressMessages(lapply(salmon_packages, require, character.only = TRUE))
  
  setwd(mainDir)
  
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
  
  files <- file.path(mainDir, "salmon", quant_directories, "aux_info", "meta_info.json") # Set path to the quant.sf files of all listed samples
  
  processed_fragments <- c()
  mapped_fragments <- c()
  mapping_rates <- c()
  
  for(i in files) {
    result <- fromJSON(file = i)
    processed_fragments <- append(processed_fragments, result$num_processed)
    mapped_fragments <- append(mapped_fragments, result$num_mapped)
    mapping_rates <- append(mapping_rates, result$percent_mapped)
  }
  
  sample_names <- c("donor1_untreated",
                    "donor2_untreated",
                    "donor3_untreated",
                    "donor1_MCSF",
                    "donor2_MCSF",
                    "donor3_MCSF",
                    "donor1_GMCSF",
                    "donor2_GMCSF",
                    "donor3_GMCSF")
  
  # create a dataset
  sample <- c(rep(sample_names , 2))
  sample <-  factor(sample, level=c(sample_names))
  type <- c(rep(c("Unmapped" ) , 9), rep(c("Mapped" ) , 9))
  type <- factor(type, levels = c("Unmapped", "Mapped"))
  donor <- c(rep(c("One", "Two", "Three"), 6))
  donor <- factor(donor, levels = c("One", "Two", "Three"))
  treatment <- rep(c(rep("Untreated", 3), rep("M-CSF", 3), rep("GM-CSF", 3)), 2)
  treatment <- factor(treatment, levels = c("Untreated", "M-CSF", "GM-CSF"))
  Fragment_count <- mapped_fragments
  Fragment_count <- append(Fragment_count, x = (processed_fragments - mapped_fragments))
  unmapping_rates <- 100 - mapping_rates
  map_percentages <- append(unmapping_rates, mapping_rates)
  data <- data.frame(sample, donor, treatment, type, Fragment_count, map_percentages)
  
  # Stacked
  salmon_bar_plot <- ggplot(data, aes(fill = type, y = map_percentages, x = donor)) + 
    labs(fill = "Fragment Type") + 
    xlab("Donor number") + 
    ylab("Percentage (%)") + 
    geom_bar(position="stack", stat="identity") + 
    facet_grid(~ treatment) +
    scale_fill_manual(values=c("#371ea3", "#46bac2")) +
    theme(axis.title = element_text(size = 10), text = element_text(size = 10))
  
  salmon_bar_plot
  
  # Save plot to file
  ggsave(filename = sprintf("%s/plots/salmon_results_plot.png", mainDir),  # name of the image file
         dpi = 300,
         width = 2000,
         height = 1200,
         units = "px")
  
  # Save the R data to file
  saveRDS(salmon_bar_plot, file = sprintf("%s/plots/plotdata/salmon_results_plot.rds", mainDir))
  
}

##### DESeq2 DEG Venn Diagram #####
# recommended path for mainDir: mainDir <- "~/Matthew_Masters/Work"
DEG_venn <- function(mainDir # path to the 'Work' directory
) {
  # Read in necessary packages
  deseq_venn_packages <- c("ggplot2", "dplyr", "tidyverse", "ggVennDiagram")
  suppressMessages(lapply(deseq_venn_packages, require, character.only = TRUE))
  
  setwd(mainDir)
  
  # Read in the lists of DEGs
  DEGs <- readRDS("DEG_lists/GMCSF_and_MCSF_comparison.rds")
  
  GMCSF_up <- subset(DEGs$GMCSF, `log2FoldChange` > 0)
  GMCSF_down <- subset(DEGs$GMCSF, `log2FoldChange` < 0)
  MCSF_up <- subset(DEGs$MCSF, `log2FoldChange` > 0)
  MCSF_down <- subset(DEGs$MCSF, `log2FoldChange` < 0)
  
  # List of items
  x1 <- list(DEGs$GMCSF$Gene_name, DEGs$MCSF$Gene_name)
  x2 <- list(GMCSF_up$Gene_name, MCSF_up$Gene_name)
  x3 <- list(GMCSF_down$Gene_name, MCSF_down$Gene_name)
  
  # Plot 3 2D Venn diagrams
  p1 <- ggVennDiagram(x1, category.names = c("GM-CSF", "M-CSF")) + 
    theme(legend.position = "none") + ggtitle("All Differentially Expressed Genes") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p2 <- ggVennDiagram(x2, category.names = c("GM-CSF", "M-CSF"), label = "count") + 
    theme(legend.position = "none") + ggtitle("Upregulated DEGs") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p3 <- ggVennDiagram(x3, category.names = c("GM-CSF", "M-CSF"), label = "count") + 
    theme(legend.position = "none") + ggtitle("Downregulated DEGs") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p23 <- cowplot::plot_grid(p2, p3, ncol = 2)
  
  DEG_venn_plot <- cowplot::plot_grid(p1, p23, nrow = 2, rel_heights = c(3,2))
  
  ggsave(sprintf("%s/plots/DESeq2_plots/DEGs_venn.png", mainDir), 
         bg = "white",
         dpi = 300)
  
  saveRDS(DEG_venn_plot, sprintf("%s/plots/plotdata/DEGs_venn.rds", mainDir))
}

##### DESeq2 boxplots of interesting genes #####
# recommended path for mainDir: mainDir <- "~/Matthew_Masters/Work"
DEG_gene_plots <- function(mainDir # path to the 'Work' directory
) {
  # Read in necessary packages
  deseq_box_packages <- c("ggplot2", "DESeq2", "dplyr", "tidyverse", "data.table")
  suppressMessages(lapply(deseq_box_packages, require, character.only = TRUE))
  
  # Set the working directory
  setwd(mainDir)
  
  # Read in the two lists of DEGs
  GMCSF_DEGs <- read.csv("DEG_lists/treatment_GMCSF_vs_none_DEGs.csv", row.names = 1)
  MCSF_DEGs <- read.csv("DEG_lists/treatment_MCSF_vs_none_DEGs.csv", row.names = 1)
  
  # Identify the common DEGs
  common <- left_join(GMCSF_DEGs, MCSF_DEGs, by = "Gene_name")
  common <- na.omit(common)
  common <- common[,c(1,2,4)]
  colnames(common) <- c("Gene_name", "GMCSF_l2fc", "MCSF_l2fc")
  
  # Determine which DEGs have a different log2 fold change sign based on treatment
  f <- function(a, b) 
  {
    ifelse(a == 0 | b == 0, as.logical("FALSE"),!xor(sign(a)+1,sign(b)+1))
  }
  regulation <- f(a = common$GMCSF_l2fc, b = common$MCSF_l2fc)
  
  common <- data.frame(common$Gene_name, common$GMCSF_l2fc, common$MCSF_l2fc, regulation)
  
  # Extract each row containing a mismatch in lfc direction
  common <- common[common$regulation %like% FALSE, ]
  genes <- common$common.Gene_name
  
  # Annotate the list of interesting genes (ENSG ID to gene name)
  source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
  generate_ensembl(mainDir)
  annotate(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = genes # character vector containing the 'initial ID' type that you want to convert
  )
  gene_names <- object_final
  
  # Read in dds object with normalised counts for boxplots
  dds <- readRDS("DESeq2/dds.rds")
  
  # function to generate count data based on an ENSG ID
  gen_plot_data <- function(gene_number
  ){
    plotCounts(dds, 
               gene = gene_names[gene_number, 1], 
               intgroup="treatment", returnData=TRUE)
  }
  
  d1 <- gen_plot_data(1)
  d2 <- gen_plot_data(2)
  d3 <- gen_plot_data(3)
  d4 <- gen_plot_data(4)
  d5 <- gen_plot_data(5)
  d6 <- gen_plot_data(6)
  
  # function to create plots
  plot_DEGs <- function(DEG_data,
                        gene_number
  ) {
    ggplot(DEG_data, aes(x = treatment, y = count, fill = treatment)) + 
      geom_boxplot() + 
      theme(legend.position="none") +
      ggtitle(gene_names[gene_number, 2]) +
      theme(plot.title = element_text(hjust = 0.5)) + 
      scale_fill_brewer(palette="Dark2")
  }
  
  p1 <- plot_DEGs(d1, gene_number = 1) + xlab(NULL) 
  p2 <- plot_DEGs(d2, gene_number = 2) + ylab(NULL) + xlab(NULL)
  p3 <- plot_DEGs(d3, gene_number = 3) + ylab(NULL)
  p4 <- plot_DEGs(d4, gene_number = 4)
  p5 <- plot_DEGs(d5, gene_number = 5) + ylab(NULL)
  p6 <- plot_DEGs(d6, gene_number = 6) + ylab(NULL) + xlab(NULL)
  
  gene_plot <- cowplot::plot_grid(p1, p2, p6, p4, p5, p3, # Plots are ordered to better show trends
                     nrow = 2, ncol = 3)
  
  # Save the plot to file
  ggsave(sprintf("%s/plots/DESeq2_plots/common_DEGs_diff_l2fc.png", mainDir), 
         dpi = 300, 
         width = 10, 
         height = 6, 
         units = "in")
  
  saveRDS(gene_plot, sprintf("%s/plots/plotdata/common_DEGs_diff_l2fc.rds", mainDir))
}

##### clusterProfiler emapplot #####
# recommended path for mainDir: mainDir <- "~/Matthew_Masters/Work"
clusterProfiler_emapplot <- function(mainDir
) {
  setwd(mainDir)
  
  # Set and read in packages
  clusterProfiler_plots_packages <- c("data.table", "clusterProfiler", "ggplot2", "tidyverse", "dplyr", "org.Hs.eg.db", "enrichplot")
  suppressMessages(lapply(clusterProfiler_plots_packages, require, character.only = TRUE))
  
  # Read in enrichGO results for GM-CSF and M-CSF
  GMCSF <- readRDS("~/Matthew_Masters/Work/clusterProfiler_bulk/treatment_GMCSF_vs_none_enrichment/GO_BOTH_BP.rds")
  MCSF <- readRDS("~/Matthew_Masters/Work/clusterProfiler_bulk/treatment_MCSF_vs_none_enrichment/GO_BOTH_BP.rds")
  
  # Establish pairwise termism
  GMCSF <- pairwise_termsim(GMCSF)
  MCSF <- pairwise_termsim(MCSF)
  
  p1 <- emapplot(GMCSF, showcategory = 30)
  p2 <- emapplot(MCSF, showcategory = 30)
  
  cowplot::plot_grid(p1, p2, nrow = 2, labels = c("GM-CSF", "M-CSF"), label_size = 30)
  
  # Creates the output directory for this figure
  output_dir <- "plots/clusterProfiler_plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Creating 'clusterProfiler_plots' directory.")
  } else {}
  rm(output_dir)
  
  ggsave(sprintf("%s/plots/clusterProfiler_plots/emapplots.png", mainDir), 
         bg = "white",
         dpi = 400,
         width = 13,
         height = 20,
         units = "in")
}

##### clusterProfiler cnetplots #####
# recommended path for mainDir: mainDir <- "~/Matthew_Masters/Work"
clusterProfiler_cnetplot <- function(mainDir, # The path to the 'Work' directory
                                     categories = c("regulation of inflammatory response") # A character vector containing GO terms to plot
) {
  # Set and read in packages
  clusterProfiler_plots_packages <- c("data.table", "clusterProfiler", "ggplot2", "tidyverse", "dplyr", "org.Hs.eg.db", "enrichplot")
  suppressMessages(lapply(clusterProfiler_plots_packages, require, character.only = TRUE))
  
  setwd(mainDir)
  
  # Read in DEGs
  GMCSF_DEGs <- read.csv("DEG_lists/treatment_GMCSF_vs_none_DEGs.csv", row.names = 1)
  MCSF_DEGs <- read.csv("DEG_lists/treatment_MCSF_vs_none_DEGs.csv", row.names = 1)
  
  source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
  generate_ensembl(mainDir)
  
  annotate(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = GMCSF_DEGs$Gene_name) # character vector containing the 'initial ID' type that you want to convert
  GMCSF_DEGs$Gene_name <- object_final[,2]
  GMCSF_geneList <- GMCSF_DEGs$log2FoldChange
  names(GMCSF_geneList) <- GMCSF_DEGs$Gene_name
  
  annotate(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = MCSF_DEGs$Gene_name) # character vector containing the 'initial ID' type that you want to convert
  MCSF_DEGs$Gene_name <- object_final[,2]
  MCSF_geneList <- MCSF_DEGs$log2FoldChange
  names(MCSF_geneList) <- MCSF_DEGs$Gene_name
  
  # Read in enrichGO results
  GMCSF_GO <- readRDS("clusterProfiler/GO_enrichment/GMCSF/GO_BOTH.rds")
  MCSF_GO <- readRDS("clusterProfiler/GO_enrichment/MCSF/GO_BOTH.rds")
  
  p1 <- cnetplot(GMCSF_GO, color.params = list(foldChange = GMCSF_geneList, category = 'firebrick'), showCategory = categories)
  p2 <- cnetplot(MCSF_GO, color.params = list(foldChange = MCSF_geneList, category = 'firebrick'), showCategory = categories)
  cowplot::plot_grid(p1, p2, ncol=2, labels = c("GM-CSF", "M-CSF"), label_size = 30)
  
  # Creates the output directory for this script
  output_dir <- "plots/clusterProfiler_plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Creating 'clusterProfiler_plots' directory.")
  } else {}
  rm(output_dir)
  
  ggsave(sprintf("%s/plots/clusterProfiler_plots/cnetplots.png", mainDir), 
         bg = "white",
         dpi = 300,
         width = 16,
         height = 8,
         units = "in")
}

##### clusterProfiler lollipop plots (Unique genes GO terms) #####
# recommended path for mainDir: mainDir <- "~/Matthew_Masters/Work"
clusterProfiler_loliplot <- function(mainDir
) {
  # Set and read in packages
  clusterProfiler_plots_packages <- c("data.table", "clusterProfiler", "ggplot2", "tidyverse", "dplyr", "org.Hs.eg.db", "enrichplot", "forcats")
  suppressMessages(lapply(clusterProfiler_plots_packages, require, character.only = TRUE))
  
  setwd(mainDir)
  
  # Read in the data
  GMCSF_unique <- readRDS("clusterProfiler_bulk/GMCSF_unique_enrichment/GO_BOTH_BP.rds")
  MCSF_unique <- readRDS("clusterProfiler_bulk/MCSF_unique_enrichment/GO_BOTH_BP.rds")
  
  p1 <- ggplot(GMCSF_unique, showCategory = 15, aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
    geom_segment(aes(xend = 0, yend = Description)) +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_viridis_c(guide = guide_colorbar(reverse=TRUE)) +
    scale_size_continuous(range =c (1, 7)) +
    theme_minimal() + 
    xlab("Gene Ratio") +
    ylab(NULL) + 
    ggtitle("GM-CSF")
  
  p2 <- ggplot(MCSF_unique, showCategory = 15, aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
    geom_segment(aes(xend = 0, yend = Description)) +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_viridis_c(guide = guide_colorbar(reverse=TRUE)) +
    scale_size_continuous(range =c (1, 7)) +
    theme_minimal() + 
    xlab("Gene Ratio") +
    ylab(NULL) + 
    ggtitle("M-CSF")
  
  cowplot::plot_grid(p1, p2, ncol = 2)
  
  ggsave(sprintf("%s/plots/clusterProfiler_plots/unique_GOs.png", mainDir), 
         bg = "white",
         dpi = 300,
         width = 14,
         height = 6,
         units = "in")
}



##### Pathview combined plot #####
# recommended path for mainDir: mainDir <- "~/Matthew_Masters/Work"
KEGG_compare_pathways <- function(mainDir,
                                  specific_ID) 
  {
  setwd(mainDir)
  
  KEGG_packages <- c("clusterProfiler", 
                     "AnnotationDbi", 
                     "org.Hs.eg.db",
                     "R.utils", 
                     "biomaRt",
                     "dplyr",
                     "pathview")
  
  suppressMessages(lapply(KEGG_packages, require, character.only = TRUE))
  
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
  output_dir <- "plots/KEGG_plots/combined_plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'KEGG_plots/combined_plots' directory. Figures generated from this pathview run can be found there.")
  } else {
    print("Figures generated this pathview can be found in the 'plots/KEGG_plots/combined_plots' directory.")
  }
  
  rm(output_dir)
  
  print("All the necessary directories have been created or already exist, proceeding with KEGG pathway visualisation.")
  print("------------------------------")
  
  # Read in data
  
  GMCSF_genes <- read.csv("DEG_lists/treatment_GMCSF_vs_none_DEGs.csv", row.names = 1)
  MCSF_genes <- read.csv("DEG_lists/treatment_MCSF_vs_none_DEGs.csv", row.names = 1)
  
  GMCSF_entrez <- mapIds(org.Hs.eg.db,
                         keys = GMCSF_genes$Gene_name, #Column containing Ensembl gene ids
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  GMCSF_genes$Gene_name <- GMCSF_entrez
  GMCSF_genes <- na.omit(GMCSF_genes)
  
  MCSF_entrez <- mapIds(org.Hs.eg.db,
                        keys = MCSF_genes$Gene_name, #Column containing Ensembl gene ids
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")
  MCSF_genes$Gene_name <- MCSF_entrez
  MCSF_genes <- na.omit(MCSF_genes)
  
  combined <- full_join(GMCSF_genes, MCSF_genes, by = "Gene_name")
  gene_data <- data.frame(combined[,2], combined[,4])
  rownames(gene_data) <- combined[,1]
  colnames(gene_data) <- c("GMCSF", "MCSF")
  
  # Run pathview
  # Create the sub-directory for the specified pathway if necessary.
  output_dir <- sprintf("plots/KEGG_plots/combined_plots/%s", specific_ID)
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  # As the pathview function creates the image in the current directory, we must switch to where it must be saved.
  setwd(sprintf("%s/plots/KEGG_plots/combined_plots/%s", mainDir, specific_ID))
  
  hsa <- pathview(gene.data = gene_data, 
                  pathway.id = specific_ID, 
                  species = "hsa",
                  keys.align = "y",
                  kegg.native = T,
                  match.data = F, 
                  multi.state = T, 
                  same.layer = T,
                  limit = list(gene = 10, cpd = 1), 
                  low = list(gene = "steelblue4"),
                  mid = list(gene = "grey90"),
                  high = list(gene = "coral"))
  
  saveRDS(hsa, file = sprintf("%s/plots/KEGG_plots/combined_plots/%s/%s_pathview.rds", mainDir, specific_ID, specific_ID))
  
  setwd(mainDir)
}

##### Plot x #####
# put function here{}

##### Plot x #####
# put function here{}

##### Plot x #####
# put function here{}

