##### Plots Script #####
# This script contains functions to produce plots for the following sections of the data pipeline:
# - FastQC/MultiQC (Summary heatmap)
# - Salmon (Barplot of mapping rates)
# - DESeq2 DEGs (Venn Diagram)
# - DESeq2 DEGs which are common yet change L2FC sign.
# - Numerous clusterProfiler plots including: Enrichment maps (emapplot), Gene-concept networks (cnetplot), lollipop plots

##### Generate all plots #####
# This function is capable of calling all others within this script, providing a quick alternative for generating all figures.
gen_all_plots <- function(mainDir) {
  multiqc_heatmap(mainDir)
  
  salmon_bar(mainDir)
  
  marker_gene_plots(mainDir)
  
  circular_marker_heatmap(mainDir)
  
  DEG_gene_plots(mainDir)
  
  DEG_venn(mainDir)
  
  GO_venn(mainDir)
}

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
  venn_packages <- c("ggplot2", "dplyr", "tidyverse", "VennDiagram", "hrbrthemes", "gridExtra", "png", "cowplot")
  suppressMessages(lapply(venn_packages, require, character.only = TRUE))
  
  setwd(mainDir)
  
  # Create necessary directories
  # Create the 'plots' directory if necessary.
  output_dir <- "plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created 'plots' directory.")
  }
  
  # Create the 'KEGG_plots' directory if necessary.
  output_dir <- "plots/Venn_plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'Venn_plots' directory. Venn diagrams can be found there.")
  } else {
    print("Venn diagrams can be found in the 'plots/Venn_plots' directory.")
  }
  
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
  
  venn_plot <- function(list_data,
                        file_name,
                        plot_title) 
  {
    venn.diagram(x = list_data,
                 filename = file_name,
                 category.names = c("", ""),
                 main = plot_title,
                 output = TRUE ,
                 imagetype = "png" ,
                 height = 2000, 
                 width = 2000, 
                 resolution = 1000,
                 compression = "lzw",
                 lwd = 1,
                 col = c('coral', "steelblue3"),
                 fill = c(alpha("coral",0.3), alpha("steelblue3",0.3)),
                 disable.logging = TRUE,
                 cex = 0.5,
                 fontfamily = "sans",
                 cat.cex = 0.3,
                 cat.default.pos = "outer",
                 cat.pos = c(-27, 27),
                 cat.dist = c(0.055, 0.055),
                 cat.fontfamily = "sans")
  }
  
  # change to the plot directory
  setwd("plots/Venn_plots")
  
  # Make the plots
  venn_plot(list_data = x1, file_name = "all_DEGs_venn.png", plot_title = "A")
  
  venn_plot(list_data = x2, file_name = "up_DEGs_venn.png", plot_title = "B")
  
  venn_plot(list_data = x3, file_name = "down_DEGs_venn.png", plot_title = "C")
  
  # Import PNG files
  img1 <- readPNG("all_DEGs_venn.png")
  img2 <- readPNG("up_DEGs_venn.png")
  img3 <- readPNG("down_DEGs_venn.png")
  
  # Create a legend
  legend_labels <- c("GM-CSF", "M-CSF")
  legend_colors <- c("coral", "steelblue3")
  
  # Create a blank ggplot with the legend
  dummy_legend_plot <- ggplot() +
    geom_point(aes(x = 0, y = 0, color = factor(legend_labels))) +
    scale_color_manual(values = legend_colors, name = "Treatment") +  # Set the legend title
    theme_void() +
    theme(legend.position = "bottom",  # Change to "bottom" for a horizontal legend
          legend.direction = "vertical",  # Change to "horizontal" for a horizontal legend
          legend.box = "vertical")  # Place legend labels under one another
  
  # Extract the legend
  legend_grob <- cowplot::get_legend(dummy_legend_plot)
  
  # Create rasterGrob objects for each image
  grob_img1 <- rasterGrob(img1, interpolate = TRUE)
  grob_img2 <- rasterGrob(img2, interpolate = TRUE)
  grob_img3 <- rasterGrob(img3, interpolate = TRUE)
  
  # Arrange the plots and legend into a composite image
  composite_plot <- grid.arrange(
    grob_img1, 
    grid.arrange(grob_img2, grob_img3, ncol = 1),  # Stack smaller Venn diagrams
    legend_grob,
    ncol = 3,
    widths = c(1.5, 1, 0.5)  # Adjust the widths as needed
  )
  
  # Save the composite image with legend
  png("DEG_venns", width = 12, height = 4, units = "in", res = 400)
  cowplot::save_plot("DEG_venns.png", composite_plot)
  dev.off()  # Close the graphics device
  
  # return to the working directory
  setwd(mainDir)
}

##### DESeq2 boxplots of marker genes #####
# recommended path for mainDir: mainDir <- "~/Matthew_Masters/Work"
marker_gene_plots <- function(mainDir # path to the 'Work' directory
) {
  # Read in necessary packages
  deseq_box_packages <- c("ggplot2", "DESeq2", "dplyr", "tidyverse", "data.table", "ggpubr")
  suppressMessages(lapply(deseq_box_packages, require, character.only = TRUE))
  
  # Set the working directory
  setwd(mainDir)
  
  # Read in lists of DEGs
  DEGs <- readRDS("DEG_lists/GMCSF_and_MCSF_comparison.rds")
  
  # Convert ENSG ID to gene name
  source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
  generate_ensembl(mainDir)
  
  # Manually provde the set of genes (M0 markers) to visualize.
  gene_names <- c(#"ENSG00000174837", # ADGRE1
                  #"ENSG00000121807", # CCR2
                  #"ENSG00000170458", # CD14
                  #"ENSG00000129226", # CD68
                  #"ENSG00000177575", # CD163
                  "ENSG00000182578", # CSF1R
                  #"ENSG00000168329", # CX3CR1
                  "ENSG00000203747", # FCGR3A
                  #"ENSG00000169896", # ITGAM
                  "ENSG00000260314", # MRC1
                  #"ENSG00000088827", # SIGLEC1
                  "ENSG00000118785"  # SPP1
  )
  
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "external_gene_name", # The ID type you want in your final object.
             object = gene_names # character vector containing the 'initial ID' type that you want to convert
  )
  gene_names <- data.frame(object_final)
  colnames(gene_names) <- c("Gene_name", "Gene_symbol")
  
  GMDEGs <- left_join(gene_names, DEGs$GMCSF, by = "Gene_name")
  GMDEGs <- na.omit(GMDEGs)
  
  MDEGs <- left_join(gene_names, DEGs$MCSF, by = "Gene_name")
  MDEGs <- na.omit(MDEGs)
  
  # Read in dds object with normalised counts for boxplots
  dds <- readRDS("DESeq2/dds.rds")
  
  # function to generate count data based on an ENSG ID
  gen_plot_data <- function(gene_number
  ){
    plotCounts(dds, 
               gene = as.character(gene_names[gene_number, 1]), 
               intgroup= "treatment", returnData = TRUE)
  }
  
  d1 <- gen_plot_data(1)
  d2 <- gen_plot_data(2)
  d3 <- gen_plot_data(3)
  d4 <- gen_plot_data(4)
  #d5 <- gen_plot_data(5)
  #d6 <- gen_plot_data(6)
  #d7 <- gen_plot_data(7)
  #d8 <- gen_plot_data(8)
  #d9 <- gen_plot_data(9)
  #d10 <- gen_plot_data(10)
  #d11 <- gen_plot_data(11)
  #d12 <- gen_plot_data(12)
  
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
  
  p1 <- plot_DEGs(d1, gene_number = 1)
  p2 <- plot_DEGs(d2, gene_number = 2) + ylab(NULL) 
  p3 <- plot_DEGs(d3, gene_number = 3) + xlab(NULL)
  p4 <- plot_DEGs(d4, gene_number = 4) + ylab(NULL) + xlab(NULL)
  #p5 <- plot_DEGs(d5, gene_number = 5) + ylab(NULL)
  #p6 <- plot_DEGs(d6, gene_number = 6) + xlab(NULL)
  #p7 <- plot_DEGs(d7, gene_number = 7) + ylab(NULL) + xlab(NULL)
  #p8 <- plot_DEGs(d8, gene_number = 8) + ylab(NULL)
  #p9 <- plot_DEGs(d9, gene_number = 9) + ylab(NULL) + xlab(NULL)
  #p10 <- plot_DEGs(d10, gene_number = 10) + ylab(NULL) + xlab(NULL)
  #p11 <- plot_DEGs(d11, gene_number = 11)
  #p12 <- plot_DEGs(d12, gene_number = 12) + xlab(NULL)
  
  #cowplot::plot_grid(p12, p2, p4, p9, p6, p10, p1, p7, p11, p3, p5, p8, # Plots are ordered to better show trends
  #                   nrow = 3, ncol = 4)
  cowplot::plot_grid(p3, p4, p1, p2, # Plots are ordered to better show trend
                     nrow = 2, ncol = 2)
  
  # Save the plot to file
  ggsave(sprintf("%s/plots/surface_markers.png", mainDir), 
         dpi = 300, 
         width = 12, 
         height = 8, 
         units = "in")
}

##### Circular heatmap of macrophage subtype marker genes #####
# recommended path for mainDir: mainDir <- "~/Matthew_Masters/Work"
circular_marker_heatmap <- function(mainDir # path to the 'Work' directory
) {
  # Read in necessary packages
  deseq_box_packages <- c("ggplot2", "DESeq2", "dplyr", "tidyverse", "data.table", "ggpattern", "ComplexHeatmap", "circlize")
  suppressMessages(lapply(deseq_box_packages, require, character.only = TRUE))
  
  # Set the working directory
  setwd(mainDir)
  
  # Read in lists of DEGs
  DEGs <- readRDS("DEG_lists/GMCSF_and_MCSF_comparison.rds")
  
  # Convert ENSG ID to gene name
  source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
  generate_ensembl(mainDir)
  
  # M1 markers
  M1_markers <- c("ENSG00000121594", # CD80
                  "ENSG00000114013", # CD86 = B7-2
                  #"ENSG00000172243", # CLEC7A = Dectin-1/CD369
                  "ENSG00000164400", # CSF2 = GM-CSF
                  "ENSG00000081041", # CXCL2 = MIP-2a
                  #"ENSG00000150337", # FCGR1A = CD64
                  #"ENSG00000143226", # FCGR2A = CD32
                  #"ENSG00000072694", # FCGR2B = CD32
                  #"ENSG00000244682", # FCGR2C = CD32
                  #"ENSG00000204287", # HLA-DRA
                  #"ENSG00000111537", # IFNG = IFN gamma
                  #"ENSG00000027697", # IFNGR1
                  "ENSG00000125538", # IL1B
                  "ENSG00000115594", # IL1R1
                  "ENSG00000136244", # IL6
                  #"ENSG00000113302", # IL12B
                  "ENSG00000110944", # IL23A
                  "ENSG00000128604", # IRF5
                  #"ENSG00000019169", # MARCO
                  "ENSG00000109320", # NFKB1
                  "ENSG00000077150", # NFKB2
                  "ENSG00000007171", # NOS2
                  "ENSG00000197329", # PELI1
                  "ENSG00000115415", # STAT1
                  #"ENSG00000137462", # TLR2 = CD282
                  #"ENSG00000136869", # TLR4 = CD284
                  "ENSG00000232810" # TNF = TNF alpha
  )
  
  # M2 markers
  M2_markers <- c("ENSG00000118520", # ARG1
                  "ENSG00000102962", # CCL22
                  "ENSG00000177575", # CD163
                  "ENSG00000090659", # CD209
                  "ENSG00000132514", # CLEC10A = mgl2/CD301
                  #"ENSG00000182578", # CSF1R = CD115
                  #"ENSG00000179639", # FCER1A
                  #"ENSG00000203747", # FCGR3A = CD16
                  "ENSG00000165457", # FOLR2
                  "ENSG00000017427", # IGF1
                  #"ENSG00000113520", # IL4
                  "ENSG00000077238", # IL4R
                  "ENSG00000136634", # IL10
                  "ENSG00000137265", # IRF4 = MUM1
                  #"ENSG00000260314", # MRC1 = CD206/mannose receptor
                  "ENSG00000038945", # MSR1
                  "ENSG00000197646", # PDCD1LG2 = pdl2
                  #"ENSG00000100311", # PDGFB
                  "ENSG00000112033", # PPARD
                  "ENSG00000132170", # PPARG
                  "ENSG00000166888" # STAT6
                  #"ENSG00000105329", # TGFB1
                  #"ENSG00000095970", # TREM2
                  #"ENSG00000155659" # VSIG4
  )
  
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "external_gene_name", # The ID type you want in your final object.
             object = M1_markers # character vector containing the 'initial ID' type that you want to convert
  )
  M1_markers <- data.frame(object_final)
  colnames(M1_markers) <- c("Gene_name", "Gene_symbol")
  
  GMCSF_M1 <- left_join(M1_markers, DEGs$GMCSF, by = "Gene_name")
  MCSF_M1 <- left_join(M1_markers, DEGs$MCSF, by = "Gene_name")
  
  # Create a vector of gene symbols which are differentially expressed in at least 1 treatment.
  M1_DEGs <- rbind(na.omit(GMCSF_M1), na.omit(MCSF_M1))
  M1_DEGs <- data.frame(Gene_symbol = unique(M1_DEGs$Gene_symbol))
  
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "external_gene_name", # The ID type you want in your final object.
             object = M2_markers # character vector containing the 'initial ID' type that you want to convert
  )
  
  M2_markers <- data.frame(object_final)
  colnames(M2_markers) <- c("Gene_name", "Gene_symbol")
  
  GMCSF_M2 <- left_join(M2_markers, DEGs$GMCSF, by = "Gene_name")
  MCSF_M2 <- left_join(M2_markers, DEGs$MCSF, by = "Gene_name")
  
  # Create a vector of gene symbols which are differentially expressed in at least 1 treatment.
  M2_DEGs <- rbind(na.omit(GMCSF_M2), na.omit(MCSF_M2))
  M2_DEGs <- data.frame(Gene_symbol = unique(M2_DEGs$Gene_symbol))
  
  # Read in dds object with normalised counts for calculating log2 fold changes manually
  dds <- readRDS("DESeq2/dds.rds")
  
  # Manually calculate and replace log2 fold changes for the M1 markers.
  for (i in seq_len(nrow(M1_markers))) {
    print(M1_markers[i, 1])
    
    gene_id <- as.character(M1_markers[i, 1])
    
    # Check if the gene is present in the count matrix
    if (gene_id %in% rownames(assay(dds))) {
      gene_counts <- plotCounts(dds, gene = gene_id, intgroup = "treatment", returnData = TRUE)
      
      GMCSF_M1[i, 3] <- log2((gene_counts[7, 1] + gene_counts[8, 1] + gene_counts[9, 1]) / (gene_counts[1, 1] + gene_counts[2, 1] + gene_counts[3, 1]))
      MCSF_M1[i, 3] <- log2((gene_counts[4, 1] + gene_counts[5, 1] + gene_counts[6, 1]) / (gene_counts[1, 1] + gene_counts[2, 1] + gene_counts[3, 1]))
    } else {
      # Gene not present, skip and assign NA to log2 fold changes
      GMCSF_M1[i, 3] <- 0
      MCSF_M1[i, 3] <- 0
      print("Gene not found in count matrix, skipping.")
    }
  }
  
  GMCSF_M1_notDEGs <- data.frame(Gene_symbol = setdiff(GMCSF_M1$Gene_symbol, M1_DEGs$Gene_symbol))
  GMCSF_M1_notDEGs <- left_join(GMCSF_M1_notDEGs, GMCSF_M1, by = "Gene_symbol")
  GMCSF_M1_DEGs <- left_join(M1_DEGs, GMCSF_M1, by = "Gene_symbol")
  
  MCSF_M1_notDEGs <- data.frame(Gene_symbol = setdiff(MCSF_M1$Gene_symbol, M1_DEGs$Gene_symbol))
  MCSF_M1_notDEGs <- left_join(MCSF_M1_notDEGs, MCSF_M1, by = "Gene_symbol")
  MCSF_M1_DEGs <- left_join(M1_DEGs, MCSF_M1, by = "Gene_symbol")
  
  
  # Manually calculate and replace log2 fold changes for the M2 markers.
  for (i in seq_len(nrow(M2_markers))) {
    print(M2_markers[i, 1])
    
    gene_id <- as.character(M2_markers[i, 1])
    
    # Check if the gene is present in the count matrix
    if (gene_id %in% rownames(assay(dds))) {
      gene_counts <- plotCounts(dds, gene = gene_id, intgroup = "treatment", returnData = TRUE)
      
      GMCSF_M2[i, 3] <- log2((gene_counts[7, 1] + gene_counts[8, 1] + gene_counts[9, 1]) / (gene_counts[1, 1] + gene_counts[2, 1] + gene_counts[3, 1]))
      MCSF_M2[i, 3] <- log2((gene_counts[4, 1] + gene_counts[5, 1] + gene_counts[6, 1]) / (gene_counts[1, 1] + gene_counts[2, 1] + gene_counts[3, 1]))
    } else {
      # Gene not present, skip and assign 0 to log2 fold changes
      GMCSF_M2[i, 3] <- 0
      MCSF_M2[i, 3] <- 0
      print("Gene not found in count matrix, skipping.")
    }
  }
  
  GMCSF_M2_notDEGs <- data.frame(Gene_symbol = setdiff(GMCSF_M2$Gene_symbol, M2_DEGs$Gene_symbol))
  GMCSF_M2_notDEGs <- left_join(GMCSF_M2_notDEGs, GMCSF_M2, by = "Gene_symbol")
  GMCSF_M2_DEGs <- left_join(M2_DEGs, GMCSF_M2, by = "Gene_symbol")
  
  MCSF_M2_notDEGs <- data.frame(Gene_symbol = setdiff(MCSF_M2$Gene_symbol, M2_DEGs$Gene_symbol))
  MCSF_M2_notDEGs <- left_join(MCSF_M2_notDEGs, MCSF_M2, by = "Gene_symbol")
  MCSF_M2_DEGs <- left_join(M2_DEGs, MCSF_M2, by = "Gene_symbol")
  
  # Create numeric matrices for each group of differentially expressed genes
  mat_GMCSF_M1_DE <- as.matrix(GMCSF_M1_DEGs[, c("log2FoldChange")])
  rownames(mat_GMCSF_M1_DE) <- GMCSF_M1_DEGs$Gene_symbol
  
  mat_GMCSF_M2_DE <- as.matrix(GMCSF_M2_DEGs[, c("log2FoldChange")])
  rownames(mat_GMCSF_M2_DE) <- GMCSF_M2_DEGs$Gene_symbol
  
  mat_MCSF_M1_DE <- as.matrix(MCSF_M1_DEGs[, c("log2FoldChange")])
  rownames(mat_MCSF_M1_DE) <- MCSF_M1_DEGs$Gene_symbol
  
  mat_MCSF_M2_DE <- as.matrix(MCSF_M2_DEGs[, c("log2FoldChange")])
  rownames(mat_MCSF_M2_DE) <- MCSF_M2_DEGs$Gene_symbol
  
  # Create numeric matrices for each group of non-differentially expressed genes
  mat_GMCSF_M1_notDE <- as.matrix(GMCSF_M1_notDEGs[, c("log2FoldChange")])
  rownames(mat_GMCSF_M1_notDE) <- GMCSF_M1_notDEGs$Gene_symbol
  
  mat_GMCSF_M2_notDE <- as.matrix(GMCSF_M2_notDEGs[, c("log2FoldChange")])
  rownames(mat_GMCSF_M2_notDE) <- GMCSF_M2_notDEGs$Gene_symbol
  
  mat_MCSF_M1_notDE <- as.matrix(MCSF_M1_notDEGs[, c("log2FoldChange")])
  rownames(mat_MCSF_M1_notDE) <- MCSF_M1_notDEGs$Gene_symbol
  
  mat_MCSF_M2_notDE <- as.matrix(MCSF_M2_notDEGs[, c("log2FoldChange")])
  rownames(mat_MCSF_M2_notDE) <- MCSF_M2_notDEGs$Gene_symbol
  
  # Combine matrices
  heat_matrix <- rbind(rbind(cbind(mat_MCSF_M1_notDE, mat_GMCSF_M1_notDE), cbind(mat_MCSF_M1_DE, mat_GMCSF_M1_DE)), 
                       rbind(cbind(mat_MCSF_M2_DE, mat_GMCSF_M2_DE), cbind(mat_MCSF_M2_notDE, mat_GMCSF_M2_notDE)))
  
  colnames(heat_matrix) <- c("MCSF", "GMCSF")
  
  # Create the split factor for the heatmap
  split <- c(rep("M1 not DE", length(mat_GMCSF_M1_notDE)), rep("M1 DE", length(mat_GMCSF_M1_DE)), 
             rep("M2 DE", length(mat_GMCSF_M2_DE)), rep("M2 not DE", length(mat_GMCSF_M2_notDE)))
  split <- factor(split, levels = c("M1 not DE", "M1 DE", "M2 DE", "M2 not DE"))
  
  circos.clear()
  png("plots/subtype_heatmap.png", width = 5000, height = 5000, res = 600)
  col_fun1 = colorRamp2(c(-5, 0, 5), c("steelblue3", "grey95", "coral"))
  circos.par(start.degree = 90, gap.degree = 10)
  circos.heatmap(heat_matrix, 
                 split = split, 
                 col = col_fun1, 
                 track.height = 0.25, 
                 bg.lwd = 2, 
                 bg.lty = 2, 
                 show.sector.labels = TRUE, 
                 rownames.side = "inside",
                 rownames.cex = 0.75)
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if(CELL_META$sector.numeric.index == 4) { # the last sector
      cn = c("GM-CSF", "M-CSF")
      n = length(cn)
      circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(2, "mm"), 
                  1:n - 0.5, cn, 
                  cex = 0.5, adj = c(0, 0.5), facing = "inside")
    }
  }, bg.border = NA)
  lgd = Legend(title = "Log2 Fold Change", col_fun = col_fun1, title_position = c("topcenter"))
  grid.draw(lgd)
  dev.off()
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
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
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
               gene = as.character(gene_names[gene_number, 1]), 
               intgroup="treatment", returnData=TRUE)
  }
  
  d1 <- gen_plot_data(1)
  d2 <- gen_plot_data(2)
  d3 <- gen_plot_data(3)
  d4 <- gen_plot_data(4)
  d5 <- gen_plot_data(5)
  d6 <- gen_plot_data(6)
  d7 <- gen_plot_data(7)
  d8 <- gen_plot_data(8)
  d9 <- gen_plot_data(9)
  d10 <- gen_plot_data(10)
  d11 <- gen_plot_data(11)
  d12 <- gen_plot_data(12)
  d13 <- gen_plot_data(13)
  d14 <- gen_plot_data(14)
  d15 <- gen_plot_data(15)
  d16 <- gen_plot_data(16)
  d17 <- gen_plot_data(17)
  d18 <- gen_plot_data(18)
  d19 <- gen_plot_data(19)
  d20 <- gen_plot_data(20)
  d21 <- gen_plot_data(21)
  d22 <- gen_plot_data(22)
  d23 <- gen_plot_data(23)
  d24 <- gen_plot_data(24)
  d25 <- gen_plot_data(25)
  
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
  p7 <- plot_DEGs(d7, gene_number = 7) + ylab(NULL) + xlab(NULL)
  p8 <- plot_DEGs(d8, gene_number = 8) + ylab(NULL) + xlab(NULL)
  p9 <- plot_DEGs(d9, gene_number = 9) + ylab(NULL) + xlab(NULL)
  p10 <- plot_DEGs(d10, gene_number = 10) + ylab(NULL) + xlab(NULL)
  p11 <- plot_DEGs(d11, gene_number = 11) + ylab(NULL) + xlab(NULL)
  p12 <- plot_DEGs(d12, gene_number = 12) + ylab(NULL) + xlab(NULL)
  p13 <- plot_DEGs(d13, gene_number = 13) + ylab(NULL) + xlab(NULL)
  p14 <- plot_DEGs(d14, gene_number = 14) + ylab(NULL) + xlab(NULL)
  p15 <- plot_DEGs(d15, gene_number = 15) + ylab(NULL) + xlab(NULL)
  p16 <- plot_DEGs(d16, gene_number = 16) + ylab(NULL) + xlab(NULL)
  p17 <- plot_DEGs(d17, gene_number = 17) + ylab(NULL) + xlab(NULL)
  p18 <- plot_DEGs(d18, gene_number = 18) + ylab(NULL) + xlab(NULL)
  p19 <- plot_DEGs(d19, gene_number = 19) + ylab(NULL) + xlab(NULL)
  p20 <- plot_DEGs(d20, gene_number = 20) + ylab(NULL) + xlab(NULL)
  p21 <- plot_DEGs(d21, gene_number = 21) + ylab(NULL) + xlab(NULL)
  p22 <- plot_DEGs(d22, gene_number = 22) + ylab(NULL) + xlab(NULL)
  p23 <- plot_DEGs(d23, gene_number = 23) + ylab(NULL) + xlab(NULL)
  p24 <- plot_DEGs(d24, gene_number = 24) + ylab(NULL) + xlab(NULL)
  p25 <- plot_DEGs(d25, gene_number = 25) + ylab(NULL) + xlab(NULL)
  
  gene_plot <- cowplot::plot_grid(p1, p2, p3, p4, p5, 
                                  p6, p7, p8, p9, p10, 
                                  p11, p12, p13, p14, p15,
                                  p16, p17, p18, p19, p20,
                                  p21, p22, p23, p24, p25, # Plots are ordered to better show trends
                                  nrow = 5, ncol = 5)
  
  # Save the plot to file
  ggsave(sprintf("%s/plots/DESeq2_plots/common_DEGs_diff_l2fc.png", mainDir), 
         dpi = 300, 
         width = 10, 
         height = 10, 
         units = "in")
  
  saveRDS(gene_plot, sprintf("%s/plots/plotdata/common_DEGs_diff_l2fc.rds", mainDir))
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
  
  # Create all the necessary directories for outputs
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
  
  # Keep only the first occurrence of each row name
  unique_row_names <- !duplicated(combined$Gene_name)
  gene_data <- data.frame(combined[unique_row_names, 2], combined[unique_row_names, 4])
  
  # Set unique row names for the gene_data data frame
  rownames(gene_data) <- combined$Gene_name[unique_row_names]
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

##### GO enrichment venn diagram #####
GO_venn <- function(mainDir) {
  # Read in necessary packages
  venn_packages <- c("ggplot2", "dplyr", "tidyverse", "VennDiagram", "hrbrthemes", "gridExtra", "png", "cowplot")
  suppressMessages(lapply(venn_packages, require, character.only = TRUE))
  
  setwd(mainDir)
  
  # Create necessary directories
  # Create the 'plots' directory if necessary.
  output_dir <- "plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created 'plots' directory.")
  }
  
  # Create the 'KEGG_plots' directory if necessary.
  output_dir <- "plots/Venn_plots"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created the 'Venn_plots' directory. Venn diagrams can be found there.")
  } else {
    print("Venn diagrams can be found in the 'plots/Venn_plots' directory.")
  }
  
  GO_results <- readRDS("clusterProfiler/GO_enrichment/GO_summary.rds")
  
  # List of items
  x1 <- list(GO_results$GMCSF$GMCSF_both$ID, GO_results$MCSF$MCSF_both$ID)
  x2 <- list(GO_results$GMCSF$GMCSF_up$ID, GO_results$MCSF$MCSF_up$ID)
  x3 <- list(GO_results$GMCSF$GMCSF_down$ID, GO_results$MCSF$MCSF_down$ID)
  x4 <- list(GO_results$GMCSF$GMCSF_up$ID, GO_results$MCSF$MCSF_down$ID)
  x5 <- list(GO_results$GMCSF$GMCSF_down$ID, GO_results$MCSF$MCSF_up$ID)
  
  venn_plot <- function(list_data,
                        file_name,
                        plot_title) 
  {
    venn.diagram(x = list_data,
                 filename = file_name,
                 category.names = c("", ""),
                 main = plot_title,
                 output = TRUE ,
                 imagetype = "png" ,
                 height = 2000, 
                 width = 2000, 
                 resolution = 1000,
                 compression = "lzw",
                 lwd = 1,
                 col = c('coral', "steelblue3"),
                 fill = c(alpha("coral",0.3), alpha("steelblue3",0.3)),
                 disable.logging = TRUE,
                 cex = 0.5,
                 fontfamily = "sans",
                 cat.cex = 0.3,
                 cat.default.pos = "outer",
                 cat.pos = c(-27, 27),
                 cat.dist = c(0.055, 0.055),
                 cat.fontfamily = "sans")
  }
  
  # change to the plot directory
  setwd("plots/Venn_plots")
  
  # Make the plots
  venn_plot(list_data = x1, file_name = "all_GOs_venn.png", plot_title = "A")
  
  venn_plot(list_data = x2, file_name = "up_GOs_venn.png", plot_title = "B")
  
  venn_plot(list_data = x3, file_name = "down_GOs_venn.png", plot_title = "C")
  
  venn_plot(list_data = x4, file_name = "Gup_Mdown_GOs_venn.png", plot_title = "D")
  
  venn_plot(list_data = x5, file_name = "Gdown_Mup_GOs_venn.png", plot_title = "E")
  
  # Import PNG files
  img1 <- readPNG("all_GOs_venn.png")
  img2 <- readPNG("up_GOs_venn.png")
  img3 <- readPNG("down_GOs_venn.png")
  img4 <- readPNG("Gup_Mdown_GOs_venn.png")
  img5 <- readPNG("Gdown_Mup_GOs_venn.png")
  
  # Create a legend
  legend_labels <- c("GM-CSF", "M-CSF")
  legend_colors <- c("coral", "steelblue3")
  
  # Create a blank ggplot with the legend
  dummy_legend_plot <- ggplot() +
    geom_point(aes(x = 0, y = 0, color = factor(legend_labels))) +
    scale_color_manual(values = legend_colors, name = "Treatment") +  # Set the legend title
    theme_void() +
    theme(legend.position = "bottom",  # Change to "bottom" for a horizontal legend
          legend.direction = "vertical",  # Change to "horizontal" for a horizontal legend
          legend.box = "vertical")  # Place legend labels under one another
  
  # Extract the legend
  legend_grob <- cowplot::get_legend(dummy_legend_plot)
  
  # Create rasterGrob objects for each image
  grob_img1 <- rasterGrob(img1, interpolate = TRUE)
  grob_img2 <- rasterGrob(img2, interpolate = TRUE)
  grob_img3 <- rasterGrob(img3, interpolate = TRUE)
  grob_img4 <- rasterGrob(img4, interpolate = TRUE)
  grob_img5 <- rasterGrob(img5, interpolate = TRUE)
  
  # Arrange the plots and legend into a composite image
  composite_plot <- grid.arrange(
    grob_img1, 
    legend_grob,
    grob_img2, grob_img3, grob_img4, grob_img5, 
    ncol = 2,
    nrow = 3,
    widths = c(0.5, 0.5),  # Adjust the widths as needed
    heights = c(1, 1, 1)
  )
  
  # Save the composite image with legend
  png("GO_venns", width = 6, height = 6, units = "in", res = 600)
  cowplot::save_plot("GO_venns.png", composite_plot, base_height = 6, base_asp = 1)
  dev.off()  # Close the graphics device
  
  # return to the working directory
  setwd(mainDir)
}

##### Circular heatmap of specified GO term #####
# recommended path for mainDir: mainDir <- "~/Matthew_Masters/Work"
circular_GO_heatmap <- function(mainDir, # path to the 'Work' directory
                                GO_ID, 
                                plot_name
) {
  # Read in necessary packages
  deseq_box_packages <- c("ggplot2", "DESeq2", "dplyr", "tidyverse", "data.table", "ggpattern", "ComplexHeatmap", "circlize")
  suppressMessages(lapply(deseq_box_packages, require, character.only = TRUE))
  
  # Set the working directory
  setwd(mainDir)
  
  # Read in GO enrichment results
  GO_results <- readRDS("clusterProfiler/GO_enrichment/GO_summary.rds")
  
  # Read in DEGs
  GMCSF_DEGs <- read.csv("DEG_lists/treatment_GMCSF_vs_none_DEGs.csv", row.names = 1)
  MCSF_DEGs <- read.csv("DEG_lists/treatment_MCSF_vs_none_DEGs.csv", row.names = 1)
  
  GMCSF_GO <- GO_results[["GMCSF"]][["GMCSF_both"]]
  MCSF_GO <- GO_results[["MCSF"]][["MCSF_both"]]
  
  # Extract a specific column from a dataframe based on  a particular character in a specified column
  GMCSF_GO <- GMCSF_GO[GMCSF_GO$ID %like% GO_ID, 9]
  GMCSF_GO <- data.frame(str_split(GMCSF_GO, "/"))
  colnames(GMCSF_GO) <- "Gene_name"
  
  # Extract a specific column from a dataframe based on  a particular character in a specified column
  MCSF_GO <- MCSF_GO[MCSF_GO$ID %like% GO_ID, 9]
  MCSF_GO <- data.frame(str_split(MCSF_GO, "/"))
  colnames(MCSF_GO) <- "Gene_name"
  
  source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
  generate_ensembl(mainDir)
  
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "external_gene_name", # The ID type you want in your final object.
             object = GMCSF_DEGs$Gene_name) # character vector containing the 'initial ID' type that you want to convert
  GMCSF_DEGs[,1] <- object_final$external_gene_name
  
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "external_gene_name", # The ID type you want in your final object.
             object = MCSF_DEGs$Gene_name) # character vector containing the 'initial ID' type that you want to convert
  MCSF_DEGs[,1] <- object_final$external_gene_name
  
  GMCSF_common <- semi_join(GMCSF_GO, MCSF_GO)
  GMCSF_common <- left_join(GMCSF_common, GMCSF_DEGs, by = "Gene_name")
  GMCSF_common <- data.frame(rep("GMCSF", nrow(GMCSF_common)), GMCSF_common$Gene_name, GMCSF_common$log2FoldChange)
  colnames(GMCSF_common) <- c("Treatment", "Gene", "LFC")
  
  MCSF_common <- semi_join(MCSF_GO, GMCSF_GO)
  MCSF_common <- left_join(MCSF_common, MCSF_DEGs, by = "Gene_name")
  MCSF_common <- data.frame(rep("MCSF", nrow(MCSF_common)), MCSF_common$Gene_name, MCSF_common$log2FoldChange)
  colnames(MCSF_common) <- c("Treatment", "Gene", "LFC")
  
  # Unique GO genes
  GMCSF_unique <- setdiff(GMCSF_GO$Gene_name, MCSF_GO$Gene_name)
  MCSF_unique <- setdiff(MCSF_GO$Gene_name, GMCSF_GO$Gene_name)
  unique_genes <- c(GMCSF_unique, MCSF_unique)
  
  GMCSF_unique <- data.frame(Gene_name = unique_genes)
  GMCSF_unique <- left_join(GMCSF_unique, GMCSF_DEGs, by = "Gene_name")
  GMCSF_unique <- data.frame(rep("GMCSF", nrow(GMCSF_unique)), GMCSF_unique$Gene_name, GMCSF_unique$log2FoldChange)
  colnames(GMCSF_unique) <- c("Treatment", "Gene", "LFC")
  GMCSF_unique[is.na(GMCSF_unique)] <- 0
  
  
  MCSF_unique <- data.frame(Gene_name = unique_genes)
  MCSF_unique <- left_join(MCSF_unique, MCSF_DEGs, by = "Gene_name")
  MCSF_unique <- data.frame(rep("MCSF", nrow(MCSF_unique)), MCSF_unique$Gene_name, MCSF_unique$log2FoldChange)
  colnames(MCSF_unique) <- c("Treatment", "Gene", "LFC")
  MCSF_unique[is.na(MCSF_unique)] <- 0
  
  # Create numeric matrices for the commonly differentially expressed genes
  mat_GMCSF_common <- as.matrix(GMCSF_common[, c("LFC")])
  rownames(mat_GMCSF_common) <- GMCSF_common$Gene
  
  mat_MCSF_common <- as.matrix(MCSF_common[, c("LFC")])
  rownames(mat_MCSF_common) <- MCSF_common$Gene
  
  # Create numeric matrices for the uniquely differentially expressed genes
  mat_GMCSF_unique <- as.matrix(GMCSF_unique[, c("LFC")])
  rownames(mat_GMCSF_unique) <- GMCSF_unique$Gene
  
  mat_MCSF_unique <- as.matrix(MCSF_unique[, c("LFC")])
  rownames(mat_MCSF_unique) <- MCSF_unique$Gene
  
  # Combine matrices
  heat_matrix <- rbind(cbind(mat_MCSF_common, mat_GMCSF_common), cbind(mat_MCSF_unique, mat_GMCSF_unique))
  
  colnames(heat_matrix) <- c("MCSF", "GMCSF")
  
  # Create the split factor for the heatmap
  split <- c(rep("Common", length(mat_GMCSF_common)), rep("Unique", length(mat_GMCSF_unique)))
  split <- factor(split, levels = c("Common", "Unique"))
  
  circos.clear()
  png(sprintf("plots/%s.png", plot_name), width = 5000, height = 5000, res = 600)
  col_fun1 = colorRamp2(c(-5, 0, 5), c("steelblue3", "grey95", "coral"))
  circos.par(start.degree = 90, gap.degree = 10)
  circos.heatmap(heat_matrix, 
                 split = split, 
                 col = col_fun1, 
                 track.height = 0.25, 
                 bg.lwd = 2, 
                 bg.lty = 2, 
                 show.sector.labels = TRUE, 
                 rownames.side = "inside",
                 rownames.cex = 0.55,
                 na.col = "grey75")
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if(CELL_META$sector.numeric.index == 2) { # the last sector
      cn = c("GM-CSF", "M-CSF")
      n = length(cn)
      circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(2, "mm"), 
                  1:n - 0.5, cn, 
                  cex = 0.5, adj = c(0, 0.5), facing = "inside")
    }
  }, bg.border = NA)
  lgd = Legend(title = "Log2 Fold Change", col_fun = col_fun1, title_position = c("topcenter"))
  grid.draw(lgd)
  dev.off()
}
##### Generate Supplementary tables #####
# recommended path for mainDir: mainDir <- "~/Matthew_Masters/Work"
circular_GO_heatmap <- function(mainDir, # path to the 'Work' directory
                                ) 
  {
  # Read in necessary packages
  deseq_packages <- c("DESeq2", "dplyr", "tidyverse", "data.table", "R.utils")
  suppressMessages(lapply(deseq_packages, require, character.only = TRUE))
  
  # Set the working directory
  setwd(mainDir)
  
  ##### DEGs (Tables S1 and S2) #####
  DEGs <- readRDS("DEG_lists/GMCSF_and_MCSF_comparison.rds")
  
  # Convert ENSG ID to gene name
  source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
  generate_ensembl(mainDir)
  
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "external_gene_name", # The ID type you want in your final object.
             object = DEGs$GMCSF$Gene_name # character vector containing the 'initial ID' type that you want to convert
  )
  DEGs$GMCSF$Gene_name <- object_final$external_gene_name
  write.csv(DEGs[["GMCSF"]], file = "~/Matthew_Masters/Docs/dissertation_tables/S2_GMCSF_DEGs.csv")
  
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "external_gene_name", # The ID type you want in your final object.
             object = DEGs$MCSF$Gene_name # character vector containing the 'initial ID' type that you want to convert
  )
  DEGs$MCSF$Gene_name <- object_final$external_gene_name
  write.csv(DEGs[["MCSF"]], file = "~/Matthew_Masters/Docs/dissertation_tables/S1_MCSF_DEGs.csv")
  
  ##### Full enrichment results (Tables S3 and S4) #####
  
  GO_results <- readRDS("clusterProfiler/GO_enrichment/GO_summary.rds")
  df1 <- GO_results[["GMCSF"]][["GMCSF_both"]]
  df1 <- df1[, c(2:10)]
  write.csv(df1, file = "~/Matthew_Masters/Docs/dissertation_tables/S4_GMCSF_GO.csv", row.names = FALSE)
  
  df2 <- GO_results[["MCSF"]][["MCSF_both"]]
  df2 <- df2[, c(2:10)]
  write.csv(df2, file = "~/Matthew_Masters/Docs/dissertation_tables/S3_MCSF_GO.csv", row.names = FALSE)
  
  ##### DETFs (Tables S5 and S6) #####
  DETFs <- readRDS("DEG_lists/GMCSF_tfs_and_MCSF_tfs_comparison.rds")
  
  # Convert ENSG ID to gene name
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "external_gene_name", # The ID type you want in your final object.
             object = DETFs$GMCSF_tfs$Gene_name # character vector containing the 'initial ID' type that you want to convert
  )
  DETFs$GMCSF_tfs$Gene_name <- object_final$external_gene_name
  write.csv(DETFs[["GMCSF_tfs"]], file = "~/Matthew_Masters/Docs/dissertation_tables/S6_GMCSF_DETFs.csv")
  
  convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
             final_ID = "external_gene_name", # The ID type you want in your final object.
             object = DETFs$MCSF_tfs$Gene_name # character vector containing the 'initial ID' type that you want to convert
  )
  DETFs$MCSF_tfs$Gene_name <- object_final$external_gene_name
  write.csv(DETFs[["MCSF_tfs"]], file = "~/Matthew_Masters/Docs/dissertation_tables/S5_MCSF_DETFs.csv")
  
  ##### TFs associated with myeloid cell differentiation or the regulation of inflammatory response GO terms (Table S7) #####
  
  # Regulation of inflammatory response
  df1 <- GO_results[["GMCSF"]][["GMCSF_both"]]
  df2 <- GO_results[["MCSF"]][["MCSF_both"]]
  
  # Extract a specific column from a dataframe based on  a particular character in a specified column
  df1 <- df1[df1$ID %like% "GO:0050727", 9]
  df1 <- data.frame(str_split(df1, "/"))
  colnames(df1) <- "Gene_name"
  
  # Extract a specific column from a dataframe based on  a particular character in a specified column
  df2 <- df2[df2$ID %like% "GO:0050727", 9]
  df2 <- data.frame(str_split(df2, "/"))
  colnames(df2) <- "Gene_name"
  
  GMCSF_TFs <- na.omit(left_join(df1, DETFs[["GMCSF_tfs"]], by = "Gene_name"))
  MCSF_TFs <- na.omit(left_join(df2, DETFs[["MCSF_tfs"]], by = "Gene_name"))
  Genes <- unique(c(GMCSF_TFs$Gene_name, MCSF_TFs$Gene_name))
  inflamm_table <- data.frame(Gene_name = Genes, GO_term = c(rep("Regulation of inflammatory response", length(Genes))))
  
  
  # Myeloid cell differentiation
  df1 <- GO_results[["GMCSF"]][["GMCSF_both"]]
  df2 <- GO_results[["MCSF"]][["MCSF_both"]]
  
  # Extract a specific column from a dataframe based on  a particular character in a specified column
  df1 <- df1[df1$ID %like% "GO:0030099", 9]
  df1 <- data.frame(str_split(df1, "/"))
  colnames(df1) <- "Gene_name"
  
  # Extract a specific column from a dataframe based on  a particular character in a specified column
  df2 <- df2[df2$ID %like% "GO:0030099", 9]
  df2 <- data.frame(str_split(df2, "/"))
  colnames(df2) <- "Gene_name"
  
  GMCSF_TFs <- na.omit(left_join(df1, DETFs[["GMCSF_tfs"]], by = "Gene_name"))
  MCSF_TFs <- na.omit(left_join(df2, DETFs[["MCSF_tfs"]], by = "Gene_name"))
  Genes <- unique(c(GMCSF_TFs$Gene_name, MCSF_TFs$Gene_name))
  myeloid_table <- data.frame(Gene_name = Genes, GO_term = c(rep("Myeloid cell differentiation", length(Genes))))
  
  # Combine the two tables
  combined <- rbind(myeloid_table, inflamm_table)
  
  write.csv(combined, file = "~/Matthew_Masters/Docs/dissertation_tables/S7_GO_DETFs.csv", row.names = FALSE)
}

