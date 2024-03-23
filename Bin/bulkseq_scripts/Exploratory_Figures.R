############ Rough Script for exploratory figures ##################
############ DESeq2 results ##############
##### examine lists of DEGs #####
mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)

DEGs <- readRDS("DEG_lists/GMCSF_and_MCSF_comparison.rds")

# Convert ENSG ID to gene name
source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
generate_ensembl(mainDir)

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
         final_ID = "external_gene_name", # The ID type you want in your final object.
         object = DEGs$GMCSF$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DEGs$GMCSF$Gene_name <- object_final$external_gene_name
View(DEGs[["GMCSF"]])

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
         final_ID = "external_gene_name", # The ID type you want in your final object.
         object = DEGs$MCSF$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DEGs$MCSF$Gene_name <- object_final$external_gene_name
View(DEGs[["MCSF"]])

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
         final_ID = "external_gene_name", # The ID type you want in your final object.
         object = DEGs$GMCSF_MCSF_common$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DEGs$GMCSF_MCSF_common$Gene_name <- object_final$external_gene_name
View(DEGs[["GMCSF_MCSF_common"]])

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
         final_ID = "external_gene_name", # The ID type you want in your final object.
         object = DEGs$GMCSF_unique$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DEGs$GMCSF_unique$Gene_name <- object_final$external_gene_name
View(DEGs[["GMCSF_unique"]])

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
         final_ID = "external_gene_name", # The ID type you want in your final object.
         object = DEGs$MCSF_unique$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DEGs$MCSF_unique$Gene_name <- object_final$external_gene_name
View(DEGs[["MCSF_unique"]])

##### examine the normalised count matrix #####
mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)

DESeq_package <- c("DESeq2", "dplyr", "ggplot2", "RColorBrewer", "AnnotationDbi", "org.Hs.eg.db")
suppressMessages(lapply(DESeq_package, require, character.only = TRUE))

dds <- readRDS("DESeq2/dds.rds")
cts <- counts(dds, normalized=TRUE)

# Convert ENSG ID to gene name
source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
generate_ensembl(mainDir)

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
         final_ID = "external_gene_name", # The ID type you want in your final object.
         object = rownames(cts) # character vector containing the 'initial ID' type that you want to convert
)
match1 <- object_final

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
         final_ID = "description", # The ID type you want in your final object.
         object = rownames(cts) # character vector containing the 'initial ID' type that you want to convert
)
match2 <- object_final

cts <- data.frame(match1$ensembl_gene_id, match1$external_gene_name, match2$description, cts[ ,1:9])


##### plot counts for specified set of genes #####
# I need to add significance bars to each plot. See this vignette: https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/
mainDir <- "~/Matthew_Masters/Work"

# Read in necessary packages
deseq_box_packages <- c("ggplot2", "DESeq2", "dplyr", "tidyverse", "data.table", "ggpubr")
suppressMessages(lapply(deseq_box_packages, require, character.only = TRUE))

# Set the working directory
setwd(mainDir)

# Read in lists of DEGs
DEGs <- readRDS("~/Work_old/L2FC_1_results/DEG_lists/GMCSF_and_MCSF_comparison.rds")

# Convert ENSG ID to gene name
source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
generate_ensembl(mainDir)

# Manually provde the set of genes (M0 markers) to visualize.
gene_names <- c("ENSG00000174837", # ADGRE1
                "ENSG00000121807", # CCR2
                "ENSG00000170458", # CD14
                "ENSG00000129226", # CD68
                "ENSG00000177575", # CD163
                "ENSG00000182578", # CSF1R
                "ENSG00000168329", # CX3CR1
                "ENSG00000203747", # FCGR3A
                "ENSG00000169896", # ITGAM
                "ENSG00000260314", # MRC1
                "ENSG00000088827", # SIGLEC1
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
d5 <- gen_plot_data(5)
d6 <- gen_plot_data(6)
d7 <- gen_plot_data(7)
d8 <- gen_plot_data(8)
d9 <- gen_plot_data(9)
d10 <- gen_plot_data(10)
d11 <- gen_plot_data(11)
d12 <- gen_plot_data(12)

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

p1 <- plot_DEGs(d1, gene_number = 1) + ylab(NULL) + xlab(NULL)
p2 <- plot_DEGs(d2, gene_number = 2) + ylab(NULL) + xlab(NULL)
p3 <- plot_DEGs(d3, gene_number = 3) + ylab(NULL)
p4 <- plot_DEGs(d4, gene_number = 4) + ylab(NULL) + xlab(NULL)
p5 <- plot_DEGs(d5, gene_number = 5) + ylab(NULL)
p6 <- plot_DEGs(d6, gene_number = 6) + xlab(NULL)
p7 <- plot_DEGs(d7, gene_number = 7) + ylab(NULL) + xlab(NULL)
p8 <- plot_DEGs(d8, gene_number = 8) + ylab(NULL)
p9 <- plot_DEGs(d9, gene_number = 9) + ylab(NULL) + xlab(NULL)
p10 <- plot_DEGs(d10, gene_number = 10) + ylab(NULL) + xlab(NULL)
p11 <- plot_DEGs(d11, gene_number = 11)
p12 <- plot_DEGs(d12, gene_number = 12) + xlab(NULL)

cowplot::plot_grid(p12, p2, p4, p9, p6, p10, p1, p7, p11, p3, p5, p8, # Plots are ordered to better show trends
                   nrow = 3, ncol = 4)

# Save the plot to file
ggsave(sprintf("%s/plots/surface_markers.png", mainDir), 
       dpi = 300, 
       width = 12, 
       height = 8, 
       units = "in")

##### plot log2 fold changes for specified set of genes #####
mainDir <- "~/Matthew_Masters/Work"

# Read in necessary packages
deseq_box_packages <- c("ggplot2", "DESeq2", "dplyr", "tidyverse", "data.table", "ggpattern")
suppressMessages(lapply(deseq_box_packages, require, character.only = TRUE))

# Set the working directory
setwd(mainDir)

# Read in lists of DEGs
DEGs <- readRDS("~/Work_old/L2FC_1_results/DEG_lists/GMCSF_and_MCSF_comparison.rds")

# Convert ENSG ID to gene name
source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
generate_ensembl(mainDir)

# Alternative set of genes (updated).
gene_names <- c("ENSG00000169896", # ITGAM
                "ENSG00000129226", # CD68
                "ENSG00000203747", # FCG3RA
                "ENSG00000170458", # CD14
                "ENSG00000182578", # CSFR1
                "ENSG00000177575", # CD163
                "ENSG00000168329", # CX3CR1
                "ENSG00000260314", # MRC1
                "ENSG00000121807", # CCR2
                "ENSG00000088827"  # SIGLEC1
                )

# Second alternative set
gene_names <- c("ENSG00000169896", # ITGAM
                "ENSG00000129226", # CD68
                "ENSG00000174837", # ADGRE1
                "ENSG00000170458", # CD14
                "ENSG00000153208", # MERTK
                "ENSG00000168329", # CX3CR1
                "ENSG00000121807", # CCR2
                "ENSG00000163606"  # CD200R1
)


convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = gene_names # character vector containing the 'initial ID' type that you want to convert
)
gene_names <- data.frame(object_final)
colnames(gene_names) <- c("Gene_name", "Gene_symbol")

GMCSF_DEGs <- left_join(gene_names, DEGs$GMCSF, by = "Gene_name")

MCSF_DEGs <- left_join(gene_names, DEGs$MCSF, by = "Gene_name")

# Read in dds object with normalised counts for calculating log2 fold changes manually
dds <- readRDS("DESeq2/dds.rds")

# Manually calculate and replace log2 fold changes.
for (i in c(1:length(gene_names$Gene_name))) {
  print(gene_names[i, 1])
  gene_counts <- plotCounts(dds, 
                            gene = as.character(gene_names[i, 1]), 
                            intgroup= "treatment", returnData = TRUE)
  
  GMCSF_DEGs[i,3] <- log2((gene_counts[7,1] + gene_counts[8,1] + gene_counts[9,1])/(gene_counts[1,1] + gene_counts[2,1] + gene_counts[3,1]))*
    abs((((gene_counts[7,1] + gene_counts[8,1] + gene_counts[9,1])/3)-((gene_counts[1,1] + gene_counts[2,1] + gene_counts[3,1])/3))/1000)
  MCSF_DEGs[i,3] <- log2((gene_counts[4,1] + gene_counts[5,1] + gene_counts[6,1])/(gene_counts[1,1] + gene_counts[2,1] + gene_counts[3,1]))*
    abs((((gene_counts[4,1] + gene_counts[5,1] + gene_counts[6,1])/3)-((gene_counts[1,1] + gene_counts[2,1] + gene_counts[3,1])/3))/1000)
}

data <- bind_rows(GMCSF_DEGs, MCSF_DEGs, .id = 'treatment')
data$treatment <- c(rep("GMCSF", length(gene_names$Gene_name)), rep("MCSF", length(gene_names$Gene_name)))

ggplot(data, aes(x = Gene_symbol, y = log2FoldChange, fill = padj)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  facet_wrap(~ treatment, scales = "free_y", ncol = 2) +
  labs(x = "Gene symbol",
    y = "log2 Fold Change",
    fill = "Adjusted P-value") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "bottom"  # Move legend to the bottom
  ) +
  scale_fill_gradient(low = "#371ea3", high = "#46bac2") +  # Adjust fill scale colors
  facet_wrap(~ treatment, scales = "fixed", ncol = 2)  # Set fixed scales for both facets

# Save the plot to file
ggsave(sprintf("%s/plots/surface_markers.png", mainDir), 
       dpi = 300, 
       width = 10, 
       height = 8, 
       units = "in")
##### Linear heatmaps for specified set of genes #####
mainDir <- "~/Matthew_Masters/Work"

# Read in necessary packages
deseq_box_packages <- c("ggplot2", "DESeq2", "dplyr", "tidyverse", "data.table", "ggpattern", "ComplexHeatmap", "circlize")
suppressMessages(lapply(deseq_box_packages, require, character.only = TRUE))

# Set the working directory
setwd(mainDir)

# Read in lists of DEGs
DEGs <- readRDS("~/Work_old/L2FC_1_results/DEG_lists/GMCSF_and_MCSF_comparison.rds")

# Convert ENSG ID to gene name
source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
generate_ensembl(mainDir)

# M1 markers
M1_markers <- c("ENSG00000007171", # NOS2
                "ENSG00000137462", # TLR2 = CD282
                "ENSG00000121594", # CD80
                "ENSG00000164400", # CSF2 = GM-CSF
                "ENSG00000232810", # TNF = TNF alpha
                "ENSG00000125538", # IL1B
                "ENSG00000136244", # IL6
                "ENSG00000136869", # TLR4 = CD284
                "ENSG00000114013", # CD86 = B7-2
                "ENSG00000081041", # CXCL2 = MIP-2a
                "ENSG00000111537", # IFNG = IFN gamma
                "ENSG00000115594"  # IL1R1
)

# M2 markers
M2_markers <- c("ENSG00000177575", # CD163
                "ENSG00000182578", # CSF1R = CD115
                "ENSG00000260314", # MRC1 = CD206/mannose receptor
                "ENSG00000132170", # PPARG
                "ENSG00000118520", # ARG1
                "ENSG00000132514", # CLEC10A = mgl2/CD301
                "ENSG00000172243", # CLEC7A = Dectin-1/CD369
                "ENSG00000197646", # PDCD1LG2 = pdl2
                "ENSG00000102962", # CCL22
                "ENSG00000136634", # IL10
                "ENSG00000203747", # FCGR3A = CD16
                "ENSG00000113520", # IL4
                "ENSG00000137265", # IRF4 = MUM1
                "ENSG00000100311", # PDGFB
                "ENSG00000166888" # STAT6
)

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = M1_markers # character vector containing the 'initial ID' type that you want to convert
)
M1_markers <- data.frame(object_final)
colnames(M1_markers) <- c("Gene_name", "Gene_symbol")

GMCSF_M1_DEGs <- left_join(M1_markers, DEGs$GMCSF, by = "Gene_name")
MCSF_M1_DEGs <- left_join(M1_markers, DEGs$MCSF, by = "Gene_name")

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = M2_markers # character vector containing the 'initial ID' type that you want to convert
)

M2_markers <- data.frame(object_final)
colnames(M2_markers) <- c("Gene_name", "Gene_symbol")

GMCSF_M2_DEGs <- left_join(M2_markers, DEGs$GMCSF, by = "Gene_name")
MCSF_M2_DEGs <- left_join(M2_markers, DEGs$MCSF, by = "Gene_name")

# Read in dds object with normalised counts for calculating log2 fold changes manually
dds <- readRDS("DESeq2/dds.rds")

# Create empty vectors for log2 fold changes
M1_DEGs <- data.frame(Gene_name = M1_markers$Gene_name, Gene_symbol = M1_markers$Gene_symbol, GMCSF_log2FC = NA, MCSF_log2FC = NA)

# Manually calculate and replace log2 fold changes for the M1 markers.
for (i in seq_len(nrow(M1_markers))) {
  print(M1_markers[i, 1])
  
  gene_id <- as.character(M1_markers[i, 1])
  
  # Check if the gene is present in the count matrix
  if (gene_id %in% rownames(assay(dds))) {
    gene_counts <- plotCounts(dds, gene = gene_id, intgroup = "treatment", returnData = TRUE)
    
    GMCSF_M1_DEGs[i, 3] <- log2((gene_counts[7, 1] + gene_counts[8, 1] + gene_counts[9, 1]) / (gene_counts[1, 1] + gene_counts[2, 1] + gene_counts[3, 1]))
    MCSF_M1_DEGs[i, 3] <- log2((gene_counts[4, 1] + gene_counts[5, 1] + gene_counts[6, 1]) / (gene_counts[1, 1] + gene_counts[2, 1] + gene_counts[3, 1]))
  } else {
    # Gene not present, skip and assign NA to log2 fold changes
    GMCSF_M1_DEGs[i, 3] <- NA
    MCSF_M1_DEGs[i, 3] <- NA
    print("Gene not found in count matrix, skipping.")
  }
}

M2_DEGs <- data.frame(Gene_name = M2_markers$Gene_name, Gene_symbol = M2_markers$Gene_symbol, GMCSF_log2FC = NA, MCSF_log2FC = NA)

# Manually calculate and replace log2 fold changes for the M2 markers.
for (i in seq_len(nrow(M2_markers))) {
  print(M2_markers[i, 1])
  
  gene_id <- as.character(M2_markers[i, 1])
  
  # Check if the gene is present in the count matrix
  if (gene_id %in% rownames(assay(dds))) {
    gene_counts <- plotCounts(dds, gene = gene_id, intgroup = "treatment", returnData = TRUE)
    
    GMCSF_M2_DEGs[i, 3] <- log2((gene_counts[7, 1] + gene_counts[8, 1] + gene_counts[9, 1]) / (gene_counts[1, 1] + gene_counts[2, 1] + gene_counts[3, 1]))
    MCSF_M2_DEGs[i, 3] <- log2((gene_counts[4, 1] + gene_counts[5, 1] + gene_counts[6, 1]) / (gene_counts[1, 1] + gene_counts[2, 1] + gene_counts[3, 1]))
  } else {
    # Gene not present, skip and assign NA to log2 fold changes
    GMCSF_M2_DEGs[i, 3] <- NA
    MCSF_M2_DEGs[i, 3] <- NA
    print("Gene not found in count matrix, skipping.")
  }
}

M1_data <- bind_rows(GMCSF_M1_DEGs, MCSF_M1_DEGs, .id = 'treatment')
M1_data$treatment <- c(rep("GMCSF", nrow(M1_markers)), rep("MCSF", nrow(M1_markers)))

M2_data <- bind_rows(GMCSF_M2_DEGs, MCSF_M2_DEGs, .id = 'treatment')
M2_data$treatment <- c(rep("GMCSF", nrow(M2_markers)), rep("MCSF", nrow(M2_markers)))

tile_size <- 10

heat_M1 <- ggplot(M1_data, aes(x = Gene_symbol, y = treatment, fill = log2FoldChange)) +
  geom_tile() +
  scale_fill_gradient2(low = "steelblue2", mid = "white", high = "coral", na.value = "grey80")

heat_M2 <- ggplot(M2_data, aes(x = Gene_symbol, y = treatment, fill = log2FoldChange)) +
  geom_tile() +
  scale_fill_gradient2(low = "steelblue2", mid = "white", high = "coral", na.value = "grey80")

cowplot::plot_grid(heat_M1, heat_M2, nrow=2)

# Save the plot to file
ggsave(sprintf("%s/plots/macrophage_subtype_markers.png", mainDir), 
       dpi = 300, 
       width = 12, 
       height = 6, 
       units = "in")


##### Examine transcription factor results #####
mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)

DETFs <- readRDS("DEG_lists/GMCSF_tfs_and_MCSF_tfs_comparison.rds")

# Convert ENSG ID to gene name
source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
generate_ensembl(mainDir)

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
         final_ID = "external_gene_name", # The ID type you want in your final object.
         object = DETFs$GMCSF_tfs$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DETFs$GMCSF_tfs$Gene_name <- object_final$external_gene_name
View(DETFs[["GMCSF_tfs"]])

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
         final_ID = "external_gene_name", # The ID type you want in your final object.
         object = DETFs$MCSF_tfs$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DETFs$MCSF_tfs$Gene_name <- object_final$external_gene_name
View(DETFs[["MCSF_tfs"]])

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
         final_ID = "external_gene_name", # The ID type you want in your final object.
         object = DETFs$GMCSF_tfs_MCSF_tfs_common$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DETFs$GMCSF_tfs_MCSF_tfs_common$Gene_name <- object_final$external_gene_name
View(DETFs[["GMCSF_tfs_MCSF_tfs_common"]])

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
         final_ID = "external_gene_name", # The ID type you want in your final object.
         object = DETFs$GMCSF_tfs_unique$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DETFs$GMCSF_tfs_unique$Gene_name <- object_final$external_gene_name
View(DETFs[["GMCSF_tfs_unique"]])

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
         final_ID = "external_gene_name", # The ID type you want in your final object.
         object = DETFs$MCSF_tfs_unique$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DETFs$MCSF_tfs_unique$Gene_name <- object_final$external_gene_name
View(DETFs[["MCSF_tfs_unique"]])

### Identify genes that also play a role in myeloid cell differentiation
# Move to the working directory
setwd(mainDir)
DEG_list <- sprintf("%s/DEG_lists/treatment_MCSF_vs_none_DEGs.csv", mainDir)

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
differentiation_factors <- object_final[grepl("myeloid cell differentiation", object_final$name_1006), ]

# Further filter this list to retain unique rows and essential information
differentiation_factors <- data.frame(differentiation_factors[1])
differentiation_factors <- unique(differentiation_factors)
colnames(differentiation_factors)[1] <- "Gene_name"

# Regain information on log2 fold change.
final_df <- left_join(differentiation_factors, df, by = "Gene_name")

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = final_df$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
final_df[,1] <- object_final[,2]
GMCSF_diff_genes <- final_df
MCSF_diff_genes <- final_df

# Write the list of DE transcription factors to file
write.csv(final_df, file = sprintf("%s/DEG_lists/%s_tfs.csv", mainDir, output_name))
print("Transcription factors have been identified for the supplied set of differentially expressed genes.")
print(sprintf("The resultant file can be found at '%s/DEG_lists/%s_tfs.csv'", mainDir, output_name))
print("------------------------------")

############ ClusterProfiler results ##############
mainDir <- "~/Matthew_Masters/Work/clusterProfiler_bulk"
setwd(mainDir)

# How can I test for common GO terms between 2 given csvs? 
# lets look at downregulated GMCSF GOs which are upregulated in MCSF
GMCSF_GO <- read.csv("treatment_GMCSF_vs_none_enrichment/GO_DOWN_BP_significant.csv")
MCSF_GO <- read.csv("treatment_MCSF_vs_none_enrichment/GO_UP_BP_significant.csv")

GMCSF_DOWN_MCSF_UP <- intersect(GMCSF_GO$Description, MCSF_GO$Description)

# Now for upregulated in GMCSF and Downregualted in MCSF
GMCSF_GO <- read.csv("treatment_GMCSF_vs_none_enrichment/GO_UP_BP_significant.csv")
MCSF_GO <- read.csv("treatment_MCSF_vs_none_enrichment/GO_DOWN_BP_significant.csv")

GMCSF_UP_MCSF_DOWN <- intersect(GMCSF_GO$Description, MCSF_GO$Description)


# Looking at the whole-DEG lists, which GO terms are shared and which are unique?
GMCSF_both <- readRDS("treatment_GMCSF_vs_none_enrichment/GO_BOTH_BP.rds")
GMCSF_both <- GMCSF_both@result
GMCSF_both <- subset(GMCSF_both, `p.adjust` < 0.05)
GMCSF_GOs <- rownames(GMCSF_both)

MCSF_both <- readRDS("treatment_MCSF_vs_none_enrichment/GO_BOTH_BP.rds")
MCSF_both <- MCSF_both@result
MCSF_both <- subset(MCSF_both, `p.adjust` < 0.05)
MCSF_GOs <- rownames(MCSF_both)

common_GOs <- intersect(GMCSF_GOs, MCSF_GOs)

GMCSF_unique_GOs <- setdiff.Vector(GMCSF_GOs, common_GOs)

MCSF_unique_GOs <- setdiff.Vector(MCSF_GOs, common_GOs)

### Focusing on inflammation GO (GO:0050727)
MCSF_both$ID == "GO:0050727"
GMCSF_inflam <- GMCSF_both[71, 8]
GMCSF_inflam <- 
  
MCSF_inflam <- MCSF_both[86, 8]
common_inflam <- intersect(GMCSF_inflam, MCSF_inflam)

# Pseudocode for pulling out gene names based on a given GO term/ID
# gene names <- MCSF_both[GO ID, c(8)] GO ID is present as the row name and col 1 of the df. col 8 has the cell with all gene names.
# gene names is number a vector with 1 character. the gene names are separated by '/'. so now i must split these using that character.

############ Bulkseq Figures plan ##############


###### Tximport ######
mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)
cts <- readRDS("tximport/txi.count.matrix.rds")

# Read in dds object with normalised counts
dds <- readRDS("misc_outputs/dds.rds")

##### Filtering
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
colSums(cts_mtx)
# subset genes where test was TRUE.
cts_abu <- cts_abu[keep,]
cts_mtx <- cts_mtx[keep,]
cts_len <- cts_len[keep,]

# Apply the changes to the cts object.
cts$abundance <- cts_abu
cts$counts <- cts_mtx
cts$length <- cts_len
txi <- as.data.frame(cts$counts)
summary(txi)
colSums(txi)

count <- c()
sample <- c()
data <- data.frame(count, sample)
colnames(data) <- c("count", "sample")
for (i in 1:ncol(txi)) {
  print(i)
  test <- as.data.frame(txi[,i])
  test[,2] <- rep(colnames(txi)[i], nrow(test))
  append(data, test)
}



p <- ggplot(test, aes(x=sample, y=count)) + 
  geom_boxplot() + 
  ylim(0, 75)

p
# Boxplot tests
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
head(ToothGrowth)

# Basic box plot
p <- ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_boxplot()
p
# Rotate the box plot
p + coord_flip()
# Notched box plot
ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_boxplot(notch=TRUE)
# Change outlier, color, shape and size
ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)

##### Tximport TPM plot #####
# Plotting the TPM in replicates (density plot + box plot)
cts <- readRDS("~/Matthew_Masters/Work/tximport/txi.count.matrix.rds")
x = cts[["abundance"]]
x = data.frame(x)
colnames(x) = colnames(cts$counts)

# Converting to long dataframe and remove NULL values
x$gene = row.names(x)
x = x %>% pivot_longer(!gene, names_to = "replicate", values_to = "abundance")
x = x %>% filter(abundance != 0)

# Preparing for facet
x$rep = (do.call('rbind', strsplit(as.character(x$replicate),'_',fixed=TRUE)))[,1]
x$condition = (do.call('rbind', strsplit(as.character(x$replicate),'_',fixed=TRUE)))[,2]

ggplot(x,aes(x = log10(abundance),y = -0.02)) +
  # Horizontal boxplot
  ggstance::geom_boxploth(aes(fill = condition),width = 0.03) +
  # density plot
  geom_density(aes(x = log10(abundance)),inherit.aes = F) +
  facet_grid(cols = vars(rep), rows = vars(condition)) +
  scale_fill_discrete() +
  geom_hline(yintercept = 0)+
  labs(y = "Proportion of genes",x = "log10(TPM)")

###### DESeq2 Volcano Plot######
mainDir <- "~/Matthew_Masters/Work/DESeq2"
setwd(mainDir)

GMCSF_DEGs <- readRDS("treatment_GMCSF_vs_none.rds")
GMCSF_DEGs <- as.data.frame(GMCSF_DEGs@listData)
Gene_name <- rownames(GMCSF_DEGs)
log2FoldChange <- GMCSF_DEGs$log2FoldChange
log10padj <- -log10(GMCSF_DEGs$padj)

df <- data.frame(Gene_name, log2FoldChange, log10padj)
colnames(df) <- c("Gene_name", "log2FoldChange", "log10padj")
df <- filter(df, log10padj < 50)

# Construct the plot
GMCSF_volcano_plot <- ggplot(data = df, aes(x = log2FoldChange, y = log10padj)) +
  geom_point() + 
  geom_hline(yintercept=-log10(0.05), col="red") + 
  ggtitle("GMCSF volcano plot")

MCSF_DEGs <- readRDS("treatment_MCSF_vs_none.rds")
MCSF_DEGs <- as.data.frame(MCSF_DEGs@listData)
Gene_name <- rownames(MCSF_DEGs)
log2FoldChange <- MCSF_DEGs$log2FoldChange
log10padj <- -log10(MCSF_DEGs$padj)

df <- data.frame(Gene_name, log2FoldChange, log10padj)
colnames(df) <- c("Gene_name", "log2FoldChange", "log10padj")
df <- filter(df, log10padj < 50)

# Construct the plot
MCSF_volcano_plot <- ggplot(data = df, aes(x = log2FoldChange, y = log10padj)) +
  geom_point() + 
  geom_hline(yintercept=-log10(0.05), col="red") + 
  ggtitle("MCSF volcano plot")
MCSF_volcano_plot + GMCSF_volcano_plot

##### DESeq2 Gene count plots #####
library("ggplot2")
mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)

dds <- readRDS("misc_outputs/dds.rds")

genes <- c("ENSG00000121807")

source("~/Matthew_Masters/Bin/bulkseq_scripts_V2/biomaRt_annotate_V2.R")
generate_ensembl(mainDir)

annotate(initial_ID = "external_gene_name", # This is the starting ID that you wish to change to another type.
         final_ID = "ensembl_gene_id", # The ID type you want in your final object.
         object = genes # character vector containing the 'initial ID' type that you want to convert
)
genes <- object_final[,2]

d <- plotCounts(dds, gene = genes[1], intgroup="treatment", returnData=TRUE)


ggplot(d, aes(x=treatment, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0))

p <- ggplot(d, aes(x = treatment, y = count, color = treatment)) + 
  geom_boxplot()
p


##### Exporting tables ######
library("finalfit")

dependent = "differ.factor"

# Specify explanatory variables of interest
explanatory = c("age", "sex.factor", 
                "extent.factor", "obstruct.factor", 
                "nodes")
colon_s %>% 
  select(age, sex.factor, extent.factor, obstruct.factor, nodes) %>% 
  names() -> explanatory

colon_s %>% 
  ff_glimpse(dependent, explanatory)

colon_s %>% 
  summary_factorlist(dependent, explanatory, 
                     p=TRUE, na_include=TRUE)
explanatory = c("age", "sex.factor", 
                "extent.factor", "nodes")

colon_s %>% 
  mutate(
    nodes = ff_label(nodes, "Lymph nodes involved")
  ) %>% 
  summary_factorlist(dependent, explanatory, 
                     p=TRUE, na_include=TRUE, 
                     add_dependent_label=TRUE) -> table1
table1
save(table1, file = "out.rda")

# Misc explorations

mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)

DESeq_package <- c("DESeq2", "dplyr", "ggplot2", "RColorBrewer", "AnnotationDbi", "org.Hs.eg.db")
suppressMessages(lapply(DESeq_package, require, character.only = TRUE))

dds <- readRDS("DESeq2/dds.rds")
cts <- counts(dds, normalized=TRUE)

keep <- rowSums(cts) > 1000
# summary of test outcome: number of genes in each class:
table(keep, useNA="always")

# subset genes where test was TRUE.
cts <- cts[keep,]
