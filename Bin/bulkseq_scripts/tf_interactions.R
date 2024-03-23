##### Identify protein-protein interations of biologically-relevant transcription factors. #####
mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)
KEGG_packages <- c("clusterProfiler", 
                   "AnnotationDbi", 
                   "org.Hs.eg.db",
                   "R.utils", 
                   "biomaRt",
                   "dplyr",
                   "ggplot2",
                   "data.table",
                   "tidyverse")

suppressMessages(lapply(KEGG_packages, require, character.only = TRUE))

# Read in and annotate the two lists of DEGs
DEGs <- readRDS("DEG_lists/GMCSF_and_MCSF_comparison.rds")

# Convert ENSG ID to gene name
source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
generate_ensembl(mainDir)

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = DEGs$GMCSF$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DEGs$GMCSF$Gene_name <- object_final$external_gene_name
GMCSF_DEGs <- DEGs$GMCSF

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = DEGs$MCSF$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DEGs$MCSF$Gene_name <- object_final$external_gene_name
MCSF_DEGs <- DEGs$MCSF

TFlist <- list()

TFs_to_DEGs <- function(TF_names, GO_term, folder_name) {
  for (i in TF_names) {
    
    print(sprintf("Identifying DEGs associated with %s", i))
    
    interactions <- read_tsv(sprintf("String/%s/%s/%s_string_interactions.tsv", GO_term, folder_name, i))
    interactions <- interactions[interactions$`#node1` == sprintf("%s", i), c(1,2,9)]
    colnames(interactions) <- c("TF", "Gene_name", "Coexpression")
    interactions[nrow(interactions)+1, ] <- list(sprintf("%s", i), sprintf("%s", i), 1)
    
    GMCSF <- left_join(interactions, GMCSF_DEGs, by = "Gene_name")
    MCSF <- left_join(interactions, MCSF_DEGs, by = "Gene_name")
    
    combined <- data.frame(TF = GMCSF$TF, 
                           Gene_name = GMCSF$Gene_name, 
                           Coexpression = GMCSF$Coexpression,
                           GMCSF_L2FC = GMCSF$log2FoldChange, 
                           MCSF_L2FC = MCSF$log2FoldChange)
    
    # Remove genes which are not DE for both treatments
    combined <- combined[!with(combined,is.na(GMCSF_L2FC) & is.na(MCSF_L2FC)),]
    
    GMCSF_combined <- data.frame(Gene_name = combined$Gene_name,
                                 Treatment = c(rep("GMCSF", nrow(combined))),
                                 L2FC = combined$GMCSF_L2FC)
    
    MCSF_combined <- data.frame(Gene_name = combined$Gene_name,
                                Treatment = c(rep("MCSF", nrow(combined))),
                                L2FC = combined$MCSF_L2FC)
    
    combined2 <- rbind(GMCSF_combined, MCSF_combined)
    
    TFlist[[folder_name]][[i]] <- combined
    TFlist[["Network"]] <- rbind(TFlist[["Network"]], combined2)
    TFlist[["Full"]] <- rbind(TFlist[["Full"]], combined)
  }
  
  assign("TFlist", TFlist, envir = .GlobalEnv)
}

# Define the set of genes for examining protein-protein interaction - regulation of inflammatory response.
MCSF_TFs <- c("CLOCK", "RB1")
GMCSF_TFs <- c("NFKB1", "PPARD", "NR1D2")
Common_TFs <- c("BCL6", "CEBPB", "ESR1", "KLF4", "NFKBIZ", "NLRP3", "NR1D1", "NR1H3", "RELA", "SMAD3")

TFs_to_DEGs(TF_names = MCSF_TFs, 
            GO_term = "Regulation_of_inflammatory_response",
            folder_name = "MCSF")

TFs_to_DEGs(TF_names = GMCSF_TFs, 
            GO_term = "Regulation_of_inflammatory_response",
            folder_name = "GMCSF")

TFs_to_DEGs(TF_names = Common_TFs, 
            GO_term = "Regulation_of_inflammatory_response",
            folder_name = "Common")

# Produce a file for cytoscape network construction
cytoscape_data <- TFlist[["Network"]]

cytoscape_data <- unique(cytoscape_data)

Full_TFs <- c(MCSF_TFs, 
              #Common_TFs, 
              GMCSF_TFs)
Full_TFs <- data.frame(Gene_name = Full_TFs, Status = c(rep(1, length(Full_TFs))))

cytoscape_data <- left_join(cytoscape_data, Full_TFs, by = "Gene_name")

write.csv(cytoscape_data, "String/regualtion_of_inflammatory_response_TF_data.csv")

# Produce a table summarizing the transcription factors
GMCSF <- left_join(Full_TFs, GMCSF_DEGs, by = "Gene_name")
MCSF <- left_join(Full_TFs, MCSF_DEGs, by = "Gene_name")
TF_table <- data.frame(TF = GMCSF$Gene_name, MCSF_LFC = MCSF$log2FoldChange, GMCSF_LFC = GMCSF$log2FoldChange)
write.csv(TF_table, "String/regualtion_of_inflammatory_response_TF_table.csv")

# Define the set of genes for examining protein-protein interaction - regulation of inflammatory response.
MCSF_TFs <- c("SPI1", 
              #"RB1", 
              "RUNX1"#, 
              #"RARG", 
              #"ZBTB7A"
              )

GMCSF_TFs <- c("MYC", 
               #"MEIS1", 
               #"KLF13", 
               "MAFB"
               )

TFlist <- list()

TFs_to_DEGs(TF_names = MCSF_TFs, 
            GO_term = "Myeloid_cell_differentiation",
            folder_name = "MCSF")

TFs_to_DEGs(TF_names = GMCSF_TFs, 
            GO_term = "Myeloid_cell_differentiation",
            folder_name = "GMCSF")

# Produce a file for cytoscape network construction
cytoscape_data <- TFlist[["Network"]]

cytoscape_data <- unique(cytoscape_data)

Full_TFs <- c(MCSF_TFs, GMCSF_TFs)
Full_TFs <- data.frame(Gene_name = Full_TFs, Status = c(rep(1, length(Full_TFs))))

cytoscape_data <- left_join(cytoscape_data, Full_TFs, by = "Gene_name")

write.csv(cytoscape_data, "String/myeloid_cell_differentiation_TF_data.csv")

# Produce a table summarizing the transcription factors
GMCSF <- left_join(Full_TFs, GMCSF_DEGs, by = "Gene_name")
MCSF <- left_join(Full_TFs, MCSF_DEGs, by = "Gene_name")
TF_table <- data.frame(TF = GMCSF$Gene_name, MCSF_LFC = MCSF$log2FoldChange, GMCSF_LFC = GMCSF$log2FoldChange)
write.csv(TF_table, "String/myeloid_cell_differentiation_TF_table.csv")


#combined2$Status <- combined2$Gene_name == i
#combined2$Status <- ifelse(combined2$Status, "TF", NA)

##### Identify transcription factors within the 'Regulation of inflammatory response' GO term' #####
mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)
KEGG_packages <- c("clusterProfiler", 
                   "AnnotationDbi", 
                   "org.Hs.eg.db",
                   "R.utils", 
                   "biomaRt",
                   "dplyr",
                   "ggplot2",
                   "data.table",
                   "tidyverse")

suppressMessages(lapply(KEGG_packages, require, character.only = TRUE))

source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
generate_ensembl(mainDir)

# Read in GO results summary.
GO_results <- readRDS("clusterProfiler/GO_enrichment/GO_summary.rds")

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

#Read in transcription factor results
DETFs <- readRDS("DEG_lists/GMCSF_tfs_and_MCSF_tfs_comparison.rds")

# Convert ENSG ID to gene name
convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = DETFs$GMCSF_tfs$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DETFs$GMCSF_tfs$Gene_name <- object_final$external_gene_name

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = DETFs$MCSF_tfs$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DETFs$MCSF_tfs$Gene_name <- object_final$external_gene_name

GMCSF_TFs <- na.omit(left_join(df1, DETFs[["GMCSF_tfs"]], by = "Gene_name"))
MCSF_TFs <- na.omit(left_join(df2, DETFs[["MCSF_tfs"]], by = "Gene_name"))

GMCSF_shared_TFs <- semi_join(GMCSF_TFs, MCSF_TFs, by = "Gene_name")
MCSF_shared_TFs <- semi_join(MCSF_TFs, GMCSF_TFs, by = "Gene_name")

# Determine which DEGs have a different log2 fold change sign based on treatment
f <- function(a, b) 
{
  ifelse(a == 0 | b == 0, as.logical("FALSE"),!xor(sign(a)+1,sign(b)+1))
}
regulation <- f(a = GMCSF_shared_TFs$log2FoldChange, b = MCSF_shared_TFs$log2FoldChange)

# Identify unique TFs
GMCSF_unique_TFs <- data.frame(Gene_name = setdiff(GMCSF_TFs$Gene_name, GMCSF_shared_TFs$Gene_name))
GMCSF_unique_TFs <- left_join(GMCSF_unique_TFs, GMCSF_TFs, by = "Gene_name")

MCSF_unique_TFs <- data.frame(Gene_name = setdiff(MCSF_TFs$Gene_name, MCSF_shared_TFs$Gene_name))
MCSF_unique_TFs <- left_join(MCSF_unique_TFs, MCSF_TFs, by = "Gene_name")

shared <- data.frame(Gene_name = setdiff(GMCSF_TFs$Gene_name, GMCSF_unique_TFs$Gene_name))
shared <- left_join(shared, GMCSF_TFs, by = "Gene_name")
shared <- left_join(shared, MCSF_TFs, by = "Gene_name")
shared_TFs <- data.frame(Gene_name = shared$Gene_name, 
                         GMCSF_L2FC = shared$log2FoldChange.x,
                         MCSF_L2FC = shared$log2FoldChange.y,
                         L2FC_diff = abs(shared$log2FoldChange.x - shared$log2FoldChange.y))

# Plot the transcription factors as a heatmap
GMCSF <- data.frame(Gene = GMCSF_TFs$Gene_name, log2foldchange = GMCSF_TFs$log2FoldChange, 
                    Treatment = c(rep("GMCSF", nrow(GMCSF_TFs))))
MCSF <- data.frame(Gene = MCSF_TFs$Gene_name, log2foldchange = MCSF_TFs$log2FoldChange, 
                   Treatment = c(rep("MCSF", nrow(MCSF_TFs))))
heat_data <- rbind(GMCSF, MCSF)

# Produce the heatmap
heat_shared <- ggplot(heat_data, aes(Treatment, Gene, fill = log2foldchange)) + 
  geom_tile() + 
  ggtitle("Transcription factors DE within the regulation of inflammatory response BP") + 
  scale_fill_gradient2(low="steelblue3", mid = "grey90", high="coral") + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
heat_shared

# Write the plot to file
ggsave(sprintf("%s/plots/inflamm_response_TF_heatmap.png", mainDir), 
       bg = "white", 
       dpi = 300, 
       width = 6, 
       height = 8, 
       units = "in")

##### Identify transcription factors within the 'Myeloid cell differentiation' GO term' #####

mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)
KEGG_packages <- c("clusterProfiler", 
                   "AnnotationDbi", 
                   "org.Hs.eg.db",
                   "R.utils", 
                   "biomaRt",
                   "dplyr",
                   "ggplot2",
                   "data.table",
                   "tidyverse")

suppressMessages(lapply(KEGG_packages, require, character.only = TRUE))

source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
generate_ensembl(mainDir)

# Read in GO results summary.
GO_results <- readRDS("clusterProfiler/GO_enrichment/GO_summary.rds")

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

#Read in transcription factor results
DETFs <- readRDS("DEG_lists/GMCSF_tfs_and_MCSF_tfs_comparison.rds")

# Convert ENSG ID to gene name
convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = DETFs$GMCSF_tfs$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DETFs$GMCSF_tfs$Gene_name <- object_final$external_gene_name

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = DETFs$MCSF_tfs$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
DETFs$MCSF_tfs$Gene_name <- object_final$external_gene_name

GMCSF_TFs <- na.omit(left_join(df1, DETFs[["GMCSF_tfs"]], by = "Gene_name"))
MCSF_TFs <- na.omit(left_join(df2, DETFs[["MCSF_tfs"]], by = "Gene_name"))

GMCSF_shared_TFs <- semi_join(GMCSF_TFs, MCSF_TFs, by = "Gene_name")
MCSF_shared_TFs <- semi_join(MCSF_TFs, GMCSF_TFs, by = "Gene_name")

# Determine which DEGs have a different log2 fold change sign based on treatment
f <- function(a, b) 
{
  ifelse(a == 0 | b == 0, as.logical("FALSE"),!xor(sign(a)+1,sign(b)+1))
}
regulation <- f(a = GMCSF_shared_TFs$log2FoldChange, b = MCSF_shared_TFs$log2FoldChange)

# Identify unique TFs
GMCSF_unique_TFs <- data.frame(Gene_name = setdiff(GMCSF_TFs$Gene_name, GMCSF_shared_TFs$Gene_name))
GMCSF_unique_TFs <- left_join(GMCSF_unique_TFs, GMCSF_TFs, by = "Gene_name")

MCSF_unique_TFs <- data.frame(Gene_name = setdiff(MCSF_TFs$Gene_name, MCSF_shared_TFs$Gene_name))
MCSF_unique_TFs <- left_join(MCSF_unique_TFs, MCSF_TFs, by = "Gene_name")

# Plot the transcription factors as a heatmap
GMCSF <- data.frame(Gene = GMCSF_TFs$Gene_name, log2foldchange = GMCSF_TFs$log2FoldChange, 
                    Treatment = c(rep("GMCSF", nrow(GMCSF_TFs))))
MCSF <- data.frame(Gene = MCSF_TFs$Gene_name, log2foldchange = MCSF_TFs$log2FoldChange, 
                   Treatment = c(rep("MCSF", nrow(MCSF_TFs))))
heat_data <- rbind(GMCSF, MCSF)

# Produce the heatmaps
heat_shared <- ggplot(heat_data, aes(Treatment, Gene, fill = log2foldchange)) + 
  geom_tile() + 
  ggtitle("Transcription factors DE within the myeloid cell differentiation BP") + 
  scale_fill_gradient2(low="steelblue3", mid = "grey90", high="coral") + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
heat_shared

# Write the plot to file
ggsave(sprintf("%s/plots/myeloid_diff_TF_heatmap.png", mainDir), 
       bg = "white", 
       dpi = 300, 
       width = 6, 
       height = 10, 
       units = "in")
