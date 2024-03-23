##### ClusterProfiler plots working area ######

###### Setup ######
mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)

clusterProfiler_plots_packages <- c("data.table", "clusterProfiler", "ggplot2", "tidyverse", "dplyr", "org.Hs.eg.db", "enrichplot", "forcats")

suppressMessages(lapply(clusterProfiler_plots_packages, require, character.only = TRUE))

# Creates the output directory for this script
output_dir <- "plots/clusterProfiler_plots"
if (!dir.exists(output_dir)){
   dir.create(output_dir)
   print("Creating 'clusterProfiler_plots' directory.")
} else {}
rm(output_dir)

##### General Plot testing #####
examineres_GMCSF_unique <- readRDS("clusterProfiler_bulk/GMCSF_unique_enrichment/GO_BOTH_BP.rds")
examineres_MCSF_unique <- readRDS("clusterProfiler_bulk/MCSF_unique_enrichment/GO_BOTH_BP.rds")
examineres_common <- readRDS("clusterProfiler_bulk/GMCSF_MCSF_common_enrichment/GO_BOTH_BP.rds")

examineres_GMCSF <- readRDS("~/Matthew_Masters/Work/clusterProfiler_bulk/treatment_GMCSF_vs_none_enrichment/GO_BOTH_BP.rds")
examineres_MCSF <- readRDS("~/Matthew_Masters/Work/clusterProfiler_bulk/treatment_MCSF_vs_none_enrichment/GO_BOTH_BP.rds")

goplot(examineres_MCSF, showCategory = 20)
goplot(examineres_common)

barplot(examineres_GMCSF, showCategory = 15)
barplot(examineres_MCSF, showCategory = 15)

dotplot(examineres_GMCSF, showCategory = 10) + ggtitle("dotplot for ORA of GMCSF")
dotplot(examineres_MCSF, showCategory = 10) + ggtitle("dotplot for ORA of MCSF")

cnetplot(examineres_MCSF)
edox <- setReadable(examineres_MCSF, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)




##### Produce a results table of common GO terms #####

GO_results <- readRDS("clusterProfiler/GO_enrichment/GO_summary.rds")

GMCSF_shared <- semi_join(GO_results[["GMCSF"]][["GMCSF_both"]], GO_results[["MCSF"]][["MCSF_both"]], by = "Description")
MCSF_shared <- semi_join(GO_results[["MCSF"]][["MCSF_both"]], GO_results[["GMCSF"]][["GMCSF_both"]], by = "Description")

df_combined <- left_join(x = GMCSF_shared, y = MCSF_shared, by = "Description")
df_combined <- df_combined[, c(2,3,10,19)]
colnames(df_combined) <- c("GO ID", "GO term description", "Number of GM-CSF DEGs", "Number of M-CSF DEGs")
df_combined <- df_combined[order(df_combined$`GO term description`), ]

# Next we must find how many gene are shared or unique across the common GO terms.
common_term_names <- df_combined$`GO ID`
Common_DEG_number <- c()
GMCSF_unique_number <- c()
MCSF_unique_number <- c()

for (i in common_term_names) {
   # Extract a specific column from a dataframe based on  a particular character in a specified column
   df1 <- GMCSF_shared[GMCSF_shared$ID %like% sprintf("%s", i), 9]
   df1 <- data.frame(str_split(df1, "/"))
   colnames(df1) <- "Gene_name"
   
   df2 <- MCSF_shared[MCSF_shared$ID %like% sprintf("%s", i), 9]
   df2 <- data.frame(str_split(df2, "/"))
   colnames(df2) <- "Gene_name"
   
   common_genes <- length(intersect(df1$Gene_name, df2$Gene_name))
   df1_unique <- (nrow(df1)-common_genes)
   df2_unique <- (nrow(df2)-common_genes)
   
   Common_DEG_number <- append(Common_DEG_number, common_genes)
   GMCSF_unique_number <- append(GMCSF_unique_number, df1_unique)
   MCSF_unique_number <- append(MCSF_unique_number, df2_unique)
}

df_combined <- data.frame(df_combined$`GO ID`, df_combined$`GO term description`, 
                          df_combined$`Number of GM-CSF DEGs`, df_combined$`Number of M-CSF DEGs`, 
                          Common_DEG_number, GMCSF_unique_number, MCSF_unique_number)
colnames(df_combined) <- c("GO ID", "GO term description", "Number of GM-CSF DEGs", "Number of M-CSF DEGs", "Number of common DEGs", "GM-CSF unique", "M-CSF unique")

# Creates the output directory for this script
output_dir <- "clusterProfiler_bulk/dissertation_tables"
if (!dir.exists(output_dir)){
   dir.create(output_dir)
   print("Creating 'dissertation_tables' directory.")
} else {}

write.csv(df_combined, file = sprintf("%s/%s/Common_GO_terms.csv", mainDir, output_dir))

##### cnetplots for ClusterProfiler #####
mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)

GO_results <- readRDS("clusterProfiler/GO_enrichment/GO_summary.rds")

# Read in lists of DEGs
DEGs <- readRDS("DEG_lists/GMCSF_and_MCSF_comparison.rds")

# Read in DEGs
GMCSF_DEGs <- DEGs$GMCSF
MCSF_DEGs <- DEGs$MCSF

source("~/Matthew_Masters/Bin/bulkseq_scripts/biomaRt_annotate.R")
generate_ensembl(mainDir)

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = GMCSF_DEGs$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
GMCSF_DEGs[,1] <- object_final[,2]
GMCSF_geneList <- GMCSF_DEGs$log2FoldChange
names(GMCSF_geneList) <- GMCSF_DEGs$Gene_name

convert_ID(initial_ID = "ensembl_gene_id", # This is the starting ID that you wish to change to another type.
           final_ID = "external_gene_name", # The ID type you want in your final object.
           object = MCSF_DEGs$Gene_name # character vector containing the 'initial ID' type that you want to convert
)
MCSF_DEGs[,1] <- object_final[,2]
MCSF_geneList <- MCSF_DEGs$log2FoldChange
names(MCSF_geneList) <- MCSF_DEGs$Gene_name

# Read in enrichGO results
examineres_GMCSF <- readRDS("~/Matthew_Masters/Work/clusterProfiler/GO_enrichment/GMCSF/GO_BOTH.rds")
examineres_MCSF <- readRDS("~/Matthew_Masters/Work/clusterProfiler/GO_enrichment/MCSF/GO_BOTH.rds")

# Specify hand-picked GO terms of interest
categories_MCSF <- c("leukocyte cell-cell adhesion", "mononuclear cell differentiation", "negative regulation of immune system process")
#OR
categories_GMCSF <- c("chromosome segregation", "nuclear division", "regulation of mitotic cell cycle")
#OR
categories <- c("regulation of inflammatory response", "myeloid cell differentiation")

p1 <- cnetplot(examineres_GMCSF, color.params = list(foldChange = GMCSF_geneList, category = 'firebrick'), showCategory = categories)
p2 <- cnetplot(examineres_MCSF, color.params = list(foldChange = MCSF_geneList, category = 'firebrick'), showCategory = categories)
cowplot::plot_grid(p1, p2, ncol=2, labels = c("GM-CSF", "M-CSF"), label_size = 30)

ggsave(sprintf("%s/plots/clusterProfiler_plots/test.png", mainDir), 
       bg = "white",
       dpi = 300,
       width = 20,
       height = 12,
       units = "in")

# OR call the handy function below

# Read in function for doing plots
source("~/Matthew_Masters/Bin/bulkseq_scripts/bulk_plots.R")

# Produce plot
clusterProfiler_cnetplot(mainDir = "~/Matthew_Masters/Work")


##### experimenting area #####

### test emmaplots

# Read in the data
GMCSF_unique <- readRDS("clusterProfiler_bulk/GMCSF_unique_enrichment/GO_BOTH_BP.rds")
MCSF_unique <- readRDS("clusterProfiler_bulk/MCSF_unique_enrichment/GO_BOTH_BP.rds")

# Establish pairwise termism
GMCSF <- pairwise_termsim(GMCSF_unique)
MCSF <- pairwise_termsim(MCSF_unique)

p1 <- emapplot(GMCSF, showcategory = 30, fill = "p.adjust") + scale_color_viridis()
p2 <- emapplot(MCSF, showcategory = 30) + scale_color_gradient(low = "#56B1F7", high = "#132B43")

cowplot::plot_grid(p1, p2, nrow = 2, labels = c("GM-CSF", "M-CSF"), label_size = 30)

### Treeplot
edox2 <- pairwise_termsim(examineres_GMCSF)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

x <- enrichDO(de)
x2 <- pairwise_termsim(x)
cnetplot(x2)
# use `layout` to change the layout of map
cnetplot(x2, layout = "star")
# use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
cnetplot(x2, showCategory = 10)

cnetplot(x2, showCategory = categorys)


p1 <- heatplot(examineres_GMCSF, showCategory=3)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

##### emmapplots #####
examineres_GMCSF <- pairwise_termsim(examineres_GMCSF)
examineres_MCSF <- pairwise_termsim(examineres_MCSF)

p1 <- emapplot(examineres_GMCSF, showcategory = 30)
p2 <- emapplot(examineres_MCSF, showcategory = 30)

cowplot::plot_grid(p1, p2, nrow = 2, labels = c("GM-CSF", "M-CSF"), label_size = 30)
ggsave(sprintf("%s/plots/clusterProfiler_plots/emmapplots.png", mainDir), 
       bg = "white",
       dpi = 400,
       width = 13,
       height = 20,
       units = "in")

# OR call the handy function below

# Read in function for doing plots
source("~/Matthew_Masters/Bin/bulkseq_scripts/bulk_plots.R")

# Produce plot
clusterProfiler_emapplot(mainDir = "~/Matthew_Masters/Work")

##### Examine similar GO terms between 2 lists #####

mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)

GO_results <- readRDS("clusterProfiler/GO_enrichment/GO_summary.rds")
df1 <- GO_results[["GMCSF"]][["GMCSF_both"]]
df2 <- GO_results[["MCSF"]][["MCSF_both"]]

df1_shared <- semi_join(df1, df2, by = "Description")
df2_shared <- semi_join(df2, df1, by = "Description")

df_combined <- left_join(x = df1_shared, y = df2_shared, by = "Description")
df_combined <- df_combined[, c(2,3,10,19)]
colnames(df_combined) <- c("GO ID", "GO term description", "Number of df1 DEGs", "Number of df2 DEGs")
df_combined <- df_combined[order(df_combined$`GO term description`), ]

GMCSF_up <- GO_results[["GMCSF"]][["GMCSF_up"]]
MCSF_up <- GO_results[["MCSF"]][["MCSF_up"]]
up_shared <- semi_join(GMCSF_up, MCSF_up, by = "Description")

GMCSF_down <- GO_results[["GMCSF"]][["GMCSF_down"]]
MCSF_down <- GO_results[["MCSF"]][["MCSF_down"]]
down_shared <- semi_join(GMCSF_down, MCSF_down, by = "Description")

##### KEGG exploration #####
mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)

KEGG_packages <- c("clusterProfiler", 
                   "AnnotationDbi", 
                   "org.Hs.eg.db",
                   "R.utils", 
                   "biomaRt",
                   "dplyr",
                   "pathview")

suppressMessages(lapply(KEGG_packages, require, character.only = TRUE))

GMCSF <- read.csv("clusterProfiler/KEGG_enrichment/GM-CSF_sig_KEGG.csv")
MCSF <- read.csv("clusterProfiler/KEGG_enrichment/M-CSF_sig_KEGG.csv")
Common <- read.csv("clusterProfiler/KEGG_enrichment/Common_sig_KEGG.csv")
GMCSF_unique <- read.csv("clusterProfiler/KEGG_enrichment/GM-CSF_unique_sig_KEGG.csv")
MCSF_unique <- read.csv("clusterProfiler/KEGG_enrichment/M-CSF_unique_sig_KEGG.csv")

shared_KEGG <- semi_join(GMCSF, MCSF, by = "ID")


##### Exploration of 'Regualtion of inflammatory response' GO term' #####

mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)
KEGG_packages <- c("clusterProfiler", 
                   "AnnotationDbi", 
                   "org.Hs.eg.db",
                   "R.utils", 
                   "biomaRt",
                   "dplyr",
                   "ggplot2")

suppressMessages(lapply(KEGG_packages, require, character.only = TRUE))

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

df1 <- mapIds(org.Hs.eg.db,
               keys = df1$Gene_name, #Column containing Ensembl gene ids
               column="ENSEMBL",
               keytype="SYMBOL",
               multiVals="first")
df1 <- as.data.frame(df1)
colnames(df1) <- "Gene_name"

df2 <- mapIds(org.Hs.eg.db,
              keys = df2$Gene_name, #Column containing Ensembl gene ids
              column="ENSEMBL",
              keytype="SYMBOL",
              multiVals="first")
df2 <- as.data.frame(df2)
colnames(df2) <- "Gene_name"

DEG_list <- list(df1, df2)
names(DEG_list) <- c("GM-CSF_both_inflamreg", "M-CSF_both_inflamreg")

shared_inflamm_genes <- semi_join(df1, df2, by = "Gene_name")
GMCSF_inflamm_genes <- setdiff(df1$Gene_name, shared_inflamm_genes$Gene_name)
MCSF_inflamm_genes <- setdiff(df2$Gene_name, shared_inflamm_genes$Gene_name)

# This function utilizes clusterProfiler's 'enrichKEGG' function to perform GO enrichment on individual sets of genes.
cluster_Profiler_KEGG(mainDir, # The path to the working directory
                      DESeq2_output, # The path to the .rds files output from the DESeq2 script
                      DEG_list # An object of the 'list' type containing dataframes with a column named 'Gene_name'. The Gene_name column should contain gene names in the ENSG ID format.
)

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
df1 <- df1[df1$ID %like% "GO:0050727", 9] # "GO:0050900" = leukocyte migration, "GO:0050727" = regulation of inflammatory response
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

# Plot the transcription factors as a heatmap
GMCSF <- data.frame(Gene = GMCSF_TFs$Gene_name, log2foldchange = GMCSF_TFs$log2FoldChange, 
                    Treatment = c(rep("GMCSF", nrow(GMCSF_TFs))))
MCSF <- data.frame(Gene = MCSF_TFs$Gene_name, log2foldchange = MCSF_TFs$log2FoldChange, 
                    Treatment = c(rep("MCSF", nrow(MCSF_TFs))))
heat_data <- rbind(GMCSF, MCSF)

# Produce the heatmap
heat_shared <- ggplot(heat_data, aes(Treatment, Gene, fill = log2foldchange)) + 
  geom_tile() + 
  ggtitle("Common Genes") + 
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
       width = 5, 
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
  ggtitle("Common Genes") + 
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
       width = 5, 
       height = 10, 
       units = "in")

##### Heatmap plot of a selected GO term #####
mainDir <- "~/Matthew_Masters/Work"
setwd(mainDir)
packages <- c("clusterProfiler", 
                   "AnnotationDbi", 
                   "org.Hs.eg.db",
                   "R.utils", 
                   "biomaRt",
                   "dplyr",
                   "ggplot2",
                   "data.table",
                   "tidyverse", 
                   "enrichplot", 
                   "forcats")

suppressMessages(lapply(packages, require, character.only = TRUE))

GO_heatmap <- function(GO_ID, plot_name) {
  
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
  GMCSF_unique <- data.frame(Gene_name = GMCSF_unique)
  GMCSF_unique <- left_join(GMCSF_unique, GMCSF_DEGs, by = "Gene_name")
  GMCSF_unique <- data.frame(rep("GMCSF", nrow(GMCSF_unique)), GMCSF_unique$Gene_name, GMCSF_unique$log2FoldChange)
  colnames(GMCSF_unique) <- c("Treatment", "Gene", "LFC")
  
  MCSF_unique <- setdiff(MCSF_GO$Gene_name, GMCSF_GO$Gene_name)
  MCSF_unique <- data.frame(Gene_name = MCSF_unique)
  MCSF_unique <- left_join(MCSF_unique, MCSF_DEGs, by = "Gene_name")
  MCSF_unique <- data.frame(rep("MCSF", nrow(MCSF_unique)), MCSF_unique$Gene_name, MCSF_unique$log2FoldChange)
  colnames(MCSF_unique) <- c("Treatment", "Gene", "LFC")
  
  # common genes
  data_shared <- bind_rows(GMCSF_common, MCSF_common)
  
  # Unique genes
  data_unique <- bind_rows(GMCSF_unique, MCSF_unique)
  
  # Produce the heatmaps
  heat_shared <- ggplot(data_shared, aes(Treatment, Gene, fill = LFC)) + 
    geom_tile() + 
    ggtitle("A") + 
    scale_fill_gradient2(low="steelblue3", mid = "grey90", high="coral") + 
    theme(plot.title = element_text(hjust = 0.5), 
          plot.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank())
  
  heat_unique <- ggplot(data_unique, aes(Treatment, Gene, fill = LFC)) + 
    geom_tile() + 
    ggtitle("B") + 
    scale_fill_gradient2(low="steelblue3", mid = "grey90", high="coral") + 
    theme(plot.title = element_text(hjust = 0.5), 
          plot.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank())
  
  cowplot::plot_grid(heat_shared, heat_unique, ncol=2)
  
  ggsave(filename = sprintf("%s/plots/%s.png", mainDir, plot_name),  # name of the image file
         dpi = 300,
         width = 8,
         height = 12,
         units = "in")
}

GO_heatmap(GO_ID = "GO:0050727", plot_name = "inflamm_heatmap")
GO_heatmap(GO_ID = "GO:0030099", plot_name = "myeloid_diff_heatmap")
GO_heatmap(GO_ID = "GO:0050900", plot_name = "leukocyte_migration_heatmap")
