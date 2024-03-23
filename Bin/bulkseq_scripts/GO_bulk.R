# Perform GO enrichment test
GO_enrich <- function(DEGs, # Data frame of differentially expressed genes.
                      background_genes,
                      Regulation, # choose between "UP", "DOWN", or "BOTH".
                      GO_Term, # choose between "BP" (biological process), "CC" (cellular component) or "MF" (molecular function).
                      output_dir)
{
  ##### Perform GO enrichment #####
  # Determines if up-regulated or down-regulated genes will be used for analysis. "BOTH" will prevent sub-setting the gene list.
  if(Regulation == "UP"){
    sub_data <- subset(DEGs, `log2FoldChange` > 0)
    print("Performing GO enrichment on upregulated genes.")
  }
  
  if(Regulation == "DOWN"){
    sub_data <- subset(DEGs, `log2FoldChange` < 0)
    print("Performing GO enrichment on downregulated genes.")
  }
  
  if(Regulation == "BOTH"){
    sub_data <- DEGs
    print("Performing GO enrichment on all differentially expressed genes.")
  }
  
  # Pulls out gene ids for all differentially expressed genes.
  hgnc_ids <- sub_data$Gene_name
  
  # Perform the GO enrichment analysis.
  DEGS_enriched <- enrichGO(gene = hgnc_ids, #Specify the gene set 
                            universe = background_genes, 	#Universe refers to all the genes that it will be compared to. This would be a list of all the genes in Entrez IDs from either the annotation file, alternatively may be all the genes found on the probe.
                            OrgDb = org.Hs.eg.db,	#The gene set must be annotated with the GO terms and requires a database to do so against. The correct database for the organsims which you're using is needed.
                            keyType = "SYMBOL",
                            ont = GO_Term, 				#This specifies what ontology you wish to do enrichment on (BP= Biological process, CC= cellular component and MF= Molecular function)
                            pAdjustMethod = "BH",	# this specifies the statistical method to perform on the genes. Benjamini and hochburg (BH) is prefered for GOEA (GO enrichment analysis)
                            pvalueCutoff = 0.01,		# Cut off point for significant genes
                            qvalueCutoff = 0.05, 	
                            readable = FALSE)
  
  # Call the simplify function to reduce redundancy in the final list of GO terms.
  DEGS_enriched <- clusterProfiler::simplify(DEGS_enriched, cutoff = 0.7, by = "p.adjust")
  
  # Write the enrichment object file.
  saveRDS(DEGS_enriched, file = sprintf("%s/GO_%s.rds", output_dir, Regulation))
  
  # Subset the enrichment object for the results and save them to file.
  DEGs_enriched_result <- DEGS_enriched@result
  DEGs_enriched_result <- subset(DEGs_enriched_result, `p.adjust` < 0.01)
  write.csv(DEGs_enriched_result, sprintf("%s/GO_%s_significant.csv", output_dir, Regulation))
  
  ##### Produce a basic lolipop plot to represent top GO terms #####
  
  p1 <- ggplot(DEGS_enriched, showCategory = 15, aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
    geom_segment(aes(xend = 0, yend = Description)) +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_viridis_c(guide = guide_colorbar(reverse=TRUE)) +
    scale_size_continuous(range = c(1, 7)) +
    theme_minimal() + 
    xlab("Gene Ratio") +
    ylab(NULL) + 
    ggtitle(sprintf("%s", output_dir))
  
  # Save the plot data
  saveRDS(p1, sprintf("%s/GO_%s_loliplot_data.rds", output_dir, Regulation))
  
  # Save the plot to file
  ggsave(sprintf("%s/GO_%s_loliplot.png", output_dir, Regulation), 
         bg = "white",
         dpi = 300)
  
  # Clean memory
  rm(DEGS_enriched)
  rm(p1)
}