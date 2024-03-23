This directory contains six scripts as detailed below:

bulkseq_script.R: A master script that calls the tximport_script.R, DESeq2_script.R, and clusterProfiler_bulk.R. The script requires that the user specify salmon 'quant.sf' files as input, along with multiple metadata fields.

tximport_script.R: This script runs Tximport using salmon output files and a transcripts-to-gene file. The script also performs annotation by calling the biomaRt_annotate.R script. The output is saved to file as 'txi.count.matrix.rds'.

DESeq2_script.R: This script runs DESeq2 on the output of the tximport_script.R, which must be specified when calling the 'DESeq2_DEA' function. The user must also specify a column of the metadata object which is to be used for the experimental design. If multiple factor levels are present, multiple files will be output based on said levels.

clusterProfiler_bulk.R: Performs GO enrichment by calling the 'GO_bulk.R' script. As default, this script subsets the list of significantly differentially expressed genes to perform the enrichment of only upregulated genes, only downregulated genes, and lastly all genes. The output is therefore 3 pairs of files (a .rds and .csv) for each file that was output from the prior DESeq2 script.