This directory contains twelve scripts, each of which are detailed below. Scripts are detailed in the order in which they are used within the data pipeline.

bulkseq_master.R: A master script that calls the following scripts: 
- tximport_script.R
- gen_metadata.R
- DESeq2_script.R
- gen_DEGs.R
- identify_TFs.R
- clusterProfiler_GO.R
- gen_GO_summary.R
- tf_interactions.R
- dissertation_plots.R
The script requires that the user specify paths to the 'Bin' and 'Work' directories which contain scripts and key input and outputs respectively.

tximport_script.R: This script runs Tximport using salmon output files and a transcripts-to-gene file. The script also performs annotation by calling the biomaRt_annotate.R script. The output is saved to file as 'txi.count.matrix.rds'.

biomarRt_annotate.R: This script contains three functions; generate_ensembl, convert_ID, and gather_attributes. generate_ensembl searches for and reads in (or creates if it does not exist) the ensembl_annotation.rds file. convert_ID takes two attribute types as input, as well as a vector of the first ID type, then converts the initial ID to matches of the final ID. The gather_attributes function is similar to convert_ID, although it allows multiple attributes to be specified.

gen_metadata.R: Creates the metadata file specific to the nine pairs of FastQ files used in the construction of the data pipeline, and saves it to file as 'sample_metadata.rds'.

DESeq2_script.R: This script runs DESeq2 on the output of the tximport_script.R, which must be specified when calling the 'DESeq2_DEA' function. The user must also specify a column of the metadata object which is to be used for the experimental design. If multiple factor levels are present, multiple files will be output based on said levels.

gen_DEGs.R: This script takes in an input of two lists of DEGs and creates three additional lists (common and unique to each input list). The outputs are all saved together in a .rds file.

identify_TFs.R: This script contains a single similarly named function which uses biomaRt to identify potential transcription factors from an imput list of DEGs, saving the output list to file.

clusterProfiler_GO.R: This script contains two functions; cluster_Profiler_GO and compare_cluster_GO. cluster_Profiler_GO performs individual GO enrichment (enrichGO) on each set of genes given in a list format using the 'GO_bulk' subscript. compare_cluster_GO performs collective GO enrichment (compareCluster) on a list containing two or more gene sets.

GO_bulk.R: This script contains the GO_enrich function, which subsets a supplied list of DEGs to perform enrichment on only upregulated genes, only downregulated genes, and lastly all genes.

gen_GO_summary.R: This script contains the 'summarise_GO' function, which collates information from various GO enrichment results into a single list object (GO_summary.rds).

tf_interactions.R: This script contains the 'collate_tf_DEGs' function, which produces two files of transcription factors and their protein-protein interacting partners which are also differentially expressed. The script is specifically used for a set of nine transcritpion factors belonging to the 'myeloid cell differentiation' and 'regulation of inflammatory reponse' GO terms.

dissertation_plots.R: This script contains a number of functions to produce various plots using outputs of all the aforementioned scripts.
