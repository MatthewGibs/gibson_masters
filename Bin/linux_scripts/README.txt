This directory contains six scripts as detailed below:

Scripts for running salmon, fastqc, and multiqc for the bulkseq data.

Please note: Do not open these files manually through windows, as it can introduce invisible characters into a .sh script that prevents it from functioning correctly. Furthermore, certain files such as the fastq data and bulk_index have been omitted from the github due to file size constraints.

salmon_bulk.sh: Performs salmon quantification on pairs of fastq files housed in directories with a corresponding name. Each of these directories should be located within the 'Raw/bulkseq_data' directory. The index ('bulk_index') should likewise be located in this directory. The output for each file pair is written to the 'Work/salmon' directory.

fastqc_bulk.sh: Detects any fastq files present in the Raw/bulkseq_data directory and performs FastQC. Outputs are written to the Work/FastQC directory.

multiqc_bulk.sh: Performs multiqc on all zipped fastqc output files in the Matthew_Masters/Work/FastQC directory. Outputs are written to the Work/FastQC/MultiQC directory.