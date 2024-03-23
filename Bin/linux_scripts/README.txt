This directory contains six scripts as detailed below:

Scripts for running salmon, fastqc, and multiqc for the bulkseq data.

Please note: Do not open these files manually through windows, as it can introduce invisible characters into a .sh script that prevents it from functioning correctly.

salmon_bulk.sh: Currently incomplete, but needs to be able to construct the index and then run any fastq files through salmon.

fastqc_bulk.sh: Detects any fastq files present in the Raw/Data/bulkseq_data directory and performs FastQC. Outputs are written to the Matthew_Masters/Work/FastQC directory.

multiqc_bulk.sh: Performs multiqc on all fasqc files in the Matthew_Masters/Work/FastQC directory. Outputs are written to the Matthew_Masters/Work/FastQC/MultiQC directory.