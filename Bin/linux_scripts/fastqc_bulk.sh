#! /bin/bash

echo "The following fastq files have been identified in the Raw/Data/bulkseq_data"
find /mnt/c/Users/Matthew/Documents/Matthew_Masters/Raw/Data/bulkseq_data/ -type f -name "*.fastq.gz"

echo "Creating the output directory for the FastQC analysis. It can be found at 'Matthew_Masters/Work/FastQC'"
mkdir /mnt/c/Users/Matthew/Documents/Matthew_Masters/Work/FastQC

echo "Performing fastQC on all fastq files. This may take some time."
find /mnt/c/Users/Matthew/Documents/Matthew_Masters/Raw/Data/bulkseq_data/ -type f -name "*.fastq.gz" | xargs fastqc -o /mnt/c/Users/Matthew/Documents/Matthew_Masters/Work/FastQC

echo "FastQC complete for all of the given fastq files."
