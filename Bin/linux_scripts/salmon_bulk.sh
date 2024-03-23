#! /bin/bash

echo "The following fastq files have been detected nested within the 'Matthew_Masters/Raw/Data/bulkseq_data' directory."
find /mnt/c/Users/Matthew/Documents/Matthew_Masters/Raw/Data/bulkseq_data/ -type f -name "*.fastq.gz"

salmon quant 
	--libType A 
	-1 filex/filex_1.fastq.gz -2 filex/filex_2.fastq.gz 
	-i /mnt/c/Users/Matthew/Documents/Matthew_Masters/Raw/Data/bulkseq_data/bulk_index 
	-o /mnt/c/Users/Matthew/Documents/Matthew_Masters/Work/salmon/filex 
	--gcBias 
	--seqBias 
	--posBias 
	--validateMappings

echo "This script is incomplete. I need a way to call each fair of files, and the method I used in the fastqc_bulk.sh script wont work because pairs of files must be called together"
