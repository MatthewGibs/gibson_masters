#! /bin/bash

echo "Running salmon on all fastq files within the 'Raw/bulkseq_data' directory."

for fn in /mnt/c/Users/Matthew/Documents/git/gibson_masters/Raw/bulkseq_data/*/;
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant \
	--libType A \
	-1 ${fn}/${samp}_1.fastq.gz \
	-2 ${fn}/${samp}_2.fastq.gz \
	-i /mnt/c/Users/Matthew/Documents/git/gibson_masters/Raw/bulkseq_data/bulk_index \
	-o /mnt/c/Users/Matthew/Documents/git/gibson_masters/Work/salmon/${samp} \
	--gcBias \
	--seqBias \
	--posBias \
	--validateMappings

echo "Salmon quantification successfully completed."

done
