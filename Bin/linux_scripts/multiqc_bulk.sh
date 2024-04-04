#! /bin/bash

echo "Performing MultiQC on all FastQC output files in the 'Matthew_Masters/Work/FastQC' directory."

multiqc /mnt/c/Users/Matthew/Documents/git/gibson_masters/Work/FastQC -o /mnt/c/Users/Matthew/Documents/git/gibson_mastersWork/FastQC/MultiQC
