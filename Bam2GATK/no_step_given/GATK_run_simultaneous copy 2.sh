#!/bin/bash


for filename in BAM_files/*.dedup.bam; do
	qsub -v ref=/home/lmoreira/nas3/Picoides_pubescens_reference_genome/Picoides_pubescens_ref-sorted.fa,bam=$filename $1
done
