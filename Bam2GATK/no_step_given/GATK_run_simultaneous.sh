#!/bin/bash

for filename in *merged.dedup.bam; do

echo $filename;

## submit the bam file to the script
qsub -v bam=$filename $1;

done;
