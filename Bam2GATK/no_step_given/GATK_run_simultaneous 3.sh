#!/bin/bash

for filename in *.bam; do
	qsub -v bam=$filename $1
done
