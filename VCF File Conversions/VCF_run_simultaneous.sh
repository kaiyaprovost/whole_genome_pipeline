#!/bin/bash

for filename in *.vcf; do
echo $filename;
qsub -v vcf=$filename $1;
done;
