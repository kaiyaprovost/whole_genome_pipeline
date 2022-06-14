#!/bin/bash

#for i in /vz-nas1-active/ProcessedGenomicReads/FILTERED/fus*/DEDUP/ ; do
#for i in /vz-nas1-active/RawGenomicReads/FASTQ_SCREEN_RESULTS/Parrots/BAM_dedup/*bam; do

#echo $i;
#cd $i;

for i in F*dedup.bam; do
echo $i;

findname=$(echo $i | cut -f 1 -d '.' | cut -f 3 -d "/")
replace=$(echo $i | cut -d "_" -f 2-5 | cut -f 3 -d "/"); 
if [ ! -f ${findname}*reheadered.bam ] && [ ! -f ../REHEAD*/${findname}*reheadered.bam ]; then
echo "Reheadering"
~/nas3/samtools/samtools view -H $i | sed -e "s/SM:${findname}/SM:${replace}/" | ~/nas3/samtools/samtools reheader - $i > $i.dedup.reheadered.bam
else
echo "Bam already reheadered"
fi

done; 

mv -v *rehead*bam ../REHEAD*/

#cd ~ ;
#done; 


