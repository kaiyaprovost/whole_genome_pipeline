#!/bin/bash

spp="Parrots"

#cd "/Users/kprovost/Documents/folder_for_popgenome/"
cd /vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/

i=${spp}-called.geno.vcf

#for i in ${spp}*vcf; do
echo $i 
j=${i%.vcf}
mkdir ./${j}

#/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java \
#-jar /Users/kprovost/Documents/folder_for_popgenome/snpEff_latest_core/snpEff/SnpSift.jar split $i

java -Djava.io.tmpdir=/home/kprovost/nas3/tmp/ -jar /home/kprovost/nas3/ANGSD_pipeline/SnpSift.jar split $i

gzip $i

mv ./${j}*.vcf ./${j}

done; 