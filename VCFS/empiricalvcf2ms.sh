#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=99999:00:00
#PBS -N vcf2temp
#PBS -j oe
#PBS -m ae
#PBS -M kprovost@amnh.org
#PBS -k oe

#!/bin/bash

cd "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/"

for folder in */; do 
echo $folder; 
cd $folder

for vcffile in *vcf; do 

echo "#####
starting vcf conversion"

echo $vcffile

prefix=${vcffile%.vcf}
fulltempfile=${prefix}.fulltemp
outSamp=25

if [ -f "$fulltempfile" ]
then
echo "$fulltempfile found."
else
echo "$fulltempfile not found."
perl "/home/kprovost/nas3/vcf2MS.pl" $vcffile $fulltempfile $outSamp
fi

done; 
cd "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/"
done;