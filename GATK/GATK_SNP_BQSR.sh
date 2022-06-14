#!/bin/bash

# change to the working directory
cd $PBS_O_WORKDIR
echo "pbsworkdir"
echo $PBS_O_WORKDIR
EXECDIR=`pwd`
export PATH=./:$PATH
echo $PATH

#Arguments:
## ref=reference sequence
## bam=bam file

name=$(echo $bam | cut -f 1 -d '.')

echo
echo "#######################"
echo $name
echo "#######################"

echo
echo
echo BQSR
echo
echo

time java -jar /home/lmoreira/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T BaseRecalibrator \
-R $ref \
-I $bam \
-knownSites combined.db.vcf \
-o $name.recal_data.grp \
-nct 32

echo
echo
echo Producing recalibrated VCF file
echo
echo

time java -jar /home/lmoreira/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller \
-R $ref \
-I $bam \
-o $name.recal.raw.g.vcf \
--emitRefConfidence GVCF \
-minPruning 1 \
-minDanglingBranchLength 1 \
-BQSR $name.recal_data.grp \
-nct 32
