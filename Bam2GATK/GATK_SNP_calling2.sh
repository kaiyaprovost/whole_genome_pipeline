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
echo Calling SNPs with GATK Haplotype Caller
echo
echo

time java -jar /home/lmoreira/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller \
-R $ref \
-I $bam \
-o $name.raw.g.vcf \
--emitRefConfidence GVCF \
-minPruning 1 \
-minDanglingBranchLength 1 \
-nct 16
