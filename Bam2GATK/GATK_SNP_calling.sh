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


date
time

for i in */; do
	cd $i
	name="${i///}"
	bam=$name.dedup.bam
	
	echo
	echo "#######################"
	echo $name
	echo "#######################"

	echo
	echo
	echo Calling SNPs with GATK Haplotype Caller
	echo
	echo

	time gatk HaplotypeCaller \
	-R $ref \
	-I $bam \
	-ERC GVCF \
	-O $name.raw.g.vcf \
	--min-pruning 1 \
	--min-dangling-branch-length 1 \

done

date
time
