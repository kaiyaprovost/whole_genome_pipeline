#!/bin/bash

# change to the working directory
cd $PBS_O_WORKDIR
echo "pbsworkdir"
echo $PBS_O_WORKDIR
EXECDIR=`pwd`
export PATH=./:$PATH
echo $PATH

#Arguments:
## vcf

echo
echo "#######################"
echo $vcf
echo "#######################"

echo
echo
echo Selecting Variants
echo
echo

vcftools \
--vcf $vcf \
--remove-indels \
--remove-filtered-all \
--min-meanDP 2 \
--max-meanDP 10 \
--max-maf 0.05 \
--max-missing 0.25 \
--hwe 0.01 \
--min-alleles 2 \
--max-alleles 2 \
--recode \
--recode-INFO-all \
--out SNP-only.vcf
