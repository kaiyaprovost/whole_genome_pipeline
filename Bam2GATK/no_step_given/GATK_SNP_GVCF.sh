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

echo
echo
echo Combine GVCFs
echo
echo

time java -jar /home/lmoreira/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R $ref \
--variant JGG1800.recal.raw.g.vcf \
--variant JGG1801.recal.raw.g.vcf \
--variant JGG1989.recal.raw.g.vcf \
--variant JGG2060.recal.raw.g.vcf \
--variant JGG2061.recal.raw.g.vcf \
--variant JGG2244.recal.raw.g.vcf \
--variant LRM75.recal.raw.g.vcf \
--variant LRM76.recal.raw.g.vcf \
--variant LRM77.recal.raw.g.vcf \
--variant LRM78.recal.raw.g.vcf \
--variant LRM79.recal.raw.g.vcf \
--variant LRM80.recal.raw.g.vcf \
--variant LRM81.recal.raw.g.vcf \
--variant LRM82.recal.raw.g.vcf \
--variant LRM83.recal.raw.g.vcf \
--variant LRM84.recal.raw.g.vcf \
--variant PAC2027.recal.raw.g.vcf \
--variant PAC2923.recal.raw.g.vcf \
--variant PAC3038.recal.raw.g.vcf \
--variant PRS3469.recal.raw.g.vcf \
-o allsamples.recal.raw.vcf
-nt 64
