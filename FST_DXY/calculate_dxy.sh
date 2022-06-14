#!/bin/bash
#PBS -l select=1:ncpus=16
#PBS -l walltime=5000:00:00
#PBS -N belSON_ANGSD_GL
#PBS -j oe
#PBS -m ae
#PBS -M kprovost@amnh.org
#PBS -k oe

# change to the working directory
cd $PBS_O_WORKDIR
echo "pbsworkdir"
echo $PBS_O_WORKDIR
EXECDIR=`pwd`
export PATH=./:$PATH
echo $PATH

#Arguments:
## ref=reference sequence
ref=/home/kprovost/nas3/genomeresequencingFromLucas/Corv_chrom/pseudochromosomesSHORT.fasta

## if needed tweak parameters between the runs -- double check, make sure it's doing what you want

## bamlist = list of BAM files
echo "run for son"
bamlist=SON_Vireo_bellii.indlist
time angsd -GL 2 -doMaf 11 -doCounts 1 -doMajorMinor 1 -skipTriallelic -ref $ref -sites CHROMS_FROM_DXY_RAWSIZES.txt -nThreads 16 -bam $bamlist -out Vireo-bellii-SONDXY

echo "run for CHI"
bamlist=CHI_Vireo_bellii.indlist
time angsd -GL 2 -doMaf 11 -doCounts 1 -doMajorMinor 1 -skipTriallelic -ref $ref -sites CHROMS_FROM_DXY_RAWSIZES.txt -nThreads 16 -bam $bamlist -out Vireo-bellii-CHIDXY

echo "run for whole species"
bamlist=/Vireo-bellii.bamlist
time angsd -GL 2 -doMaf 11 -doCounts 1 -doMajorMinor 1 -doGlf 2 -minMaf 0.05 -SNP_pval 0.01 -skipTriallelic -minInd 4 -minMapQ 20 -minQ 20 -nThreads 16 -bam $bamlist -out Vireo-bellii-DXY

## resulting mafs files to be run:

sonmaf=Vireo-bellii-SONDXY.mafs.RAW
chimaf=Vireo-bellii-CHIDXY.mafs.RAW
sonlen=`cat $sonmaf | wc -l`
#chilen=`cat $chimaf | wc -l` ## you don't need this if you aren't going to run them separately

module load R-3.4.1 ## or a different version
gunzip -f $sonmaf.gz
gunzip -f $chimaf.gz

echo "son maf file is: ${sonmaf}"
echo "length: ${sonlen}"
echo "chi maf file is: ${chimaf}"
echo "length: ${chilen}"

Rscript calcDxy.R --popA $sonmaf --popB $chimaf --totLen $sonlen > log_bel.txt 2>&1

mv Dxy_persite.txt bel_SON_Dxy_persite.txt

## this should give the same result no matter whether you use the $sonlen or $chilen
## so you don't have to run this below if you don't want to 
# Rscript calcDxy.R --popA $chimaf --popB $sonmaf --totLen $chilen
# mv Dxy_persite.txt bel_CHI_Dxy_persite.txt
