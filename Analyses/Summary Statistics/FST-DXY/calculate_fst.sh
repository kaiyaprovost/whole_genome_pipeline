#!/bin/bash
#PBS -l select=1:ncpus=32:mem=128gb
#PBS -l walltime=5000:00:00
#PBS -N bel_bamlist2fst
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
ref=pseudochromosomesSHORT.fasta

## bam_list = list of BAM files
bamlist=Vireo-bellii.bamlist

date
time

time angsd -GL 2 -dosaf 1 -fold 1 -minInd 4 -minMapQ 20 -minQ 20 -nThreads 16 -ref $ref -anc $ref -bam $bamlist -out Vireo-bellii-sfs1

## bam_list = list of BAM files
sonlist=SON_Vireo_bellii.indlist
chilist=CHI_Vireo_bellii.indlist

date >> bel_sfs2.logfile 2>&1
time >> bel_sfs2.logfile 2>&1

## get the global sfs
time realSFS -nSites 50000000 Vireo-bellii-sfs1.saf.idx >> bel_sfs2.logfile 2>&1

## get the by population sfs 
echo >> bel_sfs2.logfile 2>&1
echo "SON SFS"  >> bel_sfs2.logfile 2>&1
/home/kprovost/nas3/angsd/angsd -b $sonlist -anc $ref -out SON_Vireo_bellii-sfs1 -dosaf 1 -gl 1 >> bel_sfs2.logfile 2>&1

echo >> bel_sfs2.logfile 2>&1
echo "CHI SFS" >> bel_sfs2.logfile 2>&1
/home/kprovost/nas3/angsd/angsd -b $chilist -anc $ref -out CHI_Vireo_bellii-sfs1 -dosaf 1 -gl 1 >> bel_sfs2.logfile 2>&1

## now get the 2d sfs for each 
command="time realSFS -nSites 50000000 SON_Vireo_bellii-sfs1.saf.idx CHI_Vireo_bellii-sfs1.saf.idx > SON_CHI_Vireo_bellii.ml"
eval $command >> bel_sfs2.logfile 2>&1

echo "Prep for window analysis"
angsd sites index SON_Vireo_bellii-sfs1.saf.pos.gz

# prepare the fst for easy window analysis etc
time realSFS fst index SON_Vireo_bellii-sfs1.saf.idx CHI_Vireo_bellii-sfs1.saf.idx SON_CHI_Vireo_bellii.ml.1D.txt.1D -fstout SON_CHI_Vireo_bellii_FST 

echo "Get global estimate"

# get the global estimate
time realSFS fst stats SON_CHI_Vireo_bellii_FST.fst.idx 

echo "Sliding Windows" 

time realSFS  fst stats2 SON_CHI_Vireo_bellii_FST.fst.idx -win 100000 -step 10000 > SON_CHI_Vireo_bellii_FST_slidingwindow.fst