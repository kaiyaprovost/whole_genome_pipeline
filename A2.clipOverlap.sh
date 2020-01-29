#!/bin/bash
#PBS -l select=1:ncpus=2
#PBS -l walltime=99999:00:00
#PBS -N clipOverlap
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
## bam = BAM file

# for i in /vz-nas1-active/ProcessedGenomicReads/FILTERED/fu*/REHEADERED/*F073_AMN_245111_P002_WF06**reheadered.bam; do qsub -v bam=$i ~/nas3/ANGSD_pipeline/A2.clipOverlap.sh; done;
# for i in /vz-nas1-active/ProcessedGenomicReads/FILTERED/fu*/REHEADERED/NOTMERGED/F073_AMN_245111_P002_WF06**reheadered.bam; do qsub -v bam=$i ~/nas3/ANGSD_pipeline/A2.clipOverlap.sh; done;

#name=`echo $bam | cut -d '.' -f1`
name=${bam%.dedup.reheadered.bam}
name=`echo $name | rev | cut -d'/' -f 1 | rev`

echo
echo "#######################"
echo $name
echo "#######################"

echo
echo
echo Indel Realigner
echo
echo

cd /vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/CLIPPED/

# time java -Djava.io.tmpdir=/home/kprovost/nas3/tmp/ -jar /home/kprovost/nas3/genomeresequencingFromLucas/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
# -T RealignerTargetCreator \
#~/nas3/samtools/samtools view -u $bam > $bam.redo

/home/kprovost/nas3/bamUtil/bin/bam clipOverlap --in $bam --noeof --out $name.clipped.bam --stats --params --unmapped --storeOrig

