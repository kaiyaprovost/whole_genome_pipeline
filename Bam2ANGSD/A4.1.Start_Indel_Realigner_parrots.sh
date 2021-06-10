#!/bin/bash
#PBS -l select=1:ncpus=4:mem=32gb
#PBS -l walltime=99999:00:00
#PBS -N A4.1_startindel
#PBS -j oe
#PBS -m ae
#PBS -M kprovost@amnh.org
#PBS -k oe

## todo: sin for both
spp="Parrots"
sp3=`echo $spp | cut -c1-3`

for bam in /vz-nas1-active/RawGenomicReads/FASTQ_SCREEN_RESULTS/Parrots/CLIPPED/*clipped.bam; do
#bam=${prebam%.dedup.reheadered.bam}
#echo $bam
name=`echo $bam | rev | cut -d'/' -f 1 | rev`
echo $name
#inputbam=/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/${name}.dedup.reheadered.clipped.bam
inputbam=$bam
ls $inputbam
inputbai=${inputbam%.bam}.bai
if [ ! -f $inputbai ]; then
time java -Djava.io.tmpdir=/home/kprovost/nas3/tmp/ -jar /home/kprovost/nas3/genomeresequencingFromLucas/picard.jar BuildBamIndex I=$inputbam O=$inputbai
fi
qsub -v bam=$inputbam /home/kprovost/nas3/ANGSD_pipeline/A4.Indel_Realigner_Parrots.job
echo "#####"
done


#for prebam in /vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/CLIPPED/${spp}/AMN*.clipped.bam; do 
#echo $prebam
for bam in /vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/CLIPPED/${spp}/NOTMERGED/*clipped.bam; do
#bam=${prebam%.dedup.reheadered.bam}
#echo $bam
name=`echo $bam | rev | cut -d'/' -f 1 | rev`
echo $name
#inputbam=/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/${name}.dedup.reheadered.clipped.bam
inputbam=$bam
ls $inputbam
inputbai=${inputbam%.bam}.bai
if [ ! -f $inputbai ]; then
time java -Djava.io.tmpdir=/home/kprovost/nas3/tmp/ -jar /home/kprovost/nas3/genomeresequencingFromLucas/picard.jar BuildBamIndex I=$inputbam O=$inputbai
fi
qsub -v bam=$inputbam /home/kprovost/nas3/ANGSD_pipeline/A4.Indel_Realigner_${sp3}*.job
echo "#####"
done