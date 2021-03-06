#!/bin/bash
#Arguments:
## ref=reference sequence

date
time

module load bwa-0.7.15
module load fastqc-0.11.5
module load R-3.4.1

for i in */; do
	cd $i
	name="${i///}"
	read1=RAPiD-Genomics*_R1_*
	read2=RAPiD-Genomics*_R2_*
	header=`zcat $read1 | head -1`
	IFS=':' read -a header <<< "$header"
	INSTRUMENT=${header[0]}
	RUN_ID=${header[1]}
	FLOWCELL_BARCODE=${header[2]}
	LANE=${header[3]}
	ID=$FLOWCELL_BARCODE.$LANE
	PU=$FLOWCELL_BARCODE.$LANE.$name
	SM=$name
	PL=ILLUMINA
	LB=TrueSeq
        echo $name $read1 $read2 $ID $PU $SM $PL $LB
	
	echo
	echo "#######################"
	echo $name
	echo "#######################"

	echo
	echo
	echo Trimming adapters
	echo
	echo

	time java -jar /home/kprovost/nas1/genomeresequencingFromLucas/Trimmomatic-0.36/trimmomatic-0.36.jar \
	PE $read1 $read2 \
	$name.forward_paired.fq.gz \
	$name.forward_unpaired.fq.gz \
	$name.reverse_paired.fq.gz \
	$name.reverse_unpaired.fq.gz \
	ILLUMINACLIP:/home/kprovost/nas1/genomeresequencingFromLucas/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10:8:true

	echo
	echo	
	echo Quality assessment
	echo
	echo

	time fastqc $name.forward_paired.fq.gz
	time fastqc $name.reverse_paired.fq.gz

	echo
	echo	
	echo Running BWA alignment
	echo
	echo
	
	time bwa mem \
	-M \
	-t 28 \
	$ref \
	$name.forward_paired.fq.gz \
	$name.reverse_paired.fq.gz \
	> $name.sam

	echo	
	echo
	echo Converting SAM to sorted BAM
	echo
	echo	
	
	mkdir tmp
	time java -jar /home/kprovost/nas1/genomeresequencingFromLucas/picard.jar SortSam \
	INPUT=$name.sam \
	OUTPUT=$name.bam \
	SORT_ORDER=coordinate \
	TMP_DIR=`pwd`/tmp
	rm -r tmp/

	echo	
	echo
	echo Collecting alignment metrics
	echo
	echo
	
	time java -jar /home/kprovost/nas1/genomeresequencingFromLucas/picard.jar CollectAlignmentSummaryMetrics \
	INPUT=$name.bam \
	OUTPUT=$name.alignment_metrics.txt \
	R=$ref \
	METRIC_ACCUMULATION_LEVEL=SAMPLE \
	METRIC_ACCUMULATION_LEVEL=READ_GROUP \

	echo
	echo
	echo Collect Insert Size Metrics
	echo
	echo
	
	time java -jar /home/kprovost/nas1/genomeresequencingFromLucas/picard.jar CollectInsertSizeMetrics \
	INPUT=$name.bam \
	OUTPUT=$name.insert_metrics.txt \
	HISTOGRAM_FILE=$name.insert_size_histogram.pdf
	
	echo	
	echo
	echo Collect Coverage matrics
	echo
	echo
	
	time java -jar /home/kprovost/nas1/genomeresequencingFromLucas/picard.jar CollectRawWgsMetrics \
	I=$name.bam \
	O=$name.raw_wgs_metrics.txt \
	R=$ref \
	INCLUDE_BQ_HISTOGRAM=true

	echo
	echo
	echo Adding Read Group
	echo
	echo
	
	time java -jar /home/kprovost/nas1/genomeresequencingFromLucas/picard.jar AddOrReplaceReadGroups \
	I=$name.bam \
	O=$name.groups_added.bam \
	RGID=$ID \
	RGLB=$LB \
	RGPL=$PL \
	RGPU=$PU \
	RGSM=$SM

	echo
	echo
	echo Deduplication
	echo
	echo
	
	time java -jar /home/kprovost/nas1/genomeresequencingFromLucas/picard.jar MarkDuplicates \
	TMP_DIR=tmp \
	I=$name.groups_added.bam \
	O=$name.dedup.bam \
	METRICS_FILE=$name.dedup.metrics.txt \
	REMOVE_DUPLICATES=false \
	TAGGING_POLICY=All

	echo
	echo
	echo Indexing BAM
	echo
	echo
	
	time java -jar /home/kprovost/nas1/genomeresequencingFromLucas/picard.jar BuildBamIndex \
	I=$name.dedup.bam

	echo
	echo
	echo Validating BAM
	echo
	echo

	time java -jar /home/kprovost/nas1/genomeresequencingFromLucas/picard.jar ValidateSamFile \
	I=$name.dedup.bam \
	O=$name.validate.txt \
	MODE=SUMMARY	
	
	rm $name.sam

	cd ../

	echo
	echo "#######################"
	echo $name DONE!
	echo "#######################"
	echo
	echo

done

date
time
