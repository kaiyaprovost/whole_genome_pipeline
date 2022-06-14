cd /vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/

#for i in *unfilter*vcf; do
#echo $i;
#qsub -v vcf=$i ~/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP05.2-basicBCFstats.job ;
#done; 

# for i in *unfilter*vcf; do
# echo $i;
# qsub -v vcf=$i ~/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP05.2-basicVCFstats.job ;
# done; 

for i in *unfilter*vcf; do
echo $i;
qsub -v vcf=$i ~/nas3/genomeresequencingFromLucas/for_AMN_245109/STEP05.2-singletons.job ;
done; 

#source activate py36
#for i in *bcfstats; do
#/home/kprovost/nas3/bcftools/misc/plot-vcfstats $i -p plots/${i}_plot/
#done; 

cd ~ 

