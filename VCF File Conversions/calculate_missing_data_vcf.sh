#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=5000:00:00
#PBS -N generate_imiss_from_vcf
#PBS -j oe
#PBS -m ae
#PBS -M kprovost@amnh.org
#PBS -k oe

source activate py36

cd /vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/ 

for vcf in *vcf; do 
echo $vcf; 
#python3 /home/kprovost/nas3/ANGSD_pipeline/vcf_remove_2tabs.py $vcf; 

if [ ! -f ${vcf%}*imiss ]
then
/home/kprovost/opt/bin/vcftools --vcf $vcf --missing-indv --out $vcf.imiss
#else 
#echo "Already converted individual missing"
fi

if [ ! -f ${vcf%}*lmiss ]
then
/home/kprovost/opt/bin/vcftools --vcf $vcf --missing-site --out $vcf.lmiss
#else 
#echo "Already converted site missing"
fi

# if [ ! -f ${vcf%}*singletons ]
# then
# /home/kprovost/opt/bin/vcftools --vcf $vcf --singletons --out $vcf.singletons
# else 
# echo "Already converted singletons"
# fi

done;