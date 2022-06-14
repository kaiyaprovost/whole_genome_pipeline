# coding: utf-8

cd /home/kprovost/nas4/DEMOGRAPHY

ls *vcf

for i in *vcf; do
echo "###########################################"
echo $i; 
j=${i%.geno.vcf}; 
echo $j; 


head -n 2 $j.geno.vcf



## bgzip performs a blockwise compression
## The -c flag directs bgzip to leave the original vcf file 
##   untouched and create a new file for the vcf.gz
if [ -f $j.geno.vcf.gz ]; then
echo "File $j.geno.vcf.gz exists."
else
echo "No GZ file"
#bgzip -c $j.geno.vcf > $j.geno.vcf.gz
fi



## tabix indexes the file for searching
if [ -f $j.geno.vcf.gz.tbi ]; then
echo "File $j.geno.vcf.gz.tbi exists."
else
echo "No TBI file"
#tabix $j.geno.vcf.gz
fi



if [ -f $j.geno.bed ]; then
echo "File $j.geno.bed exists."
head -n 1 $j.geno.bed
else
echo "No BED file"
#vcf2bed --do-not-split --do-not-sort < $j.geno.vcf > $j.geno.bed
## Print the first 1 lines of this file
#head -n 1 $j.geno.bed
fi

if [ -f $j.allelecounts.txt ]; then
echo "File $j.allelecounts.txt exists."
else
echo "No COUNTS file"
#python -m momi.read_vcf --no_aa --verbose $j.geno.vcf.gz $j.popassignments.csv $j.allelecounts.txt --bed $j.geno.bed; 
fi

if [ -f $j.geno.sfs.txt ]; then
echo "File $j.geno.sfs.txt exists."
else
echo "No SFS file"
#python -m momi.extract_sfs $j.geno.sfs.txt 20 $j.geno.vcf
fi

echo ""
echo ""

done