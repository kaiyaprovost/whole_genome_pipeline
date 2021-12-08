#!/bin/bash

cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno"

for i in *geno/; do
echo $i
cd $i 
j=${i%/}
echo "xxx $j"
cat *vcf > ../$j.PseudoNC.all.vcf.temp
cd ../
done

for k in *PseudoNC.all.vcf.temp; do
echo $k
p=${k%.temp}
echo "yyy $p"
uniq $k > $p
done 