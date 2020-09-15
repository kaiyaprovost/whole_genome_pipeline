cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/VCFS/"


for i in */*window.vcf; do
echo $i;
temp=$i.temp;
outsamp1=`head -2 "${i}" |tail -1 |tr '\t' '\n' |wc -l`;
fixcols=9;
outsamp=`echo "$(($outsamp1-$fixcols))"`;
echo $outsamp;
perl "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/vcf2MS.pl" \
"${i}" \
"${temp}" \
"${outsamp}";
gzip -vf $i;
mv -v */*temp ../TEMPS/
done;
