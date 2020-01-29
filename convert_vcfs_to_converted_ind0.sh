chrom=""
for i in /Users/kprovost/Dropbox\ \(AMNH\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/NOT_CONVERTED_VCFS/BILINEATA/*${chrom}*.vcf.fixedchroms; do echo $i;
#gzip -vf "${i%.fixedchroms}" 
python3 "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/angsd_ind0_to_cat.py" "${i}" "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/Amphispiza-bilineata.bamlist"
#gzip -f "${i}" 
done
for i in /Users/kprovost/Dropbox\ \(AMNH\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/NOT_CONVERTED_VCFS/FLAVICEPS/*${chrom}*.vcf.fixedchroms; do echo $i;
#gzip -vf "${i%.fixedchroms}" 
python3 "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/angsd_ind0_to_cat.py" "${i}" "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/Auriparus-flaviceps.bamlist"
#gzip -f "${i}" 
done
for i in /Users/kprovost/Dropbox\ \(AMNH\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/NOT_CONVERTED_VCFS/BRUNNEICAPILLUS/*${chrom}*.vcf.fixedchroms; do echo $i;
#gzip -vf "${i%.fixedchroms}" 
python3 "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/angsd_ind0_to_cat.py" "${i}" "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/Campylorhynchus-brunneicapillus.bamlist"
#gzip -f "${i}" 
done
for i in /Users/kprovost/Dropbox\ \(AMNH\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/NOT_CONVERTED_VCFS/SINUTAUS/*${chrom}*.vcf.fixedchroms; do echo $i;
#gzip -vf "${i%.fixedchroms}" 
python3 "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/angsd_ind0_to_cat.py" "${i}" "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/Cardinalis-sinuatus.bamlist"
#gzip -f "${i}" 
done
for i in /Users/kprovost/Dropbox\ \(AMNH\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/NOT_CONVERTED_VCFS/FUSCA/*${chrom}*.vcf.fixedchroms; do echo $i;
#gzip -vf "${i%.fixedchroms}" 
python3 "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/angsd_ind0_to_cat.py" "${i}" "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/Melozone-fusca.bamlist"
#gzip -f "${i}" 
done
for i in /Users/kprovost/Dropbox\ \(AMNH\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/NOT_CONVERTED_VCFS/NITENS/*${chrom}*.vcf.fixedchroms; do echo $i;
#gzip -vf "${i%.fixedchroms}" 
python3 "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/angsd_ind0_to_cat.py" "${i}" "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/Phainopepla-nitens.bamlist"
#gzip -f "${i}" 
done
for i in /Users/kprovost/Dropbox\ \(AMNH\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/NOT_CONVERTED_VCFS/MELANURA/*${chrom}*.vcf.fixedchroms; do echo $i;
#gzip -vf "${i%.fixedchroms}" 
python3 "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/angsd_ind0_to_cat.py" "${i}" "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/Polioptila-melanura.bamlist"
#gzip -f "${i}" 
done
for i in /Users/kprovost/Dropbox\ \(AMNH\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/NOT_CONVERTED_VCFS/CRISSALE/*${chrom}*.vcf.fixedchroms; do echo $i;
#gzip -vf "${i%.fixedchroms}" 
python3 "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/angsd_ind0_to_cat.py" "${i}" "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/Toxostoma-crissale.bamlist"
#gzip -f "${i}" 
done
for i in /Users/kprovost/Dropbox\ \(AMNH\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/NOT_CONVERTED_VCFS/CURVIROSTRE/*${chrom}*.vcf.fixedchroms; do echo $i;
#gzip -vf "${i%.fixedchroms}" 
python3 "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/angsd_ind0_to_cat.py" "${i}" "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/Toxostoma-curvirostre.bamlist"
#gzip -f "${i}" 
done
for i in /Users/kprovost/Dropbox\ \(AMNH\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/NOT_CONVERTED_VCFS/BELLII/*${chrom}*.vcf.fixedchroms; do echo $i;
#gzip -vf "${i%.fixedchroms}" 
python3 "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/angsd_ind0_to_cat.py" "${i}" "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/Vireo-bellii.bamlist"
#gzip -f "${i}" 
done