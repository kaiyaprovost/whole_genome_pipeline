#!/bin/bash
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -l walltime=5000:00:00
#PBS -N lostruct_vcf_subset
#PBS -j oe
#PBS -m ae
#PBS -M kprovost@amnh.org
#PBS -k oe

#cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/"
#vcffile=test.vcf.gz
#smallfile=./ANALYSIS/Campylorhynchus.brunneicapillus.called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.nospace.coords_SMALL.csv
#python3 "/Users/kprovost/Documents/Github/whole_genome_pipeline/VCF File Conversions/vcf_extract_lostruct_windows.py" $vcffile $smallfile

cd /vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/

vcffile=/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/FOR_LOSTRUCT/Campylorhynchus-brunneicapillus-called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.nospace.vcf.gz
smallfile=/home/kprovost/nas3/LOSTRUCT/Campylorhynchus.brunneicapillus.called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.nospace.coords_SMALL.csv
python3 /home/kprovost/nas3/LOSTRUCT/vcf_extract_lostruct_windows.py $vcffile $smallfile


vcffile=/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/FOR_LOSTRUCT/Phainopepla-nitens-called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.nospace.vcf.gz
smallfile=/home/kprovost/nas3/LOSTRUCT/Phainopepla.nitens.called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.nospace.coords_SMALL.csv
python3 /home/kprovost/nas3/LOSTRUCT/vcf_extract_lostruct_windows.py $vcffile $smallfile


vcffile=/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/FOR_LOSTRUCT/Polioptila-melanura-called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.nospace.vcf.gz
smallfile=/home/kprovost/nas3/LOSTRUCT/Polioptila.melanura.called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.nospace.coords_SMALL.csv
python3 /home/kprovost/nas3/LOSTRUCT/vcf_extract_lostruct_windows.py $vcffile $smallfile


vcffile=/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/FOR_LOSTRUCT/Toxostoma-crissale-called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.nospace.vcf.gz
smallfile=/home/kprovost/nas3/LOSTRUCT/Toxostoma.crissale.called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.nospace.coords_SMALL.csv
python3 /home/kprovost/nas3/LOSTRUCT/vcf_extract_lostruct_windows.py $vcffile $smallfile


vcffile=/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/FOR_LOSTRUCT/Toxostoma-curvirostre-called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.nospace.vcf.gz
smallfile=/home/kprovost/nas3/LOSTRUCT/Toxostoma.curvirostre.called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.nospace.coords_SMALL.csv
python3 /home/kprovost/nas3/LOSTRUCT/vcf_extract_lostruct_windows.py $vcffile $smallfile


vcffile=/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/FOR_LOSTRUCT/Cardinalis-sinuatus-called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.copy.nospace.vcf.gz
smallfile=/home/kprovost/nas3/LOSTRUCT/Cardinalis.sinuatus.called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.nospace.coords_SMALL.csv
python3 /home/kprovost/nas3/LOSTRUCT/vcf_extract_lostruct_windows.py $vcffile $smallfile


