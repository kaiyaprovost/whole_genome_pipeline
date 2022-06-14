#devtools::install_github("zhengxwen/gdsfmt")
#devtools::install_github("zhengxwen/SNPRelate")

print("Loading packages")
library(gdsfmt)
library(SNPRelate)

print("Reading arguments")
args <- commandArgs(TRUE)
outputstring <- as.character(args[1])
vcfstring <- as.character(args[2])
suffix <- as.character(args[3])

## inport vcf and convert
name = paste(suffix,".gds",sep="")

if (!(file.exists(name))) {
  print("Converting VCF to GDS")
  snpgdsVCF2GDS(vcfstring, name, method="biallelic.only")
}

print("Opening GDS")
genofile = snpgdsOpen(name,readonly = F)
snpgdsSummary(name)

## set up the pops first 
print("Sampling")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
snp_id <- read.gdsn(index.gdsn(genofile, "snp.id"))

snp_sample = sample(snp_id,min(c(length(snp_id),10000)))

snpgdsGDS2PED(genofile,snp.id=snp_sample,ped.fn=name)

showfile.gds(closeall=TRUE)

