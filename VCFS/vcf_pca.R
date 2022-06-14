#devtools::install_github("zhengxwen/gdsfmt")
#devtools::install_github("zhengxwen/SNPRelate")

print("Loading packages")
library(gdsfmt)
library(SNPRelate)

print("Reading arguments")
args <- commandArgs(TRUE)
outputstring <- as.character(args[1])
#outputstring = "/Users/kprovost/Documents/Dissertation/CHAPTER1_REVIEW/SLIM/first_runs/"
vcfstring <- as.character(args[2])
#vcfstring = "/Users/kprovost/Documents/Dissertation/CHAPTER1_REVIEW/SLIM/first_runs/slim-model3_isolation-1550691782/VCFS/model3_isolation-1550691782-1-overlaid.vcf"
popstring <- as.character(args[3])
#popstring = "/Users/kprovost/Documents/Dissertation/CHAPTER1_REVIEW/SLIM/first_runs/popcodes_model3_isolation-1550691782.txt"
suffix <- as.character(args[4])
#suffix="testvcf"

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
print("Setting up Populations")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pops = as.matrix(read.csv(popstring,sep="\t"))
snp_id <- read.gdsn(index.gdsn(genofile, "snp.id"))

if (length(sample.id) != nrow(pops)) {
  print(length(sample.id))
  print(nrow(pops))
  print(sample.id)
  print(pops)
  stop("Population and samples files don't have same number of individuals!")
}

test = cbind(sample.id,pops)
if (sum(test[,1] != test[,2]) != 0) {
  stop("Samples are not in the same order!")
} else {
  print("Population file read, names match")
}

pops_only = pops[,2]

x = index.gdsn(genofile,"pops",silent=T)
if (is.null(x)) { 
  
  pops_only = pops[,2]
  add.gdsn(genofile,name="pops",val=pops_only)
  
  }

pop_code <- read.gdsn(index.gdsn(genofile, "pops"))

## prune the snps before PCA? optional but recommended

print("Calculating PCA")
pca <- snpgdsPCA(genofile, autosome.only = F,
                 remove.monosnp=F) ## failing here for color issues 
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

fulltable = data.frame(sample.id = pca$sample.id,
                       stringsAsFactors = FALSE)

fulltable = cbind (fulltable, pca$eigenvect)

write.table(fulltable,
            paste(outputstring,suffix,"_fullPcatable.txt",sep=""))

## pca without pop information
print("PCA Without Populations")
png(paste(outputstring,suffix,"_no_pop_pca.png",sep=""))
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
dev.off()

## pca with pop info
print("PCA With Populations")
png(paste(outputstring,suffix,"_with_pop_pca.png",sep=""))
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

par(xpd=T,mar=c(5,4,4,8))
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab$pop), pch=1, col=1:nlevels(tab$pop),
       bg="n",bty="n",inset=c(-0.1,0)) ## before was -0.75
dev.off()


## first four pcs

print("First Four PCs")
png(paste(outputstring,suffix,"_with_pop_firstfourpca.png",sep=""))
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab$pop, labels=lbls)
dev.off()


## remove problematic PC1 individuals? 

pc1 = fulltable[,3]
meanpc1 = mean(pc1)
stdevpc1 = sd(pc1)
zscore = abs((pc1 - meanpc1) / stdevpc1)

flagged = which(zscore>=3)
flagged_names = pca$sample.id[flagged]

print(flagged_names)
good_names = pca$sample.id[-(flagged)]

pca2 <- snpgdsPCA(genofile, autosome.only = F,
                 remove.monosnp=F, 
                 sample.id = good_names) ## failing here for color issues 
pc.percent2 <- pca2$varprop*100
head(round(pc.percent2, 2))

fulltable2 = data.frame(sample.id = pca2$sample.id,
                       stringsAsFactors = FALSE)

fulltable2 = cbind(fulltable2, pca2$eigenvect)

write.table(fulltable2,
            paste(outputstring,suffix,"_fullPcatable_redo.txt",sep=""))

print("PCA With Populations2")
png(paste(outputstring,suffix,"_with_pop_pca2.png",sep=""))
tab2 <- data.frame(sample.id = pca2$sample.id,
                  pop = factor(pop_code)[match(pca2$sample.id, sample.id)],
                  EV1 = pca2$eigenvect[,1],    # the first eigenvector
                  EV2 = pca2$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab2)
par(xpd=T,mar=c(5,4,4,8))
plot(tab2$EV2, tab2$EV1, col=as.integer(tab2$pop), xlab="eigenvector 2", ylab="eigenvector 1")

## parallel coordinates plot 
# print("Parallel coordinates")
# library(MASS)
# 
# num_eigenvect = pca$eigenval[!(is.na(pca$eigenval))]
# 
# png(paste(outputstring,suffix,"_parallel_coordinates.png",sep=""))
# datpop <- factor(pop_code)[match(pca$sample.id, sample.id)]
# parcoord(pca$eigenvect[,1:min(16,num_eigenvect)], col=datpop)
# dev.off()
# 
# ## snp correlations bt eigen and genotype -- not sure what this does 
# png(paste(outputstring,suffix,"_snp_correlations.png",sep=""))
# chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
# CORR <- snpgdsPCACorr(pca, genofile, eig.which=1:4)
# 
# savepar <- par(mfrow=c(3,1), mai=c(0.3, 0.55, 0.1, 0.25))
# for (i in 1:3)
# {
#   plot(abs(CORR$snpcorr[i,]), ylim=c(0,1), xlab="", ylab=paste("PC", i),
#        col=chr, pch="+")
# }
# dev.off()

showfile.gds(closeall=TRUE)

