orig = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION/bel/GENOME/notused/CHI_bel.recode.fixedchroms.PREDICT.txt",
                  header=T,sep="\t")
wind = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION/CHI_bel.recode.fixedchroms.PREDICT.txt_w100000_o10000_chrom1.txt",
                  header=T,sep=",")

bs_orig = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION/bel/NOTGENOME/CHI_bel.recode.fixedchroms.PREDICT.BSCORRECTED.txt",
                     header=T,sep="\t")
bs_wind = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION/CHI_bel.recode.fixedchroms.PREDICT.BSCORRECTED.txt_w100000_o10000_chrom1.txt",
                     header=T,sep=",")

bs_orig$size=(bs_orig$end-bs_orig$start)

sum(is.na(bs_orig$recombRate))

png("bootstraps.png")
plot(orig$recombRate[orig$chrom==1],ylim=c(0,12e-10))
points(bs_orig$recombRate[bs_orig$chrom==1],col="grey")
points(bs_orig$CI95LO[bs_orig$chrom==1],col="cyan")
points(bs_orig$CI95LO[bs_orig$chrom==1],col="red")
dev.off()

png("bootstrap_log.png")
plot(log10(bs_orig$recombRate))
dev.off()

dim(orig)
dim(bs_orig)

png("original_vs_bootstrap_recomb_bel.png")
plot(orig$recombRate,bs_orig$recombRate)
points(orig$recombRate,bs_orig$CI95LO,cex=0.1,col="grey")
points(orig$recombRate,bs_orig$CI95HI,cex=0.1,col="grey")
abline(a=0,b=1,col="red")
dev.off()

dim(wind)
dim(bs_wind)

png("original_vs_bootstrap_recombwindows_bel.png")
plot(wind$weighted_recomb,bs_wind$weighted_recomb)
abline(a=0,b=1,col="red")
dev.off()


orig$midPos = (orig$start+orig$end)/2
wind$midPos = (wind$windowstarts+wind$windowstops)/2

orig2 = orig[orig$chrom==1,]
orig2$midPos = orig2$midPos-2500

y=merge(orig2[,c("midPos","recombRate")],wind[,c("midPos","weighted_recomb")],all=F)

png("original_vs_windows_recomb_bel_chr1.png")
plot(y$recombRate,y$weighted_recomb)
abline(a=0,b=1,col="red")
dev.off()


bs_orig$midPos = (bs_orig$start+bs_orig$end)/2
bs_wind$midPos = (bs_wind$windowstarts+bs_wind$windowstops)/2

bs_orig2 = bs_orig[bs_orig$chrom==1,]
bs_orig2$midPos = bs_orig2$midPos-2500

bs_y=merge(bs_orig2[,c("midPos","recombRate")],bs_wind[,c("midPos","weighted_recomb")],all=F)

png("bootstrap_vs_windows_recomb_bel_chr1.png")
plot(bs_y$recombRate,bs_y$weighted_recomb)
abline(a=0,b=1,col="red")
dev.off()



## looking at two different versions of recomb data
a = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION/bel/GENOME/Vireo-bellii-called.geno.PseudoNC.all.sorted.PREDICT.txt_w100000_o10000_genome.txt",
               header=T,sep=",")
b = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION/bel/GENOME/Vireo-bellii-called.geno.PseudoNC.all.sorted.sorted.PREDICT.txt_w100000_o10000_genome.txt",
               header=T,sep=",")

chi = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION/bel/GENOME/CHI_bel.recode.fixedchroms.PREDICT.txt_w100000_o10000_genome.txt",
                 header=T,sep=",")
son = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION/bel/GENOME/SON_bel.recode.sorted.PREDICT.txt_w100000_o10000_genome.txt",
                 header=T,sep=",")

sum(a$windowstarts!=son$windowstarts)

png("dif_bel_versions_recomb.png",height=4,width=6,units="in",res=300)
par(mfrow=c(2,3))
plot(a$weighted_recomb,b$weighted_recomb)
abline(a=0,b=1,col="red")
plot(a$weighted_recomb,chi$weighted_recomb)
abline(a=0,b=1,col="red")
plot(a$weighted_recomb,son$weighted_recomb)
abline(a=0,b=1,col="red")
plot(b$weighted_recomb,chi$weighted_recomb)
abline(a=0,b=1,col="red")
plot(b$weighted_recomb,son$weighted_recomb)
abline(a=0,b=1,col="red")
plot(chi$weighted_recomb,son$weighted_recomb)
abline(a=0,b=1,col="red")
dev.off()

df = cbind(a=a$weighted_recomb,b=b$weighted_recomb,chi=chi$weighted_recomb,son=son$weighted_recomb)
corrplot::corrplot(cor(df,use="pairwise.complete.obs"),method="number")
