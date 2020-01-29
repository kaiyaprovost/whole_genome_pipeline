path = "/Users/kprovost/Documents/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/genomeresequencingFromLucas/for_AMN_245109/STATS_FILES/"
setwd(path)
file = "/Users/kprovost/Documents/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/genomeresequencingFromLucas/for_AMN_245109/STATS_FILES/vcfstats_by_sample.csv"
csv = read.csv(file)

bel = droplevels(csv[csv$SPECIES=="BILINEATA",])
bel = bel[bel$SPECIMEN!="ALL",]
summary(bel)

pairs(bel[,c(30,11:15)]) ## start 3 end 28, 13: interesting

pc1 = bel$PC1
meanpc1 = mean(pc1)
stdevpc1 = sd(pc1)
zscore = abs((pc1 - meanpc1) / stdevpc1)

flagged = which(zscore>=3)

library(corrplot)

for (i in 1:length(unique(csv$SPECIES))) {
  species = (unique(csv$SPECIES)[i])
  
  if (species != "CRI") {
  
  print(species)
  df = droplevels(csv[csv$SPECIES==species,])
  df = df[df$SPECIMEN!="ALL",]
  endcol = 29+nrow(df) ## 30 is pc1 -- also remove 23:29,16, count
  corrtable = cor(df[,c(3:8,10:15,17:23,30:endcol-1)])
  
  corrtable_small = corrtable[c(21:(19+nrow(df))),c(1:20)]
  
  png(paste("CorrelationPlot_",species,".png",sep=""))
  corrplot(corrtable_small,diag=T,method="ellipse")
  dev.off()
  }
  
}

#####
## below is old

#pdf(paste("/Users/kprovost/Documents/Dissertation/CHAPTER2_GENOMES/ANALYSIS/eigenvalues/AllSpecies_individuals_eigen.pdf",sep=""))

for (name in c("bel","bil","bru","fla","fus","cur","cri","mel","sin")) {
print(name)
beldatstring = paste("/Users/kprovost/Documents/Dissertation/CHAPTER2_GENOMES/ANALYSIS/eigenvalues/unfiltered_analysis",
                     name,"BOTHnomiss_fullPcatable.txt",sep="")
belpopstring = paste("/Users/kprovost/Documents/Dissertation/CHAPTER2_GENOMES/ANALYSIS/populations/populations_for_vcf_",
                     name,"BOTH.txt",sep="")

beldat = read.csv(beldatstring,sep=" ")
belpop = read.csv(belpopstring,sep="\t")

bel = merge(belpop,beldat,by="sample.id")

#boxplot(bel[,3:ncol(bel)])

#pairs(bel[,3:ncol(bel)],col=bel$pop_code)

side1 = ceiling(sqrt(ncol(bel)-2))
side2 = ceiling((ncol(bel)-2)/side1)

#png(paste("/Users/kprovost/Documents/Dissertation/CHAPTER2_GENOMES/ANALYSIS/eigenvalues/",name,"_individuals_eigen.png",sep=""))
par(mfrow=c(side1,side2),mar=c(2,2,1,1))
for (i in 3:ncol(bel)){
  plot(bel[,i],col=bel$pop_code,main=paste(names(bel)[i],name,sep=" ",collapse=" "),
       ylim=c(-1,1),pch=as.numeric(bel$sample.id))
}

plot(bel[,i],type="n",xaxt="n",yaxt="n")
legend("center", legend=bel$pop_code, pch=as.numeric(bel$sample.id), col=bel$pop_code,
       bg="n",bty="n",inset=c(-0.1,0)) ## before was -0.75

dev.off()
}
dev.off()
