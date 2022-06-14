## calculating the standard errors for each species 

sumstats_file = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct.temp",
                           header=T,fill=T)

## need to calculate z-score for each taxon
sumstats_file$zscoresppFST = 0
sumstats_file$zscoresppchromFST = 0

for(spp in unique(sumstats_file$species)) {
  print(spp)
  
  meanfst = mean(sumstats_file$Fst[sumstats_file$species==spp],na.rm=T)
  sdfst = sd(sumstats_file$Fst[sumstats_file$species==spp],na.rm=T)
  
  sumstats_file$zscoresppFST[sumstats_file$species==spp] = (sumstats_file$Fst[sumstats_file$species==spp] - meanfst) / sdfst
  
  png(paste(spp,"zscoretest.png",sep="_"))
  plot(sumstats_file$plotorder[sumstats_file$species==spp],sumstats_file$zscoresppFST[sumstats_file$species==spp])
  abline(h=5,col="red")
  abline(h=-5,col="red")
  dev.off()
  
  for(chrom in unique(sumstats_file$chr[sumstats_file$species==spp])) {
    print(chrom)
    meanfst = mean(sumstats_file$Fst[sumstats_file$species==spp & sumstats_file$chr==chrom],na.rm=T)
    sdfst = sd(sumstats_file$Fst[sumstats_file$species==spp & sumstats_file$chr==chrom],na.rm=T)
    
    sumstats_file$zscoresppchromFST[sumstats_file$species==spp & sumstats_file$chr==chrom] = (sumstats_file$Fst[sumstats_file$species==spp & sumstats_file$chr==chrom] - meanfst) / sdfst
  }
  png(paste(spp,"zscoretest_chroms.png",sep="_"))
  plot(sumstats_file$plotorder[sumstats_file$species==spp],sumstats_file$zscoresppchromFST[sumstats_file$species==spp])
  abline(h=5,col="red")
  abline(h=-5,col="red")
  dev.off()
  
}

pdf("outliers_by_species.pdf")
for(spp in unique(sumstats_file$species)) {
  df = sumstats_file[sumstats_file$species==spp,]
  boxplot(df$zscoresppFST~df$color,main=spp)
  abline(h=5,col="red")
  abline(h=-5,col="red")
}
dev.off()

pdf("barplot_lostruct_fst_outliers.pdf")
for(spp in unique(sumstats_file$species)) {
  df = sumstats_file[sumstats_file$species==spp,]
  
  upper_outliers = aggregate(df$zscoresppFST~df$color,FUN=function(x){
    sum(x >= 5)
  })
  lower_outliers = aggregate(df$zscoresppFST~df$color,FUN=function(x){
    sum(x <= -5)
  })
  not_outliers = aggregate(df$zscoresppFST~df$color,FUN=function(x){
    sum(x > -5 & x < 5)
  })
  
  total = upper_outliers[,2] + lower_outliers[,2] + not_outliers[,2]
  names(total) = upper_outliers[,1]
  
  temp=cbind(lower_outliers[,2],not_outliers[,2],upper_outliers[,2],total)
  rel_temp = cbind(lower_outliers[,2]/total,not_outliers[,2]/total,upper_outliers[,2]/total)
  #barplot(t(rel_temp[,c(1,3)]*100),beside=T,col=c("black","red"),main=spp,ylab="Percent of All Sites With This Color")
  #legend("topright",legend=c("Lower Outliers","Upper Outliers"),fill=c("black","red"))
  #rownames(temp) = upper_outliers[,1]
  #colnames(temp) = c("low","norm","high")
  #temp = temp[,2:4]
  #relative_temp
  #barplot(t(as.matrix(temp)))
  
  print(spp)
  print(sum(total))
  
  
  
}
dev.off()

write.table(sumstats_file,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore.temp")

upper_outliers = aggregate(sumstats_file$zscoresppFST~sumstats_file$species,FUN=function(x){
  sum(x >= 5)
})
lower_outliers = aggregate(sumstats_file$zscoresppFST~sumstats_file$species,FUN=function(x){
  sum(x <= -5)
})

fst  = aggregate(sumstats_file$Fst~sumstats_file$species,FUN=function(x){mean(x,na.rm=T)})
fst2  = aggregate(sumstats_file$Fst~sumstats_file$species,FUN=function(x){median(x,na.rm=T)})

order=c("cri","cur","fus","bel","sin","fla","mel","nit","bru","bil")

rownames(upper_outliers) = upper_outliers[,1]
rownames(lower_outliers) = lower_outliers[,1]
rownames(fst) = fst[,1]

plot(fst[order,2])

plot(fst[order,2],upper_outliers[order,2],ylim=c(0,800),col="red")
points(fst[order,2],lower_outliers[order,2],ylim=c(0,800),pch=0)
abline(lm(upper_outliers[order,2]~fst[order,2]),col="red",lty=3)
abline(lm(lower_outliers[order,2]~fst[order,2]))

## bin outliers by lostruct results

