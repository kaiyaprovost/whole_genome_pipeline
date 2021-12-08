## take windows and turn them to data that can be used in the big file
file ="/Users/kprovost/Dropbox (AMNH)/BELLII_MERGED.temp"
dffull = read.table(file,header=T)

stats=c("div_per_site","Fu.Li.D","Fu.Li.D.1","Fu.Li.F","Fu.Li.F.1","hap.diversity.within",
        "kurt_haplotype.counts","kurt_minor.allele.freqs","mean_haplotype.counts","mean_minor.allele.freqs",
        "mean_nuc.diversity.within","n.biallelic.sites","n.segregating.sites","n.segregating.sites.1","n.sites","Pi","skew_haplotype.counts",
        "skew_minor.allele.freqs","skew_nuc.diversity.within","Tajima.D","Tajima.D.1","theta_Achaz.Tajima","theta_Achaz.Tajima.1",
        "theta_Achaz.Watterson","theta_Achaz.Watterson.1","theta_Fu.Li","theta_Fu.Li.1","theta_Tajima","theta_Watterson","theta_Watterson.1",
        "var_haplotype.counts","var_minor.allele.freqs","var_nuc.diversity.within")

#df = dffull[,c("file",stats)]
df = dffull

window = sapply(df$file,FUN=function(x){
  split1 = strsplit(x,"_")[[1]]
  if(length(split1)>=7) {
    return(as.numeric(strsplit(split1[7],"\\.")[[1]][1]))
  } else {
    return(-1)
  }
})
df$window = window

chrom = sapply(df$file,FUN=function(x){
  split1 = strsplit(x,"_")[[1]]
  return(as.character(strsplit(split1[4],"\\.")[[1]][1]))
})
df$chrom = chrom

species = sapply(df$file,FUN=function(x){
  return(strsplit(x,"-")[[1]][2])
})
df$species = species

df$start = df$window*20000
df$end = df$start+99999


pdf("~/test.sumstat.pdf")
#for(colnum in 2:34) {
for(colnum in c(1:40,42:339)) {
  print(colnum)
  
  if(sum(complete.cases(as.numeric(df[,colnum])))>0) {
    plot(df$window,df[,colnum],main=colnames(df)[colnum])
  }
  
  
}
dev.off()

