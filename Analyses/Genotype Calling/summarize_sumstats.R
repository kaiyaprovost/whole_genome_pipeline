x = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/Amphispiza-bilineata-called.geno/VCFS/Test_emprical_data_1tab.txt",sep="\t")

for (i in c(2:75)) {
  png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/Amphispiza-bilineata-called.geno/VCFS/",names(x)[i],".png",sep=""))
  plot(log(x[,76]),x[,i],ylab=names(x)[i],xlab=paste("Log",names(x)[76]))
  dev.off()
}
