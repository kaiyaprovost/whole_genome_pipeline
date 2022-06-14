## take in a mlfile from the command line

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Cannot run without input ML file", call.=FALSE)
}
mlfile = args[1]
#mlfile = "/Users/kprovost/Downloads/SON_CHI_Vireo_bellii_75.ml"
mldf = read.table(mlfile,header=F,sep=" ",strip.white = T)
mlsum = colSums(mldf,na.rm=T)
mlsum = mlsum[complete.cases(mlsum)]
mlsum = rbind(mlsum)

if (mlsum[ncol(mlsum)]==0){
  mlsum=mlsum[-(ncol(mlsum))]
  mlsum = rbind(mlsum)
}

write.table(mlsum,file=paste(mlfile,".1D.txt",sep=""),row.names = F,sep=" ",col.names=F)
