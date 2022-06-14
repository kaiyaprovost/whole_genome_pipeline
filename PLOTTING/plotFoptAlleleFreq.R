library(R.utils)
library(tools)
#gunzip("file.gz", remove=FALSE)

path = "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/FOPTS"

setwd(path)

filelist = list.files(path=path,pattern="fopt.gz$",recursive = T,full.names = T) 

for(gzfile in filelist) {
  
  gunzip(gzfile,remove=FALSE)
  filename = tools::file_path_sans_ext(gzfile)
  
  #filename = "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/FOPTS/test.fopt"
  filename="~/Downloads/test.fopt.txt"
  split = strsplit(filename,"/")[[1]]
  basefile = split[length(split)]
  #K = strsplit(basefile,"\\.")[[1]][1]
  K=1
  
  allelefreq = read.table(filename,sep=" ")
  allelefreq = allelefreq[,1:K]
  palette(c("black","red","blue","goldenrod","green"))
  barplot(t(as.matrix(allelefreq)),beside=T,col=c(1:5))
  
  plot(allelefreq,col=1,lty=1,lwd=2,type="l")
  points(allelefreq$V2,col=2,lty=3,lwd=2,type="l")
  points(allelefreq$V3,col=3,lty=3,lwd=2,type="l")
  points(allelefreq$V4,col=4,lty=3,lwd=2,type="l")
  points(allelefreq$V5,col=5,lty=3,lwd=2,type="l")
}

