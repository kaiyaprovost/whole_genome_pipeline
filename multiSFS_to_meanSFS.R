#path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/"
#path="/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/"
path="/Users/kprovost/Documents/STAIRWAYPLOT"

setwd(path)

listfiles = list.files(path=path,pattern=".ml$",full.names = T,recursive = T)
overwrite=F

for (i in 1:length(listfiles)) {
  file = listfiles[i]
  newfile = paste(file,".1D.txt",sep="")
  newfile2 = paste(file,".1D.MEAN.txt",sep="")
  print(newfile)
  try({csv = read.table(file,header=F)})
  if(!(file.exists(newfile)) | overwrite==T)
  {
    outline = rbind(colSums(csv,na.rm=T)) ## changed to sums
    write.table(outline,newfile,sep=" ",col.names = F,row.names = F)
  } else {
    print("SUM FILE ALREADY EXISTS -- NOT OVERWRITING")
  }
  if(!(file.exists(newfile2)) | overwrite==T)
  {
    outline2 = rbind(colMeans(csv,na.rm=T)) 
    write.table(outline2,newfile2,sep=" ",col.names = F,row.names = F)
  } else {
    print("MEAN FILE ALREADY EXISTS -- NOT OVERWRITING")
  }
  R.utils::gzip(file)
}
