#path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/"
path="/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MLSTATS"

setwd(path)

listfiles = list.files(path=path,pattern=".ml$",full.names = T,recursive = F)
overwrite=F

for (i in 1:length(listfiles)) {
  file = listfiles[i]
  newfile = paste(file,".1D.txt",sep="")
  print(newfile)
  if(!(file.exists(newfile)) | overwrite==T)
  {
    csv = read.table(file,header=F)
    outline = rbind(colSums(csv,na.rm=T)) ## changed to sums
    write.table(outline,newfile,sep=" ",col.names = F,row.names = F)
  } else {
    print("FILE ALREADY EXISTS -- NOT OVERWRITING")
  }
  R.utils::gzip(file)
}
