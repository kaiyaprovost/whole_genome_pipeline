#path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/"
path="/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MLFILES/"

setwd(path)

listfiles = list.files(pattern="MIN1C.1D.txt")

for (i in 1:length(listfiles)) {
  file = listfiles[i]
  csv = read.table(file,header=F)
  outline = rbind(colMeans(csv,na.rm=T))
  write.table(outline,file,sep=" ",col.names = F,row.names = F)
}
