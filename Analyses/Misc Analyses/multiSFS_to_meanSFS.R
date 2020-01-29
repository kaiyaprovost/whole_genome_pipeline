path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/"
setwd(path)

listfiles = list.files(pattern="ml")

for (i in 2:length(listfiles)) {
  file = listfiles[i]
  csv = read.csv(file,sep="\t",header=F)
  outline = rbind(colMeans(csv,na.rm=T))
  write.table(outline,file,sep=" ",col.names = F,row.names = F)
}
