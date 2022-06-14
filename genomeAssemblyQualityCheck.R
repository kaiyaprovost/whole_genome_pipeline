## coverage histograms

bildir = "/Users/kprovost/Documents/Dissertation/genomeresequencingFromLucas/quality_check/raw_wgs_metrics/"
#list.files(bildir,pattern="raw_wgs_metrics.txt$",recursive=T)

#test1 = read.csv("/Users/kprovost/Documents/Dissertation/genomeresequencingFromLucas/quality_check/raw_wgs_metrics/brunneicapillus/AMN_245107_P01_WA02_brunneicapillus.raw_wgs_metrics.txt",
#                 skip=8,sep="\t")
#test2 = read.csv("/Users/kprovost/Documents/Dissertation/genomeresequencingFromLucas/quality_check/raw_wgs_metrics/brunneicapillus/AMN_245107_P01_WF01_brunneicapillus.raw_wgs_metrics.txt",
#                 skip=8,sep="\t")

#x = merge.data.frame(test1,test2,by="coverage",all=T)
#tail(x)

for(filen in list.files(bildir,pattern="raw_wgs_metrics.txt$",recursive=T)){
  test = read.csv(paste(bildir,filen,sep=""),skip=8,sep="\t")
  oldcolnames = colnames(test)
  newcolnames = (paste(oldcolnames,filen,sep="_"))
  colnames(test) = newcolnames
  
  nameonly = unlist(strsplit(unlist(strsplit(filen,split="/"))[2],split=".raw_wgs_metrics.txt"))
  
  colnames(test) = c("coverage",nameonly,nameonly)
  colnames(test)[1] = "coverage"
  
  if(exists("table_hqc")==FALSE){
    table_hqc = test[,c(1,2)]
    table_ufl = test[,c(1,3)]
  } else {
      table_hqc = merge.data.frame(table_hqc,test[,c(1,2)],by="coverage",all=T)
      table_ufl = merge.data.frame(table_ufl,test[,c(1,3)],by="coverage",all=T)
  }
  
  
}

plot(table_hqc$coverage,table_hqc$AMN_245107_P01_WA01_bilineata,
     type="b",xlim=c(0,10))
summary(table_hqc)

table_hqc[is.na(table_hqc)] = 0

for(colnum in 2:22){
  column_name = colnames(table_hqc)[colnum]
  new_name = paste("multiply",column_name,sep="_")
  mult = table_hqc[,1] * table_hqc[,colnum]
  table_hqc = cbind(table_hqc,mult)
  colnames(table_hqc)[ncol(table_hqc)] = new_name
  
}

sum(table_hqc[,43])


x = read.csv(test,skip=8,sep="\t")
x$cumulative = cumsum(x$high_quality_coverage_count)
x$cumul_percent = x$cumulative / sum(x$high_quality_coverage_count)
plot(x$cumul_percent,xlim=c(0,100))
head(x$cumul_percent)

maxmax = 1041168570
maxmax=4e08
y = read.csv("/Users/kprovost/Documents/Dissertation/genomeresequencingFromLucas/quality_check/histograms_raw_wgs.txt",sep="\t")
par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(y$coverage,y$coverage,xlim=c(0,20),ylim=c(0,maxmax),type="n",main="sin")
for(colnum in seq(2,12,2)){
  col = "red"
  pchs = "s"
  points(y$coverage,y[,colnum],xlim=c(0,15),ylim=c(0,4e08),col=col,pch=pchs)
  lines(y$coverage,y[,colnum],xlim=c(0,15),ylim=c(0,4e08),col=col,pch=pchs)
}
plot(y$coverage,y$coverage,xlim=c(0,20),ylim=c(0,maxmax),type="n",main="fus")
for(colnum in seq(14,22,2)){
  col = c(rep("blue",5))
  pchs = c(rep("f",5))
  points(y$coverage,y[,colnum],xlim=c(0,15),ylim=c(0,4e08),col=col,pch=pchs)
  lines(y$coverage,y[,colnum],xlim=c(0,15),ylim=c(0,4e08),col=col,pch=pchs)
}
plot(y$coverage,y$coverage,xlim=c(0,20),ylim=c(0,maxmax),type="n",main="bru")
for(colnum in seq(24,32,2)){
  col = c(rep("green",5))
  pchs = c(rep("u",5))
  points(y$coverage,y[,colnum],xlim=c(0,15),ylim=c(0,4e08),col=col,pch=pchs)
  lines(y$coverage,y[,colnum],xlim=c(0,15),ylim=c(0,4e08),col=col,pch=pchs)
}
plot(y$coverage,y$coverage,xlim=c(0,20),ylim=c(0,maxmax),type="n",main="bil")
for(colnum in seq(24,32,2)){
  col = c(rep("black",5))
  pchs = c(rep("l",5))
  points(y$coverage,y[,colnum],xlim=c(0,15),ylim=c(0,4e08),col=col,pch=pchs)
  lines(y$coverage,y[,colnum],xlim=c(0,15),ylim=c(0,4e08),col=col,pch=pchs)
}

