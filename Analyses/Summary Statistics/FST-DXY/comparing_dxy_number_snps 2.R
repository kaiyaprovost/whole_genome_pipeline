dxyalone=dxydf[,c("snps","chr","midPos","species")]
dxyalone$plotorder = paste(dxyalone$chr,dxyalone$midPos)
dxysnptable=as.data.frame(as.matrix(with(dxyalone, tapply(snps, list(plotorder, species), 
                                                     FUN=function(x){mean(x,na.rm=T)}))))
dim(dxysnptable)
dxysnptable=dxysnptable[,neworder]
corrplot::corrplot(cor(dxysnptable,use="pairwise.complete.obs"),
                   method="color",main="snps",
                   col=colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(200))

dxytable=as.data.frame(as.matrix(with(merged, tapply(dxymeans, list(plotorder, species), FUN=function(x){mean(x,na.rm=T)}))))
dim(dxytable)

png("testx.png")
plot(dxysnptable)
dev.off()
