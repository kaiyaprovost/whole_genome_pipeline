merged=read.table("~/rec_taj_dxy_fst_islswp_missing.temp",fill=T,header=T,stringsAsFactors = F)
pdf("plots.pdf")
numbers=c(4,5,8,9,23)
for(i in numbers){
  for(j in numbers){
    if(i < j){
      mod=lm(merged[,j]~merged[,i])
      png(paste(colnames(merged)[i],colnames(merged)[j],"png",sep = "."))
      plot(merged[,i],merged[,j],
           ylab=colnames(merged)[j],xlab=colnames(merged)[i],
           main=summary(mod)$adj.r.squared,
           sub=summary(mod)$coefficients[2,4])
      abline(mod,col="red")
      dev.off()
    }
  }
}
dev.off()


small = merged[,c("species","Fst","dxymeans","ranksppSWEEP","ranksppISLAND")]
small=small[complete.cases(small),]
small=unique(small)
small$type = 1
small$type[small$ranksppISLAND==1] = 3
small$type[small$ranksppSWEEP==1] = 2
for(spp in unique(small$species)){
  sub = small[small$species==spp,]
  png(paste(spp,"dxyfst plot.png"))
  par(col="white",bg="black",col.axis="white",
      col.main="white",col.sub="white")
  plot(sub$dxymeans,sub$Fst,col=c("white","magenta","cyan")[as.numeric(small$type)],
       cex=as.numeric(small$type),pch=16,main=spp)
  abline(h=quantile(sub$Fst,c(0.05,0.95)),col="red")
  abline(v=quantile(sub$dxymeans,c(0.05,0.95)),col="red")
  points(sub$dxymeans[sub$type==2],sub$Fst[sub$type==2],col="magenta",cex=2,pch=16)
  points(sub$dxymeans[sub$type==3],sub$Fst[sub$type==3],col="cyan",cex=3,pch=16)
  
  dev.off()
}

corrplot::corrplot(cor(merged[,c(4,5,8,9,23)],use="pairwise.complete.obs"),method="color")
