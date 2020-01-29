setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas")

bigplot2 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot2.txt")

names(bigplot2)

#bigplot2 = bigplot2[order(bigplot2$chr, bigplot2$midPos), ]
#write.csv(bigplot2,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot2.txt",
#          row.names = F)

corrplot::corrplot(cor(bigplot2[,c(5:8,10:23,25:35)]),method="ellipse",
                   order="hclust")

plotQuanAndStdev = function(variable=temp$Fst,chrom=temp$chr,
                            speciesname=spp,varname=names(variable),
                            printed=T){
  
  meanfst = mean(variable)
  stdevfst = sd(variable)
  high5sd = meanfst+(5*stdevfst)
  low5sd = meanfst-(5*stdevfst)
  quan95 = quantile(variable,0.95)
  quan5 = quantile(variable,0.05)
  
  if(printed==T){
  png(paste(speciesname,"_",varname,"_5stdev.png",sep=""),width=800,height=400)
  }
  plot(variable,col=as.numeric(as.factor(chrom)),ylab=varname,xlab="scaffold")
  abline(h=high5sd,col="black")
  abline(h=low5sd,col="black")
  abline(h=quan95,col="grey")
  abline(h=quan5,col="grey")
  if(printed==T){
  dev.off()
  }
}


for (spp in sort(unique(bigplot2$species))) {
  print(spp)
  
  temp=bigplot2[bigplot2$species==spp,]
  
  # png(paste(spp,"_fourpanels.png",sep=""),width=1600,height=1000)
  # palette(c("red","orange","goldenrod","green","cyan","blue","magenta","purple"))
  # par(mfrow=c(4,1))
  # #plot(temp$Fst,col=as.numeric(as.factor(temp$chr)),xlab="Fst")
  # #plot(temp$dxymeans,col=as.numeric(as.factor(temp$chr)),xlab="Dxy")
  # #plot(temp$tW,col=as.numeric(as.factor(temp$chr)),xlab="Theta")
  # #plot(temp$Tajima,col=as.numeric(as.factor(temp$chr)),xlab="TajD")
  # plotQuanAndStdev(variable=temp$Fst,chrom=temp$chr,speciesname=spp,varname="Fst",printed=F)
  # plotQuanAndStdev(variable=temp$dxymeans,chrom=temp$chr,speciesname=spp,varname="Dxy",printed=F)
  # plotQuanAndStdev(variable=temp$tW,chrom=temp$chr,speciesname=spp,varname="Theta",printed=F)
  # plotQuanAndStdev(variable=temp$Tajima,chrom=temp$chr,speciesname=spp,varname="TajD",printed=F)
  # dev.off()
  
  
  ## plot Taj vs Theta 
  quan95taj = quantile(temp$Tajima,0.95)
  quan5taj = quantile(temp$Tajima,0.05)
  quan95th = quantile(temp$tW,0.95)
  quan5th = quantile(temp$tW,0.05)
  
  png(paste(spp,"_taj_vs_theta.png",sep=""))
  plot(temp$tW,temp$Tajima,xlab="Theta",ylab="Tajima")
  abline(h=quan95taj,col="red")
  abline(h=quan5taj,col="red")
  abline(v=quan95th,col="red")
  abline(v=quan5th,col="red")
  dev.off()
    
}


## plot Taj vs Theta 
