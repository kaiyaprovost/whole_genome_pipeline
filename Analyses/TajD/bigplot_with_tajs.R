setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas")

bigplot = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.may2020.txt")
#bigplot2 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.may2020.temp")

names(bigplot2)

bigtaj = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigtaj.new.may2020.txt")
colnames(bigtaj) = c("windowinfo",
                     "chr",
                     "midPos",
                     "tW",
                     "tP-pi",
                     "tF",
                     "tH",
                     "tL",
                     "Tajima",
                     "fuf",
                     "fud",
                     "fayh",
                     "zeng",
                     "nSites",
                     "species"
)

megaplot = merge(bigtaj,bigplot,all=T)
megaplot=unique(megaplot)

write.csv(megaplot,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.fst.dxy.taj.may2020.txt",
          row.names = F)

#bigplot2 = bigplot2[order(bigplot2$chr, bigplot2$midPos), ]
#write.csv(bigplot2,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot2.txt",
#          row.names = F)

#corrplot::corrplot(cor(bigplot2[,c(5:8,10:23,25:35)]),method="ellipse",
#                   order="hclust")



plotQuanAndStdev = function(variable=temp$Fst,chrom=temp$chr,
                            speciesname=spp,varname=names(variable),
                            printed=T){
  
  #chrom=chrom[!(is.na(variable))]
  #variable=variable[!(is.na(variable))]
  
  ## check if all of the variable is na
  are_na = unique(is.na(variable))
  
  if ( length(are_na) <= 1 && sum(are_na) != 0) {
    
    if(printed==T){
      png(paste(speciesname,"_",varname,"_5stdev.may2020.png",sep=""),width=800,height=400)
    }
    plot(0,ylab=varname,xlab="scaffold",col="white",type="n")
    if(printed==T){
      dev.off()
    }
    
  } else {
    
    meanfst = mean(variable,na.rm=T)
    stdevfst = sd(variable,na.rm=T)
    high5sd = meanfst+(5*stdevfst)
    low5sd = meanfst-(5*stdevfst)
    quan95 = quantile(variable,0.95,na.rm=T)
    quan5 = quantile(variable,0.05,na.rm=T)
    
    if(printed==T){
      png(paste(speciesname,"_",varname,"_5stdev.may2020.png",sep=""),width=800,height=400)
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
  
  
}


for (spp in sort(as.character(unique(megaplot$species)))) {
  print(spp)
  
  #temp=bigplot2[bigplot2$species==spp,]
  temp=megaplot[megaplot$species==spp,]
  
  png(paste(spp,"_fourpanels.may2020.png",sep=""),width=1600,height=1000)
  palette(c("red","orange","goldenrod","green","cyan","blue","magenta","purple"))
  par(mfrow=c(4,1))
  #plot(temp$Fst,col=as.numeric(as.factor(temp$chr)),xlab="Fst")
  #plot(temp$dxymeans,col=as.numeric(as.factor(temp$chr)),xlab="Dxy")
  #plot(temp$tW,col=as.numeric(as.factor(temp$chr)),xlab="Theta")
  #plot(temp$Tajima,col=as.numeric(as.factor(temp$chr)),xlab="TajD")
  print(1)
  plotQuanAndStdev(variable=temp$Fst,chrom=temp$chr,speciesname=spp,varname="Fst",printed=F)
  print(2)
  plotQuanAndStdev(variable=temp$dxymeans,chrom=temp$chr,speciesname=spp,varname="Dxy",printed=F)
  print(3)
  plotQuanAndStdev(variable=temp$tW,chrom=temp$chr,speciesname=spp,varname="Theta",printed=F)
  print(4)
  plotQuanAndStdev(variable=temp$Tajima,chrom=temp$chr,speciesname=spp,varname="TajD",printed=F)
  dev.off()
  
  
  ## plot Taj vs Theta 
  quan95taj = quantile(temp$Tajima,0.95,na.rm=T)
  quan5taj = quantile(temp$Tajima,0.05,na.rm=T)
  quan95th = quantile(temp$tW,0.95,na.rm=T)
  quan5th = quantile(temp$tW,0.05,na.rm=T)
  
  are_na1 = unique(is.na(temp$Tajima))
  are_na2 = unique(is.na(temp$tW))
  
  if((length(are_na1) <= 1 && sum(are_na1) != 0) || (length(are_na2) <= 1 && sum(are_na2) != 0 )) {
    
    print("no data")
    
  } else {
  
  png(paste(spp,"_taj_vs_theta.may2020.png",sep=""))
  plot(temp$tW,temp$Tajima,xlab="Theta",ylab="Tajima")
  abline(h=quan95taj,col="red")
  abline(h=quan5taj,col="red")
  abline(v=quan95th,col="red")
  abline(v=quan5th,col="red")
  dev.off()
  }
  
}


## plot Taj vs Theta 

megaplot =read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.fst.dxy.taj.may2020.txt",
                     header = T)
megaplot$plotorder = paste(megaplot$chr,megaplot$midPos)
uniques=unique(megaplot$plotorder)
numuniques = uniques

lapply(unique(megaplot$plotorder),FUN=function(x){
  print(x)
  megaplot$numorder[megaplot$plotorder==x] = megaplot$plotorder[megaplot$plotorder==x] = which(uniques %in% x)
})

megaplot$plotorder = factor(megaplot$plotorder, levels = unique(as.character(megaplot$plotorder)))
megaplot$plotorder = as.numeric(as.factor(megaplot$plotorder))

# for(i in 1:length(uniques)){
#   if(i %% 1000 == 0){print(i/length(uniques))}
#   value=uniques[i]
#   found_index=which(uniques %in% value)
#   numuniques[i]=found_index
# }
# megaplot$plotorder = numuniques

write.csv(megaplot,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.fst.dxy.taj.may2020.txt",
                   row.names = F)

megaplot =read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.fst.dxy.taj.may2020.txt",
                   header = T,sep="\t")

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/all_species_lined_up_fst.may2020.png",height=1500,width=800)
par(mfrow=c(10,1),mar=c(4,4,0,0))
for (spp in sort(as.character(unique(megaplot$species)))) {
  print(spp)
  thisspp = megaplot[megaplot$species==spp,]
  plot(megaplot$plotorder,megaplot$Fst,col=rgb(0.9,0.9,0.9),
       xlab="window",ylab="Fst")
  points(thisspp$plotorder,thisspp$Fst,col=as.numeric(as.factor(thisspp$chr)),add=T)
}
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/all_species_lined_up_dxy.may2020.png",height=1500,width=800)
par(mfrow=c(10,1),mar=c(4,4,0,0))
for (spp in sort(as.character(unique(megaplot$species)))) {
  print(spp)
  thisspp = megaplot[megaplot$species==spp,]
  plot(megaplot$plotorder,megaplot$dxymeans,col=rgb(0.9,0.9,0.9),
       xlab="window",ylab="Dxy")
  points(thisspp$plotorder,thisspp$dxymeans,col=as.numeric(as.factor(thisspp$chr)),add=T)
}
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/all_species_lined_up_tajd.may2020.png",height=1500,width=800)
par(mfrow=c(10,1),mar=c(4,4,0,0))
for (spp in sort(as.character(unique(megaplot$species)))) {
  print(spp)
  thisspp = megaplot[megaplot$species==spp,]
  plot(megaplot$plotorder,megaplot$Tajima,col=rgb(0.9,0.9,0.9),
       xlab="window",ylab="Tajima's D")
  points(thisspp$plotorder,thisspp$Tajima,col=as.numeric(as.factor(thisspp$chr)),add=T)
}
dev.off()

## check which chroms are missing for which species and which analyses
all_chr = unique(megaplot$chr)
smallmega = unique(megaplot[,c("chr","species","Tajima","dxymeans","Fst")])
smallmega$Tajima[!is.na(smallmega$Tajima)] = 1
smallmega$dxymeans[!is.na(smallmega$dxymeans)] = 1
smallmega$Fst[!is.na(smallmega$Fst)] = 1

smallmega = unique(smallmega)

taj = unique(smallmega[,c("chr","species","Tajima")])
dxy = unique(smallmega[,c("chr","species","dxymeans")])
fst = unique(smallmega[,c("chr","species","Fst")])

(table(taj))
(table(dxy))
(table(fst))
