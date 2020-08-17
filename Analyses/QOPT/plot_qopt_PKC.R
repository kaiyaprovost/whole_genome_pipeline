library(RColorBrewer)
library(raster)

col = brewer.pal(12,"Paired")
red=col[6]
blue=col[2]
green=col[3]
yellow=col[7]
grey=col[9]

Env = raster::stack('~/Downloads/ENMS_multilayer4.tif')
bg = Env[[1]]
ext = raster::extent(c(-120, -95, 27, 35))
bg = raster::crop(bg, ext)

border = shapefile("/Users/kprovost/Downloads/enm_layers/cb_2016_us_state_500k/cb_2016_us_state_500k.shp")
border2 = shapefile("/Users/kprovost/Downloads/enm_layers/mexstates/mexstates.shp")

qopt=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/AllSpeciesMetadata_allK_9june2020.csv",
                header=T,sep="\t",fill=T,stringsAsFactors = F)
qopt=qopt[order(qopt$LONG,qopt$LAT),]
qopt2 = qopt[,c("SP","JLAT","JLONG","K2_A","K2_B","P2_A","P2_B","C2_A","C2_B",
               "LAT","LONG","STATE")]
qopt3 = qopt[,c("SP","JLAT","JLONG","K3_A","K3_B","K3_C",
                "P3_A","P3_B","P3_C",
                "C3_A","C3_B","C3_C",
               "LAT","LONG","STATE")]

#raster::plot(bg,col = "grey",colNA = "white",legend = F,xaxt="n",yaxt="n",main = "")
#plot(border,add=T,col="grey")
#plot(border2,add=T,col="grey")
#points(qopt$LONG,qopt$LAT,pch=4)

kcols=c(4,5) # K 
pcols=c(6,7) # P
ccols=c(8,9) # C 

kcols3=c(4,5,6) # K 
pcols3=c(7,8,9) # P
ccols3=c(10,11,12) # C 

palette(c("white","red","cyan"))
x=colorRampPalette(c("darkgreen","cyan"))

pdf("strplot_by_species_P_best_reverse.pdf",height=6,width=10)
par(mfrow=c(2,5))
par(mar=c(0,2,1,0))
par(bg="black",col.axis="white",col.lab="white",
    col="white",col.main="white")
spporder = c("crissale","curvirostre","fusca","bellii-NOWEIRD",
             "sinuatus","flaviceps-NOWEIRD","melanura","nitens",
             "brunneicapillus","bilineata")
for(spp in spporder) {
  if(spp %in% c("melanura","brunneicapillus")){
  thissp = (qopt3[qopt3$SP==spp,])
  num_colors=3
  these_cols=pcols3
  } else {
    thissp = (qopt2[qopt2$SP==spp,])
    num_colors=2
    these_cols=pcols
  }
  thissp=thissp[complete.cases(thissp[,these_cols]),]
  barplot(t(as.matrix(thissp[,these_cols])),xaxt="n",
          cex.names=0.5,
          col=rev(x(num_colors)),names.arg=thissp$STATE,las=2,
          main=strsplit(spp,"-")[[1]][1],xlab="P",horiz=T)
}
dev.off()

pdf("strplot_by_species_P_3.pdf",height=30,width=2)
par(mfrow=c(24,1))
par(mar=c(4,0,2,0))
spporder = c("crissale","curvirostre","fusca","bellii-NOWEIRD","bellii",
             "sinuatus","flaviceps-NOWEIRD","flaviceps","melanura","nitens",
             "brunneicapillus","bilineata")
for(spp in spporder) {
  thissp = (qopt3[qopt3$SP==spp,])
  thissp=thissp[complete.cases(thissp[,pcols]),]
  if(nrow(thissp)==0){plot.new()} else {
  barplot(t(as.matrix(thissp[,pcols3])),yaxt="n",
          col=x(3),names.arg=thissp$STATE,las=2,
          main=spp,xlab="P")
  }
  thissp = (qopt2[qopt2$SP==spp,])
  thissp=thissp[complete.cases(thissp[,pcols]),]
  if(nrow(thissp)==0){plot.new()} else {
  barplot(t(as.matrix(thissp[,pcols])),yaxt="n",
          col=x(2),names.arg=thissp$STATE,las=2,
          main=spp,xlab="P")
  }
}
dev.off()

pdf("strplot_by_species_P_2.pdf",height=8,width=5)
par(mfrow=c(5,2))
par(mar=c(4,0,2,0))
spporder = c("crissale","curvirostre","fusca","bellii-NOWEIRD",
             "sinuatus","flaviceps-NOWEIRD","melanura","nitens",
             "brunneicapillus","bilineata")
for(spp in spporder) {
  thissp = (qopt2[qopt2$SP==spp,])
  thissp=thissp[complete.cases(thissp[,pcols]),]
  barplot(t(as.matrix(thissp[,pcols])),yaxt="n",
          col=x(2),names.arg=thissp$STATE,las=2,
          main=strsplit(spp,"-")[[1]][1],xlab="P")
}
dev.off()

pdf("strplot_by_species_KPC_3.pdf",height=15,width=8.5)
par(mfrow=c(12,2))
par(mar=c(4,0,2,0))
for(spp in sort(unique(qopt$SP))) {
  thissp = (qopt3[qopt3$SP==spp,])
  barplot(t(as.matrix(thissp[,kcols3])),yaxt="n",
          col=x(3),names.arg=thissp$STATE,las=2,
          main=spp,xlab="K")
  barplot(t(as.matrix(thissp[,pcols3])),yaxt="n",
          col=x(3),names.arg=thissp$STATE,las=2,
          main=spp,xlab="P")
  #barplot(t(as.matrix(thissp[,ccols3])),yaxt="n",
  #        col=x(3),names.arg=thissp$STATE,las=2,
  #        xlab="C")
}
dev.off()

pdf("strplot_by_species_KPC_2.pdf",height=15,width=8.5)
par(mfrow=c(12,2))
par(mar=c(4,0,2,0))
for(spp in sort(unique(qopt$SP))) {
  thissp = (qopt2[qopt2$SP==spp,])
  barplot(t(as.matrix(thissp[,kcols])),yaxt="n",
          col=x(2),names.arg=thissp$STATE,las=2,
          xlab="K")
  barplot(t(as.matrix(thissp[,pcols])),yaxt="n",
          col=x(2),names.arg=thissp$STATE,las=2,
          main=spp,xlab="P")
  #barplot(t(as.matrix(thissp[,ccols])),yaxt="n",
  #        col=x(2),names.arg=thissp$STATE,las=2,
  #        xlab="C")
}
dev.off()

pdf("qopt_by_species_KPC.pdf",height=6,width=6)
for(spp in sort(unique(qopt$SP))) {
par(mfrow=c(3,1),mar=c(0,0,0,0))
  print(spp)
  for(loop in 1:3){
    thesecols = c(2*loop+2,2*loop+3)
    thisdata = c("K","P","C")[loop]
    print(thisdata)
    
    #raster::plot(bg,col = "grey",colNA = "white",legend = F,xaxt="n",yaxt="n")
    plot(border,add=F,col="grey",xlim=c(-120,-95),ylim=c(27,35))
    plot(border2,add=T,col="grey")
    legend("bottomleft",legend = paste(spp,thisdata))
    subsetqopt=qopt[qopt$SP==spp,]
    for (i in 1:nrow(subsetqopt)) {
      if(!(is.na(subsetqopt[i,thesecols[1]]))) {
        #cat(i,sep=" ")
        plotrix::floating.pie(xpos = subsetqopt$JLONG[i],ypos = subsetqopt$JLAT[i],
                              x = c(as.numeric(subsetqopt[i,c(thesecols)])),radius = 0.3,
                              col = c("#1F78B4","#B2DF8A","#E31A1C","#FDBF6F","#CAB2D6"))
      }
  }
  

  }
  
}
dev.off()


png("qopt_by_species_KPC_cur-sin-bil.png",height=4,width=12,units="in",res=600)
par(mfrow=c(3,3),mar=c(0,0,0,0))
print(spp)
m <- rbind(c(1,4,7), c(2,5,8), c(3,6,9))
#print(m)
layout(m)
#layout.show(9)
for(spp in c("curvirostre","sinuatus","bilineata")) {
  for(loop in 1:3){
    thesecols = c(2*loop+2,2*loop+3)
    thisdata = c("K","P","C")[loop]
    print(thisdata)
    
    #raster::plot(bg,col = "grey",colNA = "white",legend = F,xaxt="n",yaxt="n")
    plot(border,add=F,col="grey",xlim=c(-120,-95),ylim=c(27,35))
    plot(border2,add=T,col="grey")
    legend("bottomleft",legend = paste(spp,thisdata))
    subsetqopt=qopt[qopt$SP==spp,]
    for (i in 1:nrow(subsetqopt)) {
      if(!(is.na(subsetqopt[i,thesecols[1]]))) {
        #cat(i,sep=" ")
        plotrix::floating.pie(xpos = subsetqopt$JLONG[i],ypos = subsetqopt$JLAT[i],
                              x = c(as.numeric(subsetqopt[i,c(thesecols)])),radius = 0.3,
                              col = c("#1F78B4","#B2DF8A","#E31A1C","#FDBF6F","#CAB2D6"))
      }
    }
    
    
  }
}
dev.off()

## k=3
qopt=read.table("/Users/kprovost/Downloads/QOPT/AllSpeciesMetadata_allK.csv",
                header=T,sep="\t",fill=T,stringsAsFactors = F)
qopt = qopt[,c("SP","JLAT","JLONG","K3_A","K3_B","K3_C","P3_A","P3_B","P3_C","C3_A","C3_B","C3_C",
               "LAT","LONG")]

kcols=c(4:6) # K 
pcols=c(7:9) # P
ccols=c(10:12) # C 

pdf("qopt_by_species_KPC-3.pdf",height=6,width=6)
for(spp in sort(unique(qopt$SP))) {
  par(mfrow=c(3,1),mar=c(0,0,0,0))
  print(spp)
  for(loop in 1:3){
    thesecols = c(3*loop+1,3*loop+2,3*loop+3)
    thisdata = c("K","P","C")[loop]
    print(thisdata)
    
    #raster::plot(bg,col = "grey",colNA = "white",legend = F,xaxt="n",yaxt="n")
    plot(border,add=F,col="grey",xlim=c(-120,-95),ylim=c(27,35))
    plot(border2,add=T,col="grey")
    legend("bottomleft",legend = paste(spp,thisdata))
    subsetqopt=qopt[qopt$SP==spp,]
    for (i in 1:nrow(subsetqopt)) {
      if(!(is.na(subsetqopt[i,thesecols[1]]))) {
        #cat(i,sep=" ")
        plotrix::floating.pie(xpos = subsetqopt$JLONG[i],ypos = subsetqopt$JLAT[i],
                              x = c(as.numeric(subsetqopt[i,c(thesecols)])),radius = 0.3,
                              col = c("#1F78B4","#B2DF8A","#E31A1C","#FDBF6F","#CAB2D6"))
      }
    }
    
    
  }
  
}
dev.off()


