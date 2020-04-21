data = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/genomes_sequenced_for_map.txt",
                sep="\t")
plot(data$LONG,data$LAT)

data2=read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/genomes_sequenced_for_map_pivot.txt",
               sep="\t")
plot(data2$LONG,data2$LAT,cex=1,col="white",pch=16,
     ylim=c(28.9,35.1),
     xlim=c(-116.1,-97.9))
sizes=scales::rescale(rowSums(data2[,c(3:12)],na.rm=T),
                      to=c(0.1,0.6))

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/genomes_sequenced_for_map_species.png",width=1200)
par(mfrow=c(2,5),mar=c(2,2,0,0))
plot(border,xlim=c(-116,-98),ylim=c(29,35)); points(data2$LONG,data2$LAT,cex=data2$bellii-0.001,col="red",ylab="",xlab="")
plot(border,xlim=c(-116,-98),ylim=c(29,35)); points(data2$LONG,data2$LAT,cex=data2$bilineata-0.001,col="orange",ylab="",xlab="")
plot(border,xlim=c(-116,-98),ylim=c(29,35)); points(data2$LONG,data2$LAT,cex=data2$brunneicapillus-0.001,col="goldenrod",ylab="",xlab="")
plot(border,xlim=c(-116,-98),ylim=c(29,35)); points(data2$LONG,data2$LAT,cex=data2$crissale-0.001,col="green",ylab="",xlab="")
plot(border,xlim=c(-116,-98),ylim=c(29,35)); points(data2$LONG,data2$LAT,cex=data2$curvirostre-0.001,col="cyan",ylab="",xlab="")
plot(border,xlim=c(-116,-98),ylim=c(29,35)); points(data2$LONG,data2$LAT,cex=data2$flaviceps-0.001,col="blue",ylab="",xlab="")
plot(border,xlim=c(-116,-98),ylim=c(29,35)); points(data2$LONG,data2$LAT,cex=data2$fusca-0.001,col="magenta",ylab="",xlab="")
plot(border,xlim=c(-116,-98),ylim=c(29,35)); points(data2$LONG,data2$LAT,cex=data2$melanura-0.001,col="purple",ylab="",xlab="")
plot(border,xlim=c(-116,-98),ylim=c(29,35)); points(data2$LONG,data2$LAT,cex=data2$nitens-0.001,col="pink",ylab="",xlab="")
plot(border,xlim=c(-116,-98),ylim=c(29,35)); points(data2$LONG,data2$LAT,cex=data2$sinuatus-0.001,col="brown",ylab="",xlab="")
dev.off()

colors=c("red","orange","goldenrod","green","cyan",
         "blue","magenta","purple","pink","brown")
jitter_x = c( -1,  0,  1,-1.5,-0.5,0.5,1.5,-1,     0,   1)
jitter_y = c(1.5,1.5,1.5,   0,   0,  0,  0,-1.5,-1.5,-1.5)
for (i in 1:nrow(data2)) {
  for(j in 1:10) {
    points(x=data2$LONG[i]+jitter_x[j]/10,
           y=data2$LAT[i]+jitter_y[j]/10,
           cex=data2[i,j+2]-0.0001,pch=j,
           col=colors[j])
    # plotrix::floating.pie(
    #   xpos = data2$LONG[i]+jitter_x[j]/10,
    #   ypos = data2$LAT[i]+jitter_y[j]/10,
    #   x = c(as.numeric(data2[i,j+2])),
    #   radius = data2[i,j+2]/10,explode=0,
    #   col=colors[j]
    #)
  }
}

segments(x0, y0, x1 = x0, y1 = y0,
         col = par("fg"), lty = par("lty"), lwd = par("lwd"),
         ...)
