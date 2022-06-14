df = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/test_missing_by_lat.txt",
                sep="\t",header=T)

df$col=round(df$F_MISS*100)

pal <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))
df$colors=pal(100)[df$col]
palette(pal(100))

pdf("missing_jlatjlong_species.pdf",height=8,width=8)
par(mfrow=c(2,2),mar=c(0,0,2,0))
for(spp in unique(df$SP)){
  print(spp)
  temp=df[df$SP==spp,]
  plot(temp$JLONG,temp$JLAT,cex=4,pch=21,
       bg=temp$colors,main=spp,xaxt="n",yaxt="n")
  text(x=temp$JLONG,y=temp$JLAT,labels = temp$col)

  plot(temp$JLONG,temp$JLAT,cex=4,pch=21,
       bg=temp$colors,main="K",type="n",xaxt="n",yaxt="n")
  for (i in 1:nrow(temp)) {
    if(!(is.na(temp[i,c(7,8,9)]))) {
      plotrix::floating.pie(xpos = temp$JLONG[i],ypos = temp$JLAT[i],
                            x = c(as.numeric(temp[i,c(7,8,9)])),radius = 0.3,
                            col = c("#1F78B4","#B2DF8A","#E31A1C","#FDBF6F","#CAB2D6"))
    }
  }
    
    plot(temp$JLONG,temp$JLAT,cex=4,pch=21,
         bg=temp$colors,main="P",type="n",xaxt="n",yaxt="n")
    for (i in 1:nrow(temp)) {
      if(!(is.na(temp[i,c(12,13,14)]))) {
        plotrix::floating.pie(xpos = temp$JLONG[i],ypos = temp$JLAT[i],
                              x = c(as.numeric(temp[i,c(12,13,14)])),radius = 0.3,
                              col = c("#1F78B4","#B2DF8A","#E31A1C","#FDBF6F","#CAB2D6"))
      }
    }
    
    plot(temp$JLONG,temp$JLAT,cex=4,pch=21,
         bg=temp$colors,main="C",type="n",xaxt="n",yaxt="n")
    for (i in 1:nrow(temp)) {
      if(!(is.na(temp[i,c(17,18,19)]))) {
        plotrix::floating.pie(xpos = temp$JLONG[i],ypos = temp$JLAT[i],
                              x = c(as.numeric(temp[i,c(17,18,19)])),radius = 0.3,
                              col = c("#1F78B4","#B2DF8A","#E31A1C","#FDBF6F","#CAB2D6"))
      }
    }
  
  }
dev.off()


pdf("missing_jlatjlong_species_K2.pdf",height=8,width=8)
par(mfrow=c(2,2),mar=c(0,0,2,0))
for(spp in unique(df$SP)){
  print(spp)
  temp=df[df$SP==spp,]
  plot(temp$JLONG,temp$JLAT,cex=4,pch=21,
       bg=temp$colors,main=spp,xaxt="n",yaxt="n")
  text(x=temp$JLONG,y=temp$JLAT,labels = temp$col)
  
  plot(temp$JLONG,temp$JLAT,cex=4,pch=21,
       bg=temp$colors,main="K",type="n",xaxt="n",yaxt="n")
  for (i in 1:nrow(temp)) {
    if(!(is.na(temp[i,c(5,6)]))) {
      plotrix::floating.pie(xpos = temp$JLONG[i],ypos = temp$JLAT[i],
                            x = c(as.numeric(temp[i,c(5,6)])),radius = 0.3,
                            col = c("#1F78B4","#B2DF8A","#E31A1C","#FDBF6F","#CAB2D6"))
    }
  }
  
  plot(temp$JLONG,temp$JLAT,cex=4,pch=21,
       bg=temp$colors,main="P",type="n",xaxt="n",yaxt="n")
  for (i in 1:nrow(temp)) {
    if(!(is.na(temp[i,c(10,11)]))) {
      plotrix::floating.pie(xpos = temp$JLONG[i],ypos = temp$JLAT[i],
                            x = c(as.numeric(temp[i,c(10,11)])),radius = 0.3,
                            col = c("#1F78B4","#B2DF8A","#E31A1C","#FDBF6F","#CAB2D6"))
    }
  }
  
  plot(temp$JLONG,temp$JLAT,cex=4,pch=21,
       bg=temp$colors,main="C",type="n",xaxt="n",yaxt="n")
  for (i in 1:nrow(temp)) {
    if(!(is.na(temp[i,c(15,16)]))) {
      plotrix::floating.pie(xpos = temp$JLONG[i],ypos = temp$JLAT[i],
                            x = c(as.numeric(temp[i,c(15,16)])),radius = 0.3,
                            col = c("#1F78B4","#B2DF8A","#E31A1C","#FDBF6F","#CAB2D6"))
    }
  }
  
}
dev.off()

nona=df[complete.cases(df),]

mod=lme4::lmer(nona$F_MISS~nona$K2_A+(1|nona$SP))
summary(mod)

mod=lme4::lmer(nona$F_MISS~nona$P2_A+(1|nona$SP))
summary(mod)

mod=lme4::lmer(nona$F_MISS~nona$C2_A+(1|nona$SP))
summary(mod)
