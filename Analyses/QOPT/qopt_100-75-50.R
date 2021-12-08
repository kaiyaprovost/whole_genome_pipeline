df = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/qopts_1/AllSpeciesMetadata_lostruct_missing_19oct2021.csv",
              sep=",",header=T)
df = df[order(df$LONG),]
head(df)

wcdata = raster::getData("worldclim",var="bio",res=10)
wcdata = wcdata[[c(1)]]
wcdata = raster::crop(wcdata,raster::extent(c(min(df$LONG,na.rm=T),
                                              max(df$LONG,na.rm=T),
                                              min(df$LAT,na.rm=T),
                                              max(df$LAT,na.rm=T))
))
raster::plot(wcdata)

states = raster::shapefile("/Users/kprovost/Documents/OneDrive - The Ohio State University/Environment/Environmental_Layers_Dissertation/cb_2016_us_state_500k/cb_2016_us_state_500k.shp")
mex = raster::shapefile("/Users/kprovost/Documents/OneDrive - The Ohio State University/Environment/Environmental_Layers_Dissertation/mexstates/mexstates.shp")
# raster::plot(states,ylim=c(min(df$LAT,na.rm=T), max(df$LAT,na.rm=T)),
#              xlim=c(min(df$LONG,na.rm=T), max(df$LONG,na.rm=T)),col="grey")
# raster::plot(mex,ylim=c(min(df$LAT,na.rm=T), max(df$LAT,na.rm=T)),
#              xlim=c(min(df$LONG,na.rm=T), max(df$LONG,na.rm=T)),add=T,col="grey")


draw_pies = function(dataset,datacols,colors=c("blue","green","red"),longcol="JLONG",latcol="JLAT",radius=0.3){

  longs = dataset[,longcol]
  lats = dataset[,latcol]
  data = dataset[,datacols]
  for(i in 1:nrow(dataset)){
    plotrix::floating.pie(xpos=longs[i],ypos=lats[i],
                          x=as.numeric(data[i,]),
                          radius=radius*sum(as.numeric(data[i,]),na.rm=T),
                          col=colors)
  }
  
}


## pie charts
for(spp in sort(unique(df$Species))){
  
  temp = df[df$Species==spp,]
  temp_K2 = temp[,c("Species","Dataset","JLAT","JLONG","K2.A","K2.B")]
  temp_K3 = temp[,c("Species","Dataset","JLAT","JLONG","K3.A","K3.B","K3.C")]
  temp_P2 = temp[,c("Species","Dataset","JLAT","JLONG","P2.A","P2.B")]
  temp_P3 = temp[,c("Species","Dataset","JLAT","JLONG","P3.A","P3.B","P3.C")]

  temp_K2=temp_K2[complete.cases(temp_K2),]
  temp_K3=temp_K3[complete.cases(temp_K3),]
  temp_P2=temp_P2[complete.cases(temp_P2),]
  temp_P3=temp_P3[complete.cases(temp_P3),]
  
  if(nrow(temp_K2)>0){
   for(dataset in sort(unique(temp_K2$Dataset))){
     temp_K2_dat = temp_K2[temp_K2$Dataset==dataset,]
     
     png(paste(spp,"_",dataset,"-Missing_Pies_K2.png",sep=""))
     raster::plot(states,ylim=c(min(df$LAT,na.rm=T), max(df$LAT,na.rm=T)),
                  xlim=c(min(df$LONG,na.rm=T), max(df$LONG,na.rm=T)),col="grey")
     raster::plot(mex,ylim=c(min(df$LAT,na.rm=T), max(df$LAT,na.rm=T)),
                  xlim=c(min(df$LONG,na.rm=T), max(df$LONG,na.rm=T)),add=T,col="grey")
     draw_pies(dataset=temp_K2_dat,datacols=c("K2.A","K2.B"),colors = c("#4977b1","#bfde8e","#c32022"))
     dev.off()
     
   } 
  }
  
  if(nrow(temp_K3)>0){
    for(dataset in sort(unique(temp_K3$Dataset))){
      temp_K3_dat = temp_K3[temp_K3$Dataset==dataset,]
      
      png(paste(spp,"_",dataset,"-Missing_Pies_K3.png",sep=""))
      raster::plot(states,ylim=c(min(df$LAT,na.rm=T), max(df$LAT,na.rm=T)),
                   xlim=c(min(df$LONG,na.rm=T), max(df$LONG,na.rm=T)),col="grey")
      raster::plot(mex,ylim=c(min(df$LAT,na.rm=T), max(df$LAT,na.rm=T)),
                   xlim=c(min(df$LONG,na.rm=T), max(df$LONG,na.rm=T)),add=T,col="grey")
      draw_pies(dataset=temp_K3_dat,datacols=c("K3.A","K3.B","K3.C"),colors = c("#4977b1","#bfde8e","#c32022"))
      dev.off()
      
    } 
  }
  
  if(nrow(temp_P2)>0){
    for(dataset in sort(unique(temp_P2$Dataset))){
      temp_P2_dat = temp_P2[temp_P2$Dataset==dataset,]
      
      png(paste(spp,"_",dataset,"-Missing_Pies_P2.png",sep=""))
      raster::plot(states,ylim=c(min(df$LAT,na.rm=T), max(df$LAT,na.rm=T)),
                   xlim=c(min(df$LONG,na.rm=T), max(df$LONG,na.rm=T)),col="grey")
      raster::plot(mex,ylim=c(min(df$LAT,na.rm=T), max(df$LAT,na.rm=T)),
                   xlim=c(min(df$LONG,na.rm=T), max(df$LONG,na.rm=T)),add=T,col="grey")
      draw_pies(dataset=temp_P2_dat,datacols=c("P2.A","P2.B"),colors = c("#4977b1","#bfde8e","#c32022"))
      dev.off()
      
    } 
  }
  
  if(nrow(temp_P3)>0){
    for(dataset in sort(unique(temp_P3$Dataset))){
      temp_P3_dat = temp_P3[temp_P3$Dataset==dataset,]
      
      png(paste(spp,"_",dataset,"-Missing_Pies_P3.png",sep=""))
      raster::plot(states,ylim=c(min(df$LAT,na.rm=T), max(df$LAT,na.rm=T)),
                   xlim=c(min(df$LONG,na.rm=T), max(df$LONG,na.rm=T)),col="grey")
      raster::plot(mex,ylim=c(min(df$LAT,na.rm=T), max(df$LAT,na.rm=T)),
                   xlim=c(min(df$LONG,na.rm=T), max(df$LONG,na.rm=T)),add=T,col="grey")
      draw_pies(dataset=temp_P3_dat,datacols=c("P3.A","P3.B","P3.C"),colors = c("#4977b1","#bfde8e","#c32022"))
      dev.off()
      
    } 
  }
    
}

for(spp in sort(unique(df$Species))){
  print(spp)
  temp = df[df$Species==spp,]
  
  palette(c("red","blue","green"))
  
  # try({mod.k2.100.A = glm(temp$K2.A[temp$Dataset==100]~temp$LONG[temp$Dataset==100]); print(paste(spp,"K2.A 100 adjR2:",round(summary(mod.k2.100.A)$adj.r.squared,3)));},silent = T)
  # try({mod.k2.75.A = glm(temp$K2.A[temp$Dataset==75]~temp$LONG[temp$Dataset==75]); print(paste(spp,"K2.A 75 adjR2:",round(summary(mod.k2.75.A)$adj.r.squared,3)))},silent = T)
  # try({mod.k2.50.A = glm(temp$K2.A[temp$Dataset==50]~temp$LONG[temp$Dataset==50]); print(paste(spp,"K2.A 50 adjR2:",round(summary(mod.k2.50.A)$adj.r.squared,3)))},silent = T)
  # try({mod.k2.100.B = glm(temp$K2.B[temp$Dataset==100]~temp$LONG[temp$Dataset==100]); print(paste(spp,"K2.B 100 adjR2:",round(summary(mod.k2.100.B)$adj.r.squared,3)))},silent = T)
  # try({mod.k2.75.B = glm(temp$K2.B[temp$Dataset==75]~temp$LONG[temp$Dataset==75]); print(paste(spp,"K2.B 75 adjR2:",round(summary(mod.k2.75.B)$adj.r.squared,3)))},silent = T)
  # try({mod.k2.50.B = glm(temp$K2.B[temp$Dataset==50]~temp$LONG[temp$Dataset==50]); print(paste(spp,"K2.B 50 adjR2:",round(summary(mod.k2.50.B)$adj.r.squared,3)))},silent = T)
  # try({plot(temp$LONG,temp$K2.A, col=as.numeric(as.factor(temp$Dataset)),pch=as.numeric(as.factor(temp$Dataset))); abline(mod.k2.100.A,col="green"); abline(mod.k2.75.A,col="blue"); abline(mod.k2.50.A,col="red")})
  
  
  # try({mod.k3.100.A = lm(temp$K3.A[temp$Dataset==100]~temp$LONG[temp$Dataset==100]); print(paste(spp,"K3.A 100 adjR2:",round(summary(mod.k3.100.A)$adj.r.squared,3)))},silent = T)
  # try({mod.k3.75.A = lm(temp$K3.A[temp$Dataset==75]~temp$LONG[temp$Dataset==75]); print(paste(spp,"K3.A 75 adjR2:",round(summary(mod.k3.75.A)$adj.r.squared,3)))},silent = T)
  # try({mod.k3.50.A = lm(temp$K3.A[temp$Dataset==50]~temp$LONG[temp$Dataset==50]); print(paste(spp,"K3.A 50 adjR2:",round(summary(mod.k3.50.A)$adj.r.squared,3)))},silent = T)
  # try({mod.k3.100.B = lm(temp$K3.B[temp$Dataset==100]~temp$LONG[temp$Dataset==100]); print(paste(spp,"K3.B 100 adjR2:",round(summary(mod.k3.100.B)$adj.r.squared,3)))},silent = T)
  # try({mod.k3.75.B = lm(temp$K3.B[temp$Dataset==75]~temp$LONG[temp$Dataset==75]); print(paste(spp,"K3.B 75 adjR2:",round(summary(mod.k3.75.B)$adj.r.squared,3)))},silent = T)
  # try({mod.k3.50.B = lm(temp$K3.B[temp$Dataset==50]~temp$LONG[temp$Dataset==50]); print(paste(spp,"K3.B 50 adjR2:",round(summary(mod.k3.50.B)$adj.r.squared,3)))},silent = T)
  # try({mod.k3.100.C = lm(temp$K3.C[temp$Dataset==100]~temp$LONG[temp$Dataset==100]); print(paste(spp,"K3.C 100 adjR2:",round(summary(mod.k3.100.C)$adj.r.squared,3)))},silent = T)
  # try({mod.k3.75.C = lm(temp$K3.C[temp$Dataset==75]~temp$LONG[temp$Dataset==75]); print(paste(spp,"K3.C 75 adjR2:",round(summary(mod.k3.75.C)$adj.r.squared,3)))},silent = T)
  # try({mod.k3.50.C = lm(temp$K3.C[temp$Dataset==50]~temp$LONG[temp$Dataset==50]); print(paste(spp,"K3.C 50 adjR2:",round(summary(mod.k3.50.C)$adj.r.squared,3)))},silent = T)
  
  try({mod.p2.100.A = lm(temp$P2.A[temp$Dataset==100]~temp$LONG[temp$Dataset==100]); print(paste(spp,"P2.A 100 adjR2:",round(summary(mod.p2.100.A)$adj.r.squared,3)))},silent = T)
  try({mod.p2.75.A = lm(temp$P2.A[temp$Dataset==75]~temp$LONG[temp$Dataset==75]); print(paste(spp,"P2.A 75 adjR2:",round(summary(mod.p2.75.A)$adj.r.squared,3)))},silent = T)
  try({mod.p2.50.A = lm(temp$P2.A[temp$Dataset==50]~temp$LONG[temp$Dataset==50]); print(paste(spp,"P2.A 50 adjR2:",round(summary(mod.p2.50.A)$adj.r.squared,3)))},silent = T)
  try({mod.p2.100.B = lm(temp$P2.B[temp$Dataset==100]~temp$LONG[temp$Dataset==100]); print(paste(spp,"P2.B 100 adjR2:",round(summary(mod.p2.100.B)$adj.r.squared,3)))},silent = T)
  try({mod.p2.75.B = lm(temp$P2.B[temp$Dataset==75]~temp$LONG[temp$Dataset==75]); print(paste(spp,"P2.B 75 adjR2:",round(summary(mod.p2.75.B)$adj.r.squared,3)))},silent = T)
  try({mod.p2.50.B = lm(temp$P2.B[temp$Dataset==50]~temp$LONG[temp$Dataset==50]); print(paste(spp,"P2.B 50 adjR2:",round(summary(mod.p2.50.B)$adj.r.squared,3)))},silent = T)
  
  try({plot(temp$LONG,temp$P2.A, col=as.numeric(as.factor(temp$Dataset)),pch=as.numeric(as.factor(temp$Dataset))); abline(mod.p2.100.A,col="green"); abline(mod.p2.75.A,col="blue"); abline(mod.p2.50.A,col="red")})
  
  
  try({mod.p3.100.A = lm(temp$P3.A[temp$Dataset==100]~temp$LONG[temp$Dataset==100]); print(paste(spp,"P3.A 100 adjR2:",round(summary(mod.p3.100.A)$adj.r.squared,3)))},silent = T)
  try({mod.p3.75.A = lm(temp$P3.A[temp$Dataset==75]~temp$LONG[temp$Dataset==75]); print(paste(spp,"P3.A 75 adjR2:",round(summary(mod.p3.75.A)$adj.r.squared,3)))},silent = T)
  try({mod.p3.50.A = lm(temp$P3.A[temp$Dataset==50]~temp$LONG[temp$Dataset==50]); print(paste(spp,"P3.A 50 adjR2:",round(summary(mod.p3.50.A)$adj.r.squared,3)))},silent = T)
  try({mod.p3.100.B = lm(temp$P3.B[temp$Dataset==100]~temp$LONG[temp$Dataset==100]); print(paste(spp,"P3.B 100 adjR2:",round(summary(mod.p3.100.B)$adj.r.squared,3)))},silent = T)
  try({mod.p3.75.B = lm(temp$P3.B[temp$Dataset==75]~temp$LONG[temp$Dataset==75]); print(paste(spp,"P3.B 75 adjR2:",round(summary(mod.p3.75.B)$adj.r.squared,3)))},silent = T)
  try({mod.p3.50.B = lm(temp$P3.B[temp$Dataset==50]~temp$LONG[temp$Dataset==50]); print(paste(spp,"P3.B 50 adjR2:",round(summary(mod.p3.50.B)$adj.r.squared,3)))},silent = T)
  try({mod.p3.100.C = lm(temp$P3.C[temp$Dataset==100]~temp$LONG[temp$Dataset==100]); print(paste(spp,"P3.C 100 adjR2:",round(summary(mod.p3.100.C)$adj.r.squared,3)))},silent = T)
  try({mod.p3.75.C = lm(temp$P3.C[temp$Dataset==75]~temp$LONG[temp$Dataset==75]); print(paste(spp,"P3.C 75 adjR2:",round(summary(mod.p3.75.C)$adj.r.squared,3)))},silent = T)
  try({mod.p3.50.C = lm(temp$P3.C[temp$Dataset==50]~temp$LONG[temp$Dataset==50]); print(paste(spp,"P3.C 50 adjR2:",round(summary(mod.p3.50.C)$adj.r.squared,3)))},silent = T)
  
  par(mfrow=c(1,3))
  try({plot(temp$LONG,temp$P3.A, col=as.numeric(as.factor(temp$Dataset)),pch=as.numeric(as.factor(temp$Dataset))); abline(mod.p3.100.A,col="green"); abline(mod.p3.75.A,col="blue"); abline(mod.p3.50.A,col="red")})
  try({plot(temp$LONG,temp$P3.B, col=as.numeric(as.factor(temp$Dataset)),pch=as.numeric(as.factor(temp$Dataset))); abline(mod.p3.100.B,col="green"); abline(mod.p3.75.B,col="blue"); abline(mod.p3.50.B,col="red")})
  try({plot(temp$LONG,temp$P3.C, col=as.numeric(as.factor(temp$Dataset)),pch=as.numeric(as.factor(temp$Dataset))); abline(mod.p3.100.C,col="green"); abline(mod.p3.75.C,col="blue"); abline(mod.p3.50.C,col="red")})
  
  temp[is.na(temp)] = 0
  pdf(paste("~/",spp,"_NGSadmix_and_PCAngsd_assignments_missing_100-75-50.pdf",sep=""),
      height=4)
  par(mfrow=c(1,3),mar=c(6,2,0.5,0.5))
  k2.100 = as.matrix(temp[temp$Dataset==100,c("K2.A","K2.B")])
  k2.75 = as.matrix(temp[temp$Dataset==75,c("K2.A","K2.B")])
  k2.50 = as.matrix(temp[temp$Dataset==50,c("K2.A","K2.B")])
  if(nrow(k2.100)!=0) {barplot(t(k2.100),col=c("blue","green"),names=temp$LONG[temp$Dataset==100],las=2,main="K2 100")} else {plot(0,type="n",xaxt="n",yaxt="n",xlab="",main="K2 100")}
  if(nrow(k2.75)!=0) {barplot(t(k2.75),col=c("blue","green"),names=temp$LONG[temp$Dataset==75],las=2,main="K2 75")} else {plot(0,type="n",xaxt="n",yaxt="n",xlab="",main="K2 75")}
  if(nrow(k2.50)!=0) {barplot(t(k2.50),col=c("blue","green"),names=temp$LONG[temp$Dataset==50],las=2,main="K2 50")} else {plot(0,type="n",xaxt="n",yaxt="n",xlab="",main="K2 50")}
  
  k3.100 = as.matrix(temp[temp$Dataset==100,c("K3.A","K3.B","K3.C")])
  k3.75 = as.matrix(temp[temp$Dataset==75,c("K3.A","K3.B","K3.C")])
  k3.50 = as.matrix(temp[temp$Dataset==50,c("K3.A","K3.B","K3.C")])
  if(nrow(k3.100)!=0) {barplot(t(k3.100),col=c("blue","yellow","green"),names=temp$LONG[temp$Dataset==100],las=2,main="K3 100")} else {plot(0,type="n",xaxt="n",yaxt="n",xlab="",main="K3 100")}
  if(nrow(k3.75)!=0) {barplot(t(k3.75),col=c("blue","yellow","green"),names=temp$LONG[temp$Dataset==75],las=2,main="K3 75")} else {plot(0,type="n",xaxt="n",yaxt="n",xlab="",main="K3 75")}
  if(nrow(k3.50)!=0) {barplot(t(k3.50),col=c("blue","yellow","green"),names=temp$LONG[temp$Dataset==50],las=2,main="K3 50")} else {plot(0,type="n",xaxt="n",yaxt="n",xlab="",main="K3 50")}
  
  p2.100 = as.matrix(temp[temp$Dataset==100,c("P2.A","P2.B")])
  p2.75 = as.matrix(temp[temp$Dataset==75,c("P2.A","P2.B")])
  p2.50 = as.matrix(temp[temp$Dataset==50,c("P2.A","P2.B")])
  if(nrow(p2.100)!=0) {barplot(t(p2.100),col=c("blue","green"),names=temp$LONG[temp$Dataset==100],las=2,main="P2 100")} else {plot(0,type="n",xaxt="n",yaxt="n",xlab="",main="P2 100")}
  if(nrow(p2.75)!=0) {barplot(t(p2.75),col=c("blue","green"),names=temp$LONG[temp$Dataset==75],las=2,main="P2 75")} else {plot(0,type="n",xaxt="n",yaxt="n",xlab="",main="P2 75")}
  if(nrow(p2.50)!=0) {barplot(t(p2.50),col=c("blue","green"),names=temp$LONG[temp$Dataset==50],las=2,main="P2 50")} else {plot(0,type="n",xaxt="n",yaxt="n",xlab="",main="P2 50")}
  
  
  p3.100 = as.matrix(temp[temp$Dataset==100,c("P3.A","P3.B","P3.C")])
  p3.75 = as.matrix(temp[temp$Dataset==75,c("P3.A","P3.B","P3.C")])
  p3.50 = as.matrix(temp[temp$Dataset==50,c("P3.A","P3.B","P3.C")])
  if(nrow(p3.100)!=0) {barplot(t(p3.100),col=c("blue","yellow","green"),names=temp$LONG[temp$Dataset==100],las=2,main="P3 100")} else {plot(0,type="n",xaxt="n",yaxt="n",xlab="",main="P3 100")}
  if(nrow(p3.75)!=0) {barplot(t(p3.75),col=c("blue","yellow","green"),names=temp$LONG[temp$Dataset==75],las=2,main="P3 75")} else {plot(0,type="n",xaxt="n",yaxt="n",xlab="",main="P3 75")}
  if(nrow(p3.50)!=0) {barplot(t(p3.50),col=c("blue","yellow","green"),names=temp$LONG[temp$Dataset==50],las=2,main="P3 50")} else {plot(0,type="n",xaxt="n",yaxt="n",xlab="",main="P3 50")}
  
  dev.off()
  
  
  temp2 = df[df$Species==spp,]
  temp2$LONG = round(temp2$LONG,digits=0)
  temp2$LAT = round(temp2$LAT,digits=0)
  
  pdf(paste("~/",spp,"_NGSadmix_and_PCAngsd_assignments_missing_100-75-50_piecharts.pdf",sep=""),
      width=7,height=10)
  par(mfrow=c(3,1),mar=c(0,0,1,0))
  
  
  try({
    agg.k2 = aggregate(cbind(temp2$K2.A,temp2$K2.B)~temp2$Dataset+temp2$LONG+temp2$LAT,
                       FUN=function(x){sum(x,na.rm=T)})
    colnames(agg.k2) = c("Dataset","LONG","LAT","K2.A","K2.B")
    agg.k2$sum = rowSums(agg.k2[,c("K2.A","K2.B")])
    agg.k2= agg.k2[order(agg.k2$sum,decreasing = T),]
    
    
    agg.k2.100 = agg.k2[agg.k2$Dataset==100,]
    agg.k2.75 = agg.k2[agg.k2$Dataset==75,]
    agg.k2.50 = agg.k2[agg.k2$Dataset==50,]
  },silent = T) 
  
  
  raster::plot(wcdata,col=colorRampPalette(colors=c("black","grey"))(255),main="K2 100")
  
  try ({
    sapply(1:ncol(agg.k2.75),FUN=function(i){
      try({plotrix::floating.pie(
        xpos = agg.k2.100$LONG[i],
        ypos = agg.k2.100$LAT[i],
        x = c(as.numeric(agg.k2.100[i,c("K2.A","K2.B")])),
        radius = 0.3*sum(agg.k2.100[i,c("K2.A","K2.B")]),
        col = c("blue","green")
      )},silent = T)
    })
  },silent = T)
  
  
  raster::plot(wcdata,col=colorRampPalette(colors=c("black","grey"))(255),main="K2 75")
  
  try ({
    sapply(1:ncol(agg.k2.75),FUN=function(i){
      try({plotrix::floating.pie(
        xpos = agg.k2.75$LONG[i],
        ypos = agg.k2.75$LAT[i],
        x = c(as.numeric(agg.k2.75[i,c("K2.A","K2.B")])),
        radius = 0.3*sum(agg.k2.75[i,c("K2.A","K2.B")]),
        col = c("blue","green")
      )},silent = T)
    })
    
  },silent = T)
  
  raster::plot(wcdata,col=colorRampPalette(colors=c("black","grey"))(255),main="K2 50")
  
  try ({
    sapply(1:ncol(agg.k2.50),FUN=function(i){
      try({plotrix::floating.pie(
        xpos = agg.k2.50$LONG[i],
        ypos = agg.k2.50$LAT[i],
        x = c(as.numeric(agg.k2.50[i,c("K2.A","K2.B")])),
        radius = 0.3*sum(agg.k2.50[i,c("K2.A","K2.B")]),
        col = c("blue","green")
      )},silent = T)
    })
  },silent = T)
  
  
  try({
    agg.k3 = aggregate(cbind(temp2$K3.A,temp2$K3.B,temp2$K3.C)~temp2$Dataset+temp2$LONG+temp2$LAT,
                       FUN=function(x){sum(x,na.rm=T)})
    colnames(agg.k3) = c("Dataset","LONG","LAT","K3.A","K3.B","K3.C")
    agg.k3$sum = rowSums(agg.k3[,c("K3.A","K3.B","K3.C")])
    agg.k3= agg.k3[order(agg.k3$sum,decreasing = T),]
    
    agg.k3.100 = agg.k3[agg.k3$Dataset==100,]
    agg.k3.75 = agg.k3[agg.k3$Dataset==75,]
    agg.k3.50 = agg.k3[agg.k3$Dataset==50,]
    
  },silent = T)
  
  raster::plot(wcdata,col=colorRampPalette(colors=c("black","grey"))(255),main="K3 100")
  
  try ({
    sapply(1:ncol(agg.k3.75),FUN=function(i){
      try({plotrix::floating.pie(
        xpos = agg.k3.100$LONG[i],
        ypos = agg.k3.100$LAT[i],
        x = c(as.numeric(agg.k3.100[i,c("K3.A","K3.B","K3.C")])),
        radius = 0.3*sum(agg.k3.100[i,c("K3.A","K3.B","K3.C")]),
        col = c("blue","yellow","green")
      )},silent = T)
    })
    
  },silent = T)
  
  raster::plot(wcdata,col=colorRampPalette(colors=c("black","grey"))(255),main="K3 75")
  
  try ({
    sapply(1:ncol(agg.k3.75),FUN=function(i){
      try({plotrix::floating.pie(
        xpos = agg.k3.75$LONG[i],
        ypos = agg.k3.75$LAT[i],
        x = c(as.numeric(agg.k3.75[i,c("K3.A","K3.B","K3.C")])),
        radius = 0.3*sum(agg.k3.75[i,c("K3.A","K3.B","K3.C")]),
        col = c("blue","yellow","green")
      )},silent = T)
    })
    
  },silent = T)
  
  raster::plot(wcdata,col=colorRampPalette(colors=c("black","grey"))(255),main="K3 50")
  
  try ({
    sapply(1:ncol(agg.k3.50),FUN=function(i){
      try({plotrix::floating.pie(
        xpos = agg.k3.50$LONG[i],
        ypos = agg.k3.50$LAT[i],
        x = c(as.numeric(agg.k3.50[i,c("K3.A","K3.B","K3.C")])),
        radius = 0.3*sum(agg.k3.50[i,c("K3.A","K3.B","K3.C")]),
        col = c("blue","yellow","green")
      )},silent = T)
    })
  },silent = T)
  
  try({
    
    agg.p2 = aggregate(cbind(temp2$P2.A,temp2$P2.B)~temp2$Dataset+temp2$LONG+temp2$LAT,
                       FUN=function(x){sum(x,na.rm=T)})
    colnames(agg.p2) = c("Dataset","LONG","LAT","P2.A","P2.B")
    agg.p2$sum = rowSums(agg.p2[,c("P2.A","P2.B")])
    agg.p2= agg.p2[order(agg.p2$sum,decreasing = T),]
    
    agg.p2.100 = agg.p2[agg.p2$Dataset==100,]
    agg.p2.75 = agg.p2[agg.p2$Dataset==75,]
    agg.p2.50 = agg.p2[agg.p2$Dataset==50,]
  },silent = T)
  
  raster::plot(wcdata,col=colorRampPalette(colors=c("black","grey"))(255),main="P2 100")
  try ({
    sapply(1:ncol(agg.p2.75),FUN=function(i){
      try({plotrix::floating.pie(
        xpos = agg.p2.100$LONG[i],
        ypos = agg.p2.100$LAT[i],
        x = c(as.numeric(agg.p2.100[i,c("P2.A","P2.B")])),
        radius = 0.3*sum(agg.p2.100[i,c("P2.A","P2.B")]),
        col = c("blue","green")
      )},silent = T)
    })
  },silent = T)
  
  raster::plot(wcdata,col=colorRampPalette(colors=c("black","grey"))(255),main="P2 75")
  
  try ({
    sapply(1:ncol(agg.p2.75),FUN=function(i){
      try({plotrix::floating.pie(
        xpos = agg.p2.75$LONG[i],
        ypos = agg.p2.75$LAT[i],
        x = c(as.numeric(agg.p2.75[i,c("P2.A","P2.B")])),
        radius = 0.3*sum(agg.p2.75[i,c("P2.A","P2.B")]),
        col = c("blue","green")
      )},silent = T)
    })
  },silent = T)
  
  raster::plot(wcdata,col=colorRampPalette(colors=c("black","grey"))(255),main="P2 50")
  
  try ({
    sapply(1:ncol(agg.p2.50),FUN=function(i){
      try({plotrix::floating.pie(
        xpos = agg.p2.50$LONG[i],
        ypos = agg.p2.50$LAT[i],
        x = c(as.numeric(agg.p2.50[i,c("P2.A","P2.B")])),
        radius = 0.3*sum(agg.p2.50[i,c("P2.A","P2.B")]),
        col = c("blue","green")
      )},silent = T)
    })
  },silent = T)
  
  
  try({
    
    agg.p3 = aggregate(cbind(temp2$P3.A,temp2$P3.B,temp2$P3.C)~temp2$Dataset+temp2$LONG+temp2$LAT,
                       FUN=function(x){sum(x,na.rm=T)})
    colnames(agg.p3) = c("Dataset","LONG","LAT","P3.A","P3.B","P3.C")
    agg.p3$sum = rowSums(agg.p3[,c("P3.A","P3.B","P3.C")])
    agg.p3= agg.p3[order(agg.p3$sum,decreasing = T),]
    
    agg.p3.100 = agg.p3[agg.p3$Dataset==100,]
    agg.p3.75 = agg.p3[agg.p3$Dataset==75,]
    agg.p3.50 = agg.p3[agg.p3$Dataset==50,]
    
  },silent = T)
  
  raster::plot(wcdata,col=colorRampPalette(colors=c("black","grey"))(255),main="P3 100")
  
  try ({
    sapply(1:ncol(agg.p3.75),FUN=function(i){
      try({plotrix::floating.pie(
        xpos = agg.p3.100$LONG[i],
        ypos = agg.p3.100$LAT[i],
        x = c(as.numeric(agg.p3.100[i,c("P3.A","P3.B","P3.C")])),
        radius = 0.3*sum(agg.p3.100[i,c("P3.A","P3.B","P3.C")]),
        col = c("blue","yellow","green")
      )},silent = T)
    })
  },silent = T)
  raster::plot(wcdata,col=colorRampPalette(colors=c("black","grey"))(255),main="P3 75")
  try ({
    sapply(1:ncol(agg.p3.75),FUN=function(i){
      try({plotrix::floating.pie(
        xpos = agg.p3.75$LONG[i],
        ypos = agg.p3.75$LAT[i],
        x = c(as.numeric(agg.p3.75[i,c("P3.A","P3.B","P3.C")])),
        radius = 0.3*sum(agg.p3.75[i,c("P3.A","P3.B","P3.C")]),
        col = c("blue","yellow","green")
      )},silent = T)
    })
  },silent = T)
  raster::plot(wcdata,col=colorRampPalette(colors=c("black","grey"))(255),main="P3 50")
  try ({
    sapply(1:ncol(agg.p3.50),FUN=function(i){
      try({plotrix::floating.pie(
        xpos = agg.p3.50$LONG[i],
        ypos = agg.p3.50$LAT[i],
        x = c(as.numeric(agg.p3.50[i,c("P3.A","P3.B","P3.C")])),
        radius = 0.3*sum(agg.p3.50[i,c("P3.A","P3.B","P3.C")]),
        col = c("blue","yellow","green")
      )},silent = T)
    })
  },silent = T)
  dev.off()
  
}

library(hzar)

df$distance = df$LONG-min(df$LONG,na.rm=T)
df$nSamples = 1
run_hzar=function(qopt_ada,plot=T,index=1){
  
  if(plot==T){hzar.plot.obsData(qopt_ada)}
  
  qopt_model <-
    hzar::hzar.makeCline1DFreq(data=qopt_ada, scaling="fixed",tails="none");
  
  qopt_model <-
    hzar::hzar.model.addBoxReq(meta.model=qopt_model,
                               low=0,high=20)
  
  
  qopt_model_FitR <-
    hzar::hzar.first.fitRequest.old.ML(model=qopt_model ,
                                       obsData=qopt_ada,
                                       verbose=FALSE);
  
  qopt_model_FitR$mcmcParam$chainLength <- 1e5;
  qopt_model_FitR$mcmcParam$burnin <- 5e2;
  
  qopt_model_Fit <- hzar::hzar.doFit(qopt_model_FitR)
  if(plot==T){plot(hzar.mcmc.bindLL(qopt_model_Fit))}
  qopt_model_data <- hzar::hzar.dataGroup.add(qopt_model_Fit);
  
  qopt_model_data <-
    hzar::hzar.dataGroup.add(qopt_model_data,
                             hzar.chain.doSeq(hzar.next.fitRequest(qopt_model_Fit)));
  
  if(plot==T){
    hzar::hzar.plot.cline(qopt_model_data,xlim=c(0,20),ylim=c(0,1),col=index,
                          add=index!=1);
  }
  
  if(plot==T){hzar::hzar.plot.fzCline(qopt_model_data); ## won't work with just raw longitude
  }
  
  if(plot==T){print(hzar::hzar.getLLCutParam(qopt_model_data,c("center","width")))};
  
  qopt_model_null <- hzar::hzar.dataGroup.null(qopt_ada);
  qopt_null_ada <- list(clineModel = qopt_model_data,
                        nullModel  = qopt_model_null);
  
  qopt_null_ada_2 <- hzar::hzar.make.obsDataGroup(qopt_null_ada);
  qopt_null_ada_2 <- hzar::hzar.copyModelLabels(qopt_null_ada,qopt_null_ada_2);
  
  
  aicc=(hzar::hzar.AICc.hzar.obsDataGroup(qopt_null_ada_2))
  delta = aicc - min(aicc,na.rm=T)
  cline_delta = delta[1,]
  null_delta = delta[2,]
  
  ## second plot
  if (cline_delta > null_delta){
    ## print the null line
    hzar::hzar.plot.cline(qopt_model_null,xlim=c(0,20),ylim=c(0,1),col=index,
                          add=index!=1,lty=3,main="n.s.")
  } else if (null_delta - cline_delta > 10){
    ## print the line bold
    hzar::hzar.plot.cline(qopt_model_data,xlim=c(0,20),ylim=c(0,1),col=index,
                          add=index!=1,lwd=2,main=null_delta); ## i think this is the same as previous curve
  } else if (null_delta - cline_delta > 2){
    ## print the line regular
    hzar::hzar.plot.cline(qopt_model_data,xlim=c(0,20),ylim=c(0,1),col=index,
                          add=index!=1,main=null_delta); ## i think this is the same as previous curve
  } else {
    ## print the line dashed 
    hzar::hzar.plot.cline(qopt_model_data,xlim=c(0,20),ylim=c(0,1),col=index,
                          add=index!=1,lty=3,main=null_delta); ## i think this is the same as previous curve
    
    res<-list()
    res$mean_width<-mean(qopt_model_Fit$mcmcRaw[,2])
    # res$mean_width<-mean(qopt_model_Fit$mcmcRaw[,2])
    res$out<-qopt_model_Fit 
    res$dat<-qopt_model_data
    return(res)
    
  }
}

no_mcmc_hzar = function(qopt_ada,index=1){
  qopt_model <-    hzar.makeCline1DFreq(data=qopt_ada, scaling="fixed",tails="none");
  qopt_model <-    hzar.model.addBoxReq(meta.model=qopt_model,
                                        low=0,high=20);
  qopt_model_FitR <-
    hzar.first.fitRequest.old.ML(model=qopt_model ,
                                 obsData=qopt_ada,
                                 verbose=FALSE);
  qopt_model_FitR$mcmcParam$chainLength <- 1e5;
  qopt_model_FitR$mcmcParam$burnin <- 5e2;
  qopt_model_data2 <- hzar.dataGroup.add(qopt_model_FitR);
  qopt_model_data2 <-
    hzar.dataGroup.add(
      qopt_model_data2,
      hzar.chain.doSeq(hzar.next.fitRequest(qopt_model_FitR)));
  qopt_model_null <- hzar.dataGroup.null(qopt_ada);
  qopt_null_ada_A <- list(clineModel = qopt_model_data2,
                          nullModel  = qopt_model_null);
  qopt_null_ada_A_2 <- hzar.make.obsDataGroup(qopt_null_ada_A);
  qopt_null_ada_A_2 <- hzar.copyModelLabels(qopt_null_ada_A,qopt_null_ada_A_2);
  aicc=(hzar.AICc.hzar.obsDataGroup(qopt_null_ada_A_2))
  delta = aicc - min(aicc,na.rm=T)
  cline_delta = delta[1,]
  null_delta = delta[2,]
  if (cline_delta > null_delta){
    ## print the null line
    hzar.plot.cline(qopt_model_null,xlim=c(0,20),ylim=c(0,1),col=index,
                    add=index!=1)
  } else if (null_delta - cline_delta > 10){
    ## print the line bold
    hzar.plot.cline(qopt_model_data2,xlim=c(0,20),ylim=c(0,1),col=index,
                    add=index!=1,lwd=2); ## i think this is the same as previous curve
  } else if (null_delta - cline_delta > 2){
    ## print the line regular
    hzar.plot.cline(qopt_model_data2,xlim=c(0,20),ylim=c(0,1),col=index,
                    add=index!=1); ## i think this is the same as previous curve
  } else {
    ## print the line dashed 
    hzar.plot.cline(qopt_model_data2,xlim=c(0,20),ylim=c(0,1),col=index,
                    add=index!=1,lty=3); ## i think this is the same as previous curve
    
  }
  
}

for(spp in rev(sort(unique(df$Species)))){
  print(spp)
  pdf(paste("~/",spp,"_hzar_missing_100-75-50.pdf",sep=""),height=4)
  par(mfrow=c(1,3))
  for(dataset in unique(df$Dataset)){
    print(dataset)
    df = df[order(df$Ind),]
    add_qopt=df[df$Species==spp & df$Dataset==dataset,c("Ind","Species","LAT","LONG","distance","nSamples","P2.A","P2.B")]
    add_qopt = add_qopt[complete.cases(add_qopt),]
    
    if(nrow(add_qopt)>0){
      
      colnames(add_qopt)[7:ncol(add_qopt)]=paste("CLUSTER",seq(1:(ncol(add_qopt)-6)),sep="")
      qopt_ada = hzar::hzar.doMolecularData1DPops(add_qopt$distance,
                                                  add_qopt$CLUSTER1,
                                                  add_qopt$nSamples)
      is_positive=(cor(qopt_ada$frame[,1:2],use="pairwise.complete.obs")[2])>0
      if(is_positive==F){
        qopt_ada = hzar::hzar.doMolecularData1DPops(add_qopt$distance,
                                                    add_qopt$CLUSTER2,
                                                    add_qopt$nSamples)
      }
      
      res=run_hzar(qopt_ada,plot=F,index=1)
      
    } else {
      plot(0,type="n",xaxt="n",xlab="",ylab="",yaxt="n")
    }
  }
  dev.off()
}

speciesfiles = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/qopts_1/",
                          pattern="1.admix.Q.npy.csv",recursive = T,full.names = T)
speciesfiles = speciesfiles[!(grepl("SON",speciesfiles))]
speciesfiles = speciesfiles[!(grepl("CHI",speciesfiles))]
speciesfiles = speciesfiles[(grepl("Tgut",speciesfiles))]

df = df[order(df$Ind),]

sppdatasetout = NULL
for(spp in (sort(unique(df$Species)))){
  print(spp)
  this_speciesfiles = speciesfiles[(grepl(tolower(spp),speciesfiles))]
  
  for(dataset in c(75,50)){
    print(dataset)
    dataset_speciesfiles = this_speciesfiles[(grepl(paste(dataset,"bamlist",sep="\\."),this_speciesfiles))]
    
    meta_qopt=df[df$Species==spp & df$Dataset==dataset,c("Ind","Species","LAT","LONG","distance","nSamples")]
    
    
    for(file in dataset_speciesfiles){
      chr = basename(file)
      print(chr)
      chr=strsplit(chr,"Tgut_")[[1]][2]
      chr=strsplit(chr,"\\.")[[1]][1]
      try({
        individual_qopt = read.csv(file,header=F,sep="\t")
        colnames(individual_qopt)[1:ncol(individual_qopt)]=paste("CLUSTER",seq(1:(ncol(individual_qopt))),sep="")
        add_qopt = cbind(meta_qopt,individual_qopt)
        
        add_qopt = add_qopt[complete.cases(add_qopt),]
        
        if(nrow(add_qopt)>0){
          
          qopt_ada = hzar::hzar.doMolecularData1DPops(add_qopt$distance,
                                                      add_qopt$CLUSTER1,
                                                      add_qopt$nSamples)
          is_positive=(cor(qopt_ada$frame[,1:2],use="pairwise.complete.obs")[2])>0
          if(is_positive==F){
            qopt_ada = hzar::hzar.doMolecularData1DPops(add_qopt$distance,
                                                        add_qopt$CLUSTER2,
                                                        add_qopt$nSamples)
          }
          
          res=run_hzar(qopt_ada,plot=F,index=1)
          
          if(is.null(res)){
            mean_width=NA
            center = NA
          } else {
            mean_width=as.numeric(res$mean_width)
            center = as.numeric(mean(as.numeric(res$out$mcmcRaw[,1]),na.rm=T))
          }
          
          
          output = cbind(mean_width,center,spp,dataset,chr)
          
          #output=c(mean_width,center)
          
        } else {
          plot(0,type="n",xaxt="n",xlab="",ylab="",yaxt="n")
          
          output = cbind(mean_width=NA,center=NA,spp,dataset,chr)
        }
        
        if(is.null(sppdatasetout)){
          sppdatasetout = output
        } else {
          sppdatasetout = rbind(sppdatasetout,output)
          write.table(output,"~/sppdatasetout.temp",append = T)
        }
        
      })
      
      
    }
  }
}
print(sppdatasetout)
write.table(sppdatasetout,"~/sppdatasetout.temp",append = T)
