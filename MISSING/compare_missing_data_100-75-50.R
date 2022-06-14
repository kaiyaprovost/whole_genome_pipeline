
## FST
{
bil50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Amphispiza_bilineata_50_FST_slidingwindow.fst.FIXED.txt",data.table = F) 
bil75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Amphispiza_bilineata_75_FST_slidingwindow.fst.FIXED.txt",data.table = F) 

fla50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Auriparus_flaviceps_50_FST_slidingwindow.fst.FIXED.txt",data.table = F) 
fla75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Auriparus_flaviceps_75_FST_slidingwindow.fst.FIXED.txt",data.table = F) 

bru50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Campylorhynchus_brunneicapillus_50_FST_slidingwindow.fst.FIXED.txt",data.table = F)
bru75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Campylorhynchus_brunneicapillus_75_FST_slidingwindow.fst.FIXED.txt",data.table = F)

sin50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Cardinalis_sinuatus_50_FST_slidingwindow.fst.FIXED.txt",data.table = F) 
sin75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Cardinalis_sinuatus_75_FST_slidingwindow.fst.FIXED.txt",data.table = F) 

fus50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Melozone_fusca_50_FST_slidingwindow.fst.FIXED.txt",data.table = F) 
fus75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Melozone_fusca_75_FST_slidingwindow.fst.FIXED.txt",data.table = F) 

nit50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Phainopepla_nitens_50_FST_slidingwindow.fst.FIXED.txt",data.table = F) 
nit75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Phainopepla_nitens_75_FST_slidingwindow.fst.FIXED.txt",data.table = F) 

mel50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Polioptila_melanura_50_FST_slidingwindow.fst.FIXED.txt",data.table = F) 
mel75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Polioptila_melanura_75_FST_slidingwindow.fst.FIXED.txt",data.table = F) 

cri50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Toxostoma_crissale_50_FST_slidingwindow.fst.FIXED.txt",data.table = F) 
cri75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Toxostoma_crissale_75_FST_slidingwindow.fst.FIXED.txt",data.table = F) 

cur50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Toxostoma_curvirostre_50_FST_slidingwindow.fst.FIXED.txt",data.table = F) 
cur75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Toxostoma_curvirostre_75_FST_slidingwindow.fst.FIXED.txt",data.table = F) 

bel50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Vireo_bellii_50_FST_slidingwindow.fst.FIXED.txt",data.table = F) 
bel75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/SON_CHI_Vireo_bellii_75_FST_slidingwindow.fst.FIXED.txt",data.table = F) 


bil=merge(bil50,bil75,all=T)
fla=merge(fla50,fla75,all=T)
bel=merge(bel50,bel75,all=T)
#bel=bel75
#bel$Nsites50 = -1
#bel$Fst50 = -1
bru=merge(bru50,bru75,all=T)
cri=merge(cri50,cri75,all=T)
cur=merge(cur50,cur75,all=T)
fus=merge(fus50,fus75,all=T)
mel=merge(mel50,mel75,all=T)
nit=merge(nit50,nit75,all=T)
sin=merge(sin50,sin75,all=T)

full5070=rbind(bel,bil,bru,cri,cur,fla,fus,mel,nit,sin)

colnames(full5070)[colnames(full5070)=="spp"] = "species"

full=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/fst_COPY.temp",data.table = F)
full=full[,c("chr","midPos","Nsites","Fst","species")]

fullcombo = merge(full,full5070,all=T)

write.table(fullcombo,file="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/fst_100-75-50.temp",row.names = F)

corrs=cor(fullcombo[,c("Fst","Fst75","Fst50")],use = "pairwise.complete.obs")
print(corrs)

setwd("~")
for(spp in unique(fullcombo$species)){
  print(spp)
  subset = fullcombo[fullcombo$species==spp,]
  
  png(paste(spp,"_100-75-50.sites.corrs.png",sep=""),width=1500,height=500)
  par(mfrow=c(1,3))
  plot(subset$Nsites,subset$Nsites75,main="100-75")
  abline(b=1,a=0,col="red")
  plot(subset$Nsites,subset$Nsites50,main="100-50")
  abline(b=1,a=0,col="red")
  plot(subset$Nsites75,subset$Nsites50,main="75-50")
  abline(b=1,a=0,col="red")
  dev.off()
  
  png(paste(spp,"_100-75-50.FST.corrs.png",sep=""),width=1500,height=500)
  par(mfrow=c(1,3))
  plot(subset$Fst,subset$Fst75,main="100-75")
  abline(b=1,a=0,col="red")
  plot(subset$Fst,subset$Fst50,main="100-50")
  abline(b=1,a=0,col="red")
  plot(subset$Fst75,subset$Fst50,main="75-50")
  abline(b=1,a=0,col="red")
  dev.off()
}
}





## DXY
{
bel50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/bel_50_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)
bel75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/bel_75_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)

bil50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/bil_50_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)
bil75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/bil_75_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)

bru50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/bru_50_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)
bru75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/bru_75_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)

cri50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/cri_50_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)
cri75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/cri_75_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)

cur75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/cur_75_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)
cur50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/cur_50_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)

fla75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/fla_75_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)
#fla50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/fla_50_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)

fus75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/fus_75_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)
#fus50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/fus_50_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)

mel75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/mel_75_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)
mel50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/mel_50_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)

nit75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/nit_75_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)
nit50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/nit_50_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)

sin75=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/sin_75_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)
#sin50=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/sin_50_SON_Dxy_persite_FIXED.txt.windows.dxy.txt",data.table = F)



bil=merge(bil50,bil75,all=T)#; rm(bil50)#; rm(bil75)
bel=merge(bel50,bel75,all=T)#; rm(bel50)#; rm(bel75)
bru=merge(bru50,bru75,all=T)#; rm(bru50)#; rm(bru75)
cri=merge(cri50,cri75,all=T)#; rm(cri50)#; rm(cri75)
cur=merge(cur50,cur75,all=T)#; rm(cur50)#; rm(cur75)
# fus=merge(fus50,fus75,all=T)#; rm(fus50)#; rm(fus75)
mel=merge(mel50,mel75,all=T)#; rm(mel50)#; rm(mel75)
nit=merge(nit50,nit75,all=T)#; rm(nit50)#; rm(nit75)
# sin=merge(sin50,sin75,all=T)#; rm(sin50)#; rm(sin75)
# fla=merge(fla50,fla75,all=T)#; rm(fla50)#; rm(fla75)
fla=fla75#; rm(fla75)
fus=fus75#; rm(fus75)
sin=sin75#; rm(sin75)

full5070=gtools::smartbind(bel,bil,bru,cri,cur,fla,fus,mel,nit,sin)

full=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/dxy.temp",data.table = F)


fullcombo = merge(full,full5070,all=T,by=c("scafs","starts","species"))

write.table(fullcombo,file="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/dxy_100-75-50.temp",row.names = F)

fullcombo$dxymeans50[is.na(fullcombo$dxymeans50)] = -1

corrs=cor(fullcombo[,c("dxymeans","dxymeans75","dxymeans50")],use = "pairwise.complete.obs")
print(corrs)

for(spp in (unique(fullcombo$species))[10:1]){
  subset = fullcombo[fullcombo$species==spp,]
  
  png(paste(spp,"_100-75-50.DXY.corrs.png",sep=""),width=1500,height=500)
  par(mfrow=c(1,3))
  plot(subset$dxymeans,subset$dxymeans75,main="100-75")
  abline(b=1,a=0,col="red")
  plot(subset$dxymeans,subset$dxymeans50,main="100-50")
  abline(b=1,a=0,col="red")
  plot(subset$dxymeans75,subset$dxymeans50,main="75-50")
  abline(b=1,a=0,col="red")
  dev.off()
}
}


## TAJD BUT FIRST MAKE FULL  TAJ.TEMP FILE
{
listfiles = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas",pattern="thetasWindow.pestPG.FIXED.txt",
                       full.names = T)
listfiles = listfiles[!(grepl("NOWEIRD",listfiles))]

son=listfiles[grepl("SON",listfiles)]
chi=listfiles[grepl("CHI",listfiles)]
norm = listfiles[!(grepl("SON",listfiles))]
norm=norm[!(grepl("CHI",norm))]

dfs = lapply(son,FUN=function(x){
  data.table::fread(x,data.table = F)
})
dfc = lapply(chi,FUN=function(x){
  data.table::fread(x,data.table = F)
})
dfn = lapply(norm,FUN=function(x){
  data.table::fread(x,data.table = F)
})

dfss = do.call(gtools::smartbind,dfs)
dfcc = do.call(gtools::smartbind,dfc)
dfnn = do.call(gtools::smartbind,dfn)

merged = merge(dfss,dfcc,all=T)
merged = merge(merged,dfnn,all=T)

write.table(merged,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas/taj_tajpop.temp",row.names = F,quote = F)
}
## TAJD
{
  
  all50.CHI = data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ALL_50_CHI_thetasWindow.pestPG.FIXED.txt",data.table=F)
  all75.CHI = data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ALL_75_CHI_thetasWindow.pestPG.FIXED.txt",data.table=F)
  all75.SON = data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ALL_75_SON_thetasWindow.pestPG.FIXED.txt",data.table=F)
  
  #bel50.SON=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-SON/bel-50-SON-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  bil50.SON=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-SON/bil-50-SON-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  bru50.SON=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-SON/bru-50-SON-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  cri50.SON=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-SON/cri-50-SON-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  cur50.SON=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-SON/cur-50-SON-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  fla50.SON=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-SON/fla-50-SON-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  fus50.SON=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-SON/fus-50-SON-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  mel50.SON=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-SON/mel-50-SON-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  nit50.SON=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-SON/nit-50-SON-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  sin50.SON=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-SON/sin-50-SON-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)

  all50.SON=rbind(#bel50.SON,
                  bil50.SON,bru50.SON,cri50.SON,
                  cur50.SON,fla50.SON,fus50.SON,
                  mel50.SON,nit50.SON,sin50.SON)
  
  bel50.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-ALL/bel-50-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  bil50.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-ALL/bil-50-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  bru50.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-ALL/bru-50-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  cri50.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-ALL/cri-50-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  cur50.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-ALL/cur-50-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  #fla50.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-ALL/fla-50-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  fus50.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-ALL/fus-50-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  mel50.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-ALL/mel-50-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  nit50.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-ALL/nit-50-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  sin50.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/50-ALL/sin-50-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  
  all50.ALL=rbind(
    bel50.ALL,  bil50.ALL,bru50.ALL,cri50.ALL,cur50.ALL,
    #fla50.ALL,
    fus50.ALL,mel50.ALL,nit50.ALL,sin50.ALL)
  
  bel75.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/75-ALL/bel-75-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  bil75.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/75-ALL/bil-75-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  bru75.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/75-ALL/bru-75-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  cri75.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/75-ALL/cri-75-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  cur75.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/75-ALL/cur-75-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  #fla75.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/75-ALL/fla-75-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  #fus75.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/75-ALL/fus-75-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  mel75.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/75-ALL/mel-75-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  nit75.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/75-ALL/nit-75-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  sin75.ALL=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/tajima_1/75-ALL/sin-75-taj2.thetasWindow.pestPG.FIXED.txt",data.table= F)
  
  all75.ALL=rbind(
    bel75.ALL,  bil75.ALL,bru75.ALL,cri75.ALL,cur75.ALL,
    #fla75.ALL,fus75.ALL,
    mel75.ALL,nit75.ALL,sin75.ALL)
  
  
  all.CHI = merge(all50.CHI,all75.CHI,all=T)
  all.SON = merge(all50.SON,all75.SON,all=T)
  all.ALL = merge(all50.ALL,all75.ALL,all=T)
  all = merge(all.CHI,all.SON,all=T)
  full5070 = merge(all,all.ALL,all=T)
  
  write.table(full5070,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/taj_tajpop_75-50.temp")
  
  full=data.table::fread("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas/taj_tajpop.temp",data.table = F)
  
  
  fullcombo = merge(full,full5070,all=T,by=c("chr","midPos","species"))
  
  
  for(i in 1:ncol(fullcombo)){
    fullcombo[is.infinite(fullcombo[,i]),i] = NA
    if(i > 3){
      fullcombo[,i] = as.numeric(fullcombo[,i])
    }
  }
  
  write.table(fullcombo,file="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/taj_tajpop_100-75-50.temp",row.names = F)
  
  
  corrs=cor(fullcombo[,sort(colnames(fullcombo)[4:ncol(fullcombo)])[37:45]],use = "pairwise.complete.obs")
  print(corrs)
  corrplot::corrplot(corrs,method="color")
  
  #fullcombo[is.na(fullcombo)] = -99
  #fullcombo$Tajima.50.CHI[is.infinite(fullcombo$Tajima.50.CHI)] = -99
  
  for(spp in (unique(fullcombo$species))[1:10]){
    print(spp)
    subset = fullcombo[fullcombo$species==spp,sort(colnames(fullcombo)[4:ncol(fullcombo)])[37:45]]
    
    for(i in 1:ncol(subset)){
      subset[is.infinite(subset[,i]),i] = -99
    }
    
    subset[,which(colSums(!(is.na(subset)))==0)] = -99

    png(paste(spp,"_100-75-50.TAJ.corrs.png",sep=""),width=1500,height=1500)
    par(mfrow=c(3,3))
    plot(subset$Tajima,subset$Tajima.75,main="ALL 100-75")
    abline(b=1,a=0,col="red")
    plot(subset$Tajima,subset$Tajima.50,main="ALL 100-50")
    abline(b=1,a=0,col="red")
    plot(subset$Tajima.75,subset$Tajima.50,main="ALL 75-50")
    abline(b=1,a=0,col="red")
    
    plot(subset$Tajima.CHI,subset$Tajima.75.CHI,main="CHI 100-75")
    abline(b=1,a=0,col="red")
    plot(subset$Tajima.CHI,subset$Tajima.50.CHI ,main="CHI 100-50")
    abline(b=1,a=0,col="red")
    plot(subset$Tajima.50.CHI ,subset$Tajima.75.CHI,main="CHI 75-50")
    abline(b=1,a=0,col="red")
    
    plot(subset$Tajima.SON,subset$Tajima.75.SON,main="SON 100-75")
    abline(b=1,a=0,col="red")
    plot(subset$Tajima.SON,subset$Tajima.50.SON ,main="SON 100-50")
    abline(b=1,a=0,col="red")
    plot(subset$Tajima.50.SON ,subset$Tajima.75.SON,main="SON 75-50")
    abline(b=1,a=0,col="red")
    
    dev.off()
  }
}

