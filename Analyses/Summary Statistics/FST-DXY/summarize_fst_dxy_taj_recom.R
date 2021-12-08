chromsdf = NULL

files=list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/FST/3.slidingWindowJobs/",
                 pattern="fst$",full.names = T)
for(fstfile in files){
  #fstfile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/FST/3.slidingWindowJobs/SON_CHI_Amphispiza_bilineata_FST_slidingwindow_chrfix.fst"
  spp=substr(strsplit(basename(fstfile),"_")[[1]][4],1,3)
  cat(spp,sep = "\n")
  
  
  df = read.table(fstfile,header=T,stringsAsFactors = F,strip.white = T,sep="\t")
  
  chroms=unique(df$chr)
  addchroms = cbind(spp,"fst",chroms)
  
  if(is.null(chromsdf)){
    chromsdf=addchroms
  } else {
    chromsdf=rbind(chromsdf,addchroms)
  }
  
  wholemean=mean(df$Fst,na.rm=T)
  wholesd=sd(df$Fst,na.rm=T)
  numwindows = nrow(df)
  #cat(paste(signif(wholemean,3),"±",signif(wholesd,3)," (",numwindows,")",sep=""),sep="\n")
  
  aggmean = aggregate(df$Fst~df$chr,FUN=function(x){mean(x,na.rm=T)})
  
  print(aggmean)
  
  #aggsd = aggregate(df$Fst~df$chr,FUN=function(x){sd(x,na.rm=T)})
  sppmean = mean(aggmean$`df$Fst`,na.rm=T)
  sppsd = sd(aggmean$`df$Fst`,na.rm=T)
  #cat(paste(signif(sppmean,3),"±",signif(sppsd,3)," (",nrow(aggmean),")",sep=""),sep="\n")
}

dxyfiles = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY",
                      pattern="_4_SON_Dxy_WINDOWS_chrfix_1-ALL.txt$",full.names = T,recursive = T)
for(dxyfile in dxyfiles){
  #fstfile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/FST/3.slidingWindowJobs/SON_CHI_Amphispiza_bilineata_FST_slidingwindow_chrfix.fst"
  df = read.table(dxyfile,header=T,strip.white = T,sep=" ",fill = T,stringsAsFactors = F)
  
  spp=strsplit(basename(dxyfile),"_")[[1]][1]
  cat(spp,sep = "\n")
  
  chroms=unique(df$scafs)
  addchroms = cbind(spp,"dxy",chroms)
  
  if(is.null(chromsdf)){
    chromsdf=addchroms
  } else {
    chromsdf=rbind(chromsdf,addchroms)
  }
  
  wholemean=mean(df$means,na.rm=T)
  wholesd=sd(df$means,na.rm=T)
  numwindows = nrow(df)
  #cat(paste(signif(wholemean,3),"±",signif(wholesd,3)," (",numwindows,")",sep=""),sep="\n")
  
  
  aggmean = aggregate(df$means~df$scafs,FUN=function(x){mean(x,na.rm=T)})
  #aggsd = aggregate(df$Fst~df$chr,FUN=function(x){sd(x,na.rm=T)})
  sppmean = mean(aggmean$`df$means`,na.rm=T)
  sppsd = sd(aggmean$`df$means`,na.rm=T)
  print(aggmean)
  #cat(paste(signif(sppmean,3),"±",signif(sppsd,3)," (",nrow(aggmean),")",sep=""),sep="\n")
}

duprecom = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigrecom_duplicates.may2020.txt"
df = read.table(duprecom,header=T,fill = T,stringsAsFactors = F)
dim(df)
for(spp in sort(unique(df$species))){
  cat(spp,sep="\n")
  small=df[df$species==spp,]
  
  wholemean=mean(small$weighted_recomb,na.rm=T)
  wholesd=sd(small$weighted_recomb,na.rm=T)
  numwindows = nrow(small)
  cat(paste(signif(wholemean,3),"±",signif(wholesd,3)," (",numwindows,")",sep=""),sep="\n")
  
  aggmean = aggregate(small$weighted_recomb~small$chr,FUN=function(x){mean(x,na.rm=T)})
  #aggsd = aggregate(df$Fst~df$chr,FUN=function(x){sd(x,na.rm=T)})
  sppmean = mean(aggmean$`small$weighted_recomb`,na.rm=T)
  sppsd = sd(aggmean$`small$weighted_recomb`,na.rm=T)
  cat(paste(signif(sppmean,3),"±",signif(sppsd,3)," (",nrow(aggmean),")",sep=""),sep="\n")
  
  
}

singrecom = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigrecom.may2020.txt"
df = read.table(singrecom,header=T,fill = T,stringsAsFactors = F)
dim(df)
for(spp in sort(unique(df$species))){
  cat(spp,sep="\n")
  small=df[df$species==spp,]
  
  chroms=unique(df$chr)
  addchroms = cbind(spp,"rec",chroms)
  
  if(is.null(chromsdf)){
    chromsdf=addchroms
  } else {
    chromsdf=rbind(chromsdf,addchroms)
  }
  
  wholemean=mean(small$weighted_recomb,na.rm=T)
  wholesd=sd(small$weighted_recomb,na.rm=T)
  numwindows = nrow(small)
  #cat(paste(signif(wholemean,3),"±",signif(wholesd,3)," (",numwindows,")",sep=""),sep="\n")
  
  aggmean = aggregate(small$weighted_recomb~small$chr,FUN=function(x){mean(x,na.rm=T)})
  #aggsd = aggregate(df$Fst~df$chr,FUN=function(x){sd(x,na.rm=T)})
  sppmean = mean(aggmean$`small$weighted_recomb`,na.rm=T)
  sppsd = sd(aggmean$`small$weighted_recomb`,na.rm=T)
  #cat(paste(signif(sppmean,3),"±",signif(sppsd,3)," (",nrow(aggmean),")",sep=""),sep="\n")
  print(aggmean)
  
}

chromsdf = as.data.frame(unique(chromsdf))
colnames(chromsdf) =c("spp","data","chr")
allunique = unique(chromsdf[,c("spp","chr")])
table(chromsdf[,c(3,2)])
table(allunique$spp)




listfiles = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas/",
             pattern="-taj2.thetas.idx_CHRFIX.pestPG",recursive=T,full.names = T)
listfiles = listfiles[!(grepl("CHI",listfiles))]
listfiles = listfiles[!(grepl("SON",listfiles))]

for(file in listfiles){
  
}
