path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/FST/3.slidingWindowJobs/"
setwd(path)

species = c("Amphispiza_bilineata",
            "Auriparus_flaviceps-NOWEIRDMIN1C",
            "Campylorhynchus_brunneicapillus",
            "Cardinalis_sinuatus",
            "Melozone_fusca",
            "Phainopepla_nitens",
            "Polioptila_melanura",
            "Toxostoma_crissale",
            "Toxostoma_curvirostre",
            "Vireo_bellii-NOWEIRD"
)



## draw the window-by-window images
for (spp in species) {
  print(spp)
  file = paste(path,"SON_CHI_",spp,"_FST_slidingwindow_chrfix.fst",sep="")
  
  lines = readLines(file)
  end = substr(lines[1],start=nchar(lines[1])-2,stop=nchar(lines[1]))
  
  if (end != "Fst") {
    newline1 = paste(lines[1],"\tFst",sep="")
    lines[1] = newline1
    writeLines(lines,file)
  }
  
  windows = read.csv(file,sep="\t",row.names=NULL)
  head(windows)
  
  windows = (windows[order(-(windows$Nsites)),])
  
  #plot(windows$Fst,col=as.numeric(windows$chr))
  
  final = data.frame(Fst=character(0),Chrom=character(0),Freq=character(0))
  
  for (i in 1:length(unique(windows$chr))) {
    print(paste(spp, i,"of",length(unique(windows$chr))))
    chrom = unique(windows$chr)[[i]]
    temp = windows[windows$chr==chrom,]
    frequ = nrow(temp)
    fsts = temp$Fst
    
    eyes = rep(chrom,frequ)
    freqs = rep(frequ,frequ)
    
    toadd = data.frame(Fst=fsts,Chrom=eyes,Freq=freqs)
    #colnames(toadd) = c("Fst","Chrom","Freq")
    
    final = rbind(toadd,final)
    
  }
  
  final=final[order(-final$Freq,final$Chrom),] 
  final$Chrom = factor(final$Chrom,levels(final$Chrom)[unique(final$Chrom)])
  
  head(final)
  
  lvls = levels(final$Chrom)
  
  palette(c("grey","black"))
  
  png(paste(spp,"_windows_chrfix.png",sep=""),
      width=700,height=350)
  par(bg=NA,col.axis="white",fg="white",
      col.lab="white",col.main="white")
  plot(final$Fst,col=as.numeric(as.factor(final$Chrom)),cex=0.2,
       main=spp,xlab="Window (Scaffold)",ylab="FST",
       ylim=c(0,1))
  dev.off()
  
  palette(c("red","orange","goldenrod","green","blue","purple",
            "cyan","grey","brown","magenta"))
  
  png(paste(spp,"_windows_chrfix_zoomed.png",sep=""),
      width=700,height=350)
  par(bg=NA,col.axis="white",fg="white",
      col.lab="white",col.main="white")
  plot(final$Fst,col=as.numeric(as.factor(final$Chrom)),cex=0.2,
       main=spp,xlab="Window (Scaffold)",ylab="FST")
  dev.off()
  
  
  shortchr = substr(final$Chrom ,1,8)
  #only_chr = final[which(shortchr %in% "PseudoNC"),]
  only_chr = final
  
#  chrnames = substr(unique(as.factor(only_chr$Chrom)),24,50)
  chrnames = (unique(as.factor(only_chr$Chrom)))
  
  
  png(paste(spp,"_windows_chrfix_chrs.png",sep=""),
      width=700,height=350)
  plot(only_chr$Fst,col=as.numeric(as.factor(only_chr$Chrom)),cex=0.2,
       main=spp,xlab="Window (Scaffold)",ylab="FST",
       ylim=c(0,1))
  text(y=0.8,x=seq(1,100000,round(100000/length(chrnames))),labels=chrnames,cex=0.5)
  
  dev.off()
  
  
  png(paste(spp,"_windows_chrfix_chrs_zoon.png",sep=""),
      width=700,height=350)
  plot(only_chr$Fst,col=as.numeric(as.factor(only_chr$Chrom)),cex=0.2,
       main=spp,xlab="Window (Scaffold)",ylab="FST")
  text(y=mean(only_chr$Fst),x=seq(1,100000,round(100000/length(chrnames))),labels=chrnames,cex=0.5)
  
  
  dev.off()
  
}


## calculate for chromosomes


outputdf = c()
for (spp in species) {
  print(spp)
  file = paste(path,"SON_CHI_",spp,"_FST_slidingwindow_chrfix.fst",sep="")
  
  
  windows = read.csv(file,sep="\t",row.names=NULL)
  
  quantiles = quantile(windows$Fst,probs=c(0.05,0.95))
  mean = mean(windows$Fst)
  sd = sd(windows$Fst)
  sdbounds = c(mean-5*sd,mean+5*sd)
  plotmin = max(min(c(windows$Fst,quantiles,sdbounds),na.rm=T),0,na.rm=T)
  plotmax = min(max(c(windows$Fst,quantiles,sdbounds),na.rm=T),1,na.rm=T)
  
  windows$color = "black"
  windows$color[windows$Fst>quantiles[2]] = "darkred"
  windows$color[windows$Fst<quantiles[1]] = "blue"
  windows$color[windows$Fst>sdbounds[2]] = "red"
  windows$color[windows$Fst<sdbounds[1]] = "cyan"
  
  windows$chrcolor = "black"
  numsites = nrow(windows)
  
  ## get the number in each category
  num_upper_quan = sum(windows$color=="darkred")
  num_lower_quan = sum(windows$color=="blue")
  num_upper_sdv5 = sum(windows$color=="red")
  num_lower_sdv5 = sum(windows$color=="cyan")
  
  chr="ALL"
  
  inds = which(windows$color!="black")
  cols = windows$color[inds]    
  uniquepeaks=c(cols[1])
  
  if(length(inds)>1){
  
  for(i in 2:length(inds)){
    dif = inds[i] - inds[i-1]
    if(dif != 1) {newbreak=T} else { if(cols[i] != cols[i-1]) {newbreak=T} else {newbreak=F}}
    if(newbreak==T){
      uniquepeaks = c(uniquepeaks,cols[i])
    }
  }
  }
  
  pks_upper_quan = sum(uniquepeaks=="darkred")
  pks_lower_quan = sum(uniquepeaks=="blue")
  pks_upper_sdv5 = sum(uniquepeaks=="red")
  pks_lower_sdv5 = sum(uniquepeaks=="cyan")
  
  todf = cbind(spp,chr,numsites,
               num_upper_quan,num_lower_quan,num_upper_sdv5,num_lower_sdv5,
               pks_upper_quan,pks_lower_quan,pks_upper_sdv5,pks_lower_sdv5)
  
  if(is.null(outputdf)) {
    outputdf = todf
  } else {
    outputdf = rbind(outputdf,todf)
  }
  
  png(paste(spp,"_windows_chrfix-chrALL.png",sep=""),
      width=700,height=700)
  par(mfrow=c(2,1))
  plot(windows$Fst,
       ylim=c(plotmin,plotmax),
       col=windows$chr,
       main=paste(spp,"ALL"))
  abline(h=quantiles,col=c("blue","darkred"),lwd=2)
  abline(h=sdbounds,col=c("cyan","red"),lwd=2)
  plot(windows$Fst,
       ylim=c(plotmin,plotmax),
       col=windows$color,
       main=paste(spp,"ALL"))
  abline(h=quantiles,col=c("blue","darkred"),lwd=2)
  abline(h=sdbounds,col=c("cyan","red"),lwd=2)
  dev.off()
  
  
  for(chr in unique(windows$chr)) {
    print(chr)
    temp = windows[windows$chr==chr,]
    
    quantilesC = quantile(temp$Fst,probs=c(0.05,0.95))
    meanC = mean(temp$Fst)
    sdC = sd(temp$Fst)
    sdboundsC = c(meanC-5*sdC,meanC+5*sdC)
    plotminC = max(min(c(temp$Fst,quantilesC,sdboundsC),na.rm=T),0,na.rm=T)
    plotmaxC = min(max(c(temp$Fst,quantilesC,sdboundsC),na.rm=T),1,na.rm=T)
    
    
    temp$chrcolor[temp$Fst>quantilesC[2]] = "goldenrod"
    temp$chrcolor[temp$Fst<quantilesC[1]] = "pink"
    temp$chrcolor[temp$Fst>sdboundsC[2]] = "darkgreen"
    temp$chrcolor[temp$Fst<sdboundsC[1]] = "green"
    
    
    numsitesC = nrow(temp)
    
    ## get the number in each category
    num_upper_quanC = sum(temp$chrcolor=="goldenrod")
    num_lower_quanC = sum(temp$chrcolor=="pink")
    num_upper_sdv5C = sum(temp$chrcolor=="darkgreen")
    num_lower_sdv5C = sum(temp$chrcolor=="green")

    inds = which(temp$chrcolor!="black")
    cols = temp$chrcolor[inds]    
    uniquepeaks=c(cols[1])
    if(length(inds)>1){
    for(i in 2:length(inds)){
      dif = inds[i] - inds[i-1]
      if(dif != 1) {newbreak=T} else { if(cols[i] != cols[i-1]) {newbreak=T} else {newbreak=F}}
      if(newbreak==T){
        uniquepeaks = c(uniquepeaks,cols[i])
      }
    }
    }
    
    pks_upper_quanC = sum(uniquepeaks=="goldenrod")
    pks_lower_quanC = sum(uniquepeaks=="pink")
    pks_upper_sdv5C = sum(uniquepeaks=="darkgreen")
    pks_lower_sdv5C = sum(uniquepeaks=="green")
    
    todf = cbind(spp,paste(chr,"NORM",sep="-"),numsitesC,
                 num_upper_quanC,num_lower_quanC,num_upper_sdv5C,num_lower_sdv5C,
                 pks_upper_quanC,pks_lower_quanC,pks_upper_sdv5C,pks_lower_sdv5C)
    
    if(is.null(outputdf)) {
      outputdf = todf
    } else {
      outputdf = rbind(outputdf,todf)
    }
    
    num_upper_quan = sum(temp$color=="darkred")
    num_lower_quan = sum(temp$color=="blue")
    num_upper_sdv5 = sum(temp$color=="red")
    num_lower_sdv5 = sum(temp$color=="cyan")
    
    inds = which(temp$color!="black")
    cols = temp$color[inds]    
    uniquepeaks=c(cols[1])
    if(length(inds)>1){
    for(i in 2:length(inds)){
      dif = inds[i] - inds[i-1]
      if(dif != 1) {newbreak=T} else { if(cols[i] != cols[i-1]) {newbreak=T} else {newbreak=F}}
      if(newbreak==T){
        uniquepeaks = c(uniquepeaks,cols[i])
      }
    }
    }
    
    pks_upper_quan = sum(uniquepeaks=="darkred")
    pks_lower_quan = sum(uniquepeaks=="blue")
    pks_upper_sdv5 = sum(uniquepeaks=="red")
    pks_lower_sdv5 = sum(uniquepeaks=="cyan")
    
    todf = cbind(spp,chr,numsitesC,
                 num_upper_quan,num_lower_quan,num_upper_sdv5,num_lower_sdv5,
                 pks_upper_quan,pks_lower_quan,pks_upper_sdv5,pks_lower_sdv5)
    
    if(is.null(outputdf)) {
      outputdf = todf
    } else {
      outputdf = rbind(outputdf,todf)
    }
    
    png(paste(spp,"_windows_chrfix-chr",chr,".png",sep=""),
        width=700,height=700)
    par(mfrow=c(2,1))
    
    plot(temp$midPos,temp$Fst,
         ylim=c(plotmin,plotmax),
         col=temp$color,
         main=paste(spp,chr))
    abline(h=quantiles,col=c("blue","darkred"),lwd=2)
    abline(h=sdbounds,col=c("cyan","red"),lwd=2)
    abline(h=quantilesC,col=c("pink","goldenrod"),lwd=1,lty=3)
    abline(h=sdboundsC,col=c("green","darkgreen"),lwd=1,lty=3)
    
    plot(temp$midPos,temp$Fst,
         ylim=c(plotminC,plotmaxC),
         col=temp$chrcolor,
         main=paste(spp,chr))
    abline(h=quantiles,col=c("blue","darkred"),lwd=1,lty=3)
    abline(h=sdbounds,col=c("cyan","red"),lwd=1,lty=3)
    abline(h=quantilesC,col=c("pink","goldenrod"),lwd=2)
    abline(h=sdboundsC,col=c("green","darkgreen"),lwd=2)
    
    dev.off()
  }
}
colnames(outputdf) = c("Species","Chromosomes","NumWindows",
                       "Above95%","Below5%","Above5SD","Below5SD",
                       "PeaksAbove95%","PeaksBelow5%","PeaksAbove5SD","PeaksBelow5SD")
write.table(outputdf,"fst_quantiles.txt",row.names = F,sep="\t",quote=F)

#####

# testfile = "/Users/kprovost/Documents/Dissertation/CHAPTER2_GENOMES/ANALYSIS/FST/3.slidingWindowJobs/SON_CHI_Toxostoma_curvirostre_FST_slidingwindow.fst"
# test = read.csv(testfile,sep="\t")
# 
# summary(test$Fst)
# avg = mean(test$Fst)
# stdev = sd(test$Fst)
# test$zscore = (test$Fst - avg) / stdev
# quan1 = quantile(test$Fst,0.05)
# quan2 = quantile(test$Fst,0.95)
# test$quantile = 1
# test$quantile[test$Fst <= quan1] = 2
# test$quantile[test$Fst >= quan2] = 3
# 
# #hist(test$zscore)
# summary(test$zscore)
# test$abszscore = abs(test$zscore)
# #plot(test$Fst,test$zscore)
# 
# #plot(test$Fst,col="black")
# plot(test$Fst,col=as.numeric(test$quantile))

