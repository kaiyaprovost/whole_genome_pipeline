windowsRecalculator=function(data,recomb_file,window_size = 100000,overlap = 20000,reset_overlap = T,overwrite=F,chrom=NULL){
  
  verbose=F
  
  if(overlap > window_size){
    stop("ERROR: Overlap cannot be greater than window size!")
    
  }
  
  start=min(data$start)
  stop=max(data$end)
  
  original_window_size = data[1,"end"]-data[1,"start"]
  
  if(overlap < original_window_size) {
    print("WARNING: Original window size is bigger than overlap!")
    verbose=T
    if(reset_overlap==T){
      print("WARNING: Resizing overlap to match original window size.")
      print(paste("Old:",overlap,"New:",original_window_size))
      overlap = original_window_size
    }
    
  }
  
  #print(overlap)
  
  if(is.null(chrom)) {
    outfilename_noext = paste(recomb_file,"_w",window_size,"_o",overlap,sep="")
  } else {
    outfilename_noext = paste(recomb_file,"_w",window_size,"_o",overlap,"_chrom",chrom,sep="")
  }
  
  if (file.exists(paste(outfilename_noext,".txt",sep="")) && overwrite==F){
    print("SKIPPING")
  } else {
    
    windowstarts=seq(start,stop,overlap)
    windowstops=windowstarts+window_size
    
    ## to do this: get all windows that overlap the range
    ## calculate how much overlap they all have
    ## do a weighted average
    
    new_data=as.data.frame(cbind(windowstarts,windowstops))
    new_data$weighted_recomb = NA
    
    #for (i in 1:10) {
    for (i in 1:length(windowstarts)) {
      if (i %% 1000 == 0) {paste(print(i),"/",length(windowstarts))}
      
      thiswindow_start = windowstarts[i]
      thiswindow_stop = windowstops[i]
      
      subset = data[data$end <= thiswindow_stop+window_size+1,]
      subset = subset[subset$start >= thiswindow_start-window_size-1,]
      
      if(verbose==T && i == 1){
        print(head(subset))
        print(thiswindow_start)
        print(thiswindow_stop)
        print("#####")
      }
      
      new_data$weighted_recomb[i] = calculateOverlap(subset,thiswindow_start,thiswindow_stop,verbose)
      verbose=F
      
    }
    
    write.csv(new_data,paste(outfilename_noext,".txt",sep=""),row.names = F)
    
    plotmax = max(data$recombRate,new_data$weighted_recomb,na.rm=T)
    
    png(paste(outfilename_noext,".png",sep=""))
    par(mfrow=c(2,1))
    plot(data$recombRate,cex=0.5,ylim=c(0,plotmax))
    plot(new_data$weighted_recomb,cex=0.5,ylim=c(0,plotmax))
    dev.off()
    
    return(new_data)
  }
}

calculateOverlap = function(subset,thiswindow_start,thiswindow_stop,verbose=F){
  ## get the ones that overlap
  window_starts_between=dplyr::between(subset$start,thiswindow_start,thiswindow_stop)
  window_stops_between=dplyr::between(subset$end,thiswindow_start,thiswindow_stop)
  keep=rownames(subset)[which(window_starts_between | window_stops_between)]
  
  if(length(keep)<=0) {
    return(NA)
  } else {
    
    entire=rownames(subset)[which(window_starts_between & window_stops_between)]
    not_entire=intersect(rownames(subset)[which(!(window_starts_between & window_stops_between))],keep)
    ## calculate how much of each window is between
    overlapping=subset[rownames(subset)  %in% keep,]
    #fulloverlap=subset[entire,]
    #partialoverlap=subset[not_entire,]
    overlapping$overlap = NA
    
    if(length(entire)>0) {
      overlapping[rownames(overlapping) %in% entire,]$overlap=1
    } 
    for (this_value in not_entire){
      this_index = (rownames(overlapping) %in% this_value)
      thisline =  overlapping[this_index,]
      start_overlap = dplyr::between(thiswindow_start,thisline$start,thisline$end)
      stop_overlap = dplyr::between(thiswindow_stop,thisline$start,thisline$end)
      ## check if the window start/stops are between start and ends for ones in notentire
      
      if(start_overlap && !(stop_overlap)){
        amount_over=(thisline$end-thiswindow_start)/window_size
        #print(1)
      } else if (!(start_overlap) && stop_overlap) {
        amount_over=(thiswindow_stop-thisline$start)/window_size
        #print(-1)
      } else if (start_overlap && stop_overlap) {
        amount_over=NA
      } else if (!(start_overlap) && !(stop_overlap)) {
        amount_over=NaN
      }
      overlapping$overlap[this_index] = amount_over
    }
    
    if(verbose==T){
      print(head(overlapping))
    }
    
    overlapping$weighted = (overlapping$recombRate * overlapping$overlap)
    weighted_average = sum(overlapping$weighted)/sum(overlapping$overlap)
    return(weighted_average)
  }
}

## THIS WONT CALCULATE CORRECTLY IF THE WINDOW SIZE IS SMALLER THAN THE WINDOWS ESTIMATED

options(scipen=10)

filelist = sort(c(#list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION/",
                  #           pattern="PREDICT.BSCORRECTED.txt$",full.names = T,recursive = T),
                  list.files(path="/Users/kprovost/Documents/TempDropbox/BS/",
                             pattern="BOTH_cri.PREDICT.BSCORRECTED.txt$",full.names = T,recursive = F)))
#filelist = filelist[grepl("geno",filelist)]
#filelist = filelist[!(grepl("sorted.sorted",filelist))]


for(j in 1:length(filelist)) {
  recomb_file = filelist[j]
  print(paste(j,"/",length(filelist)))
  print(recomb_file)
  
  #recomb_file = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION/Amphispiza-bilineata-called.geno.PseudoNC_011462.1_Tgut_1.fixedchroms.converted.PREDICT.BSCORRECTED.txt"
  #recomb_file="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION//Amphispiza-bilineata-called.geno.PseudoNC_011469.1_Tgut_5.fixedchroms.converted.PREDICT.BSCORRECTED.txt"
  
  data=read.table(recomb_file,header=T)
  #plot(data$recombRate,
  #     cex=0.5)
  #plot(data[,c(4,5)])
  
  if(nrow(data) <= 0){
    print("NO DATA IN FILE -- SKIPPING")
  } else {
    
    ## force window size via a weighted average
    window_size = 100000
    overlap = 10000
    
    if(length(unique(data$chrom))==1){
      new_data = windowsRecalculator(data,recomb_file,window_size,overlap,reset_overlap=F,overwrite=F)
      new_data2 = windowsRecalculator(data,recomb_file,window_size,overlap,reset_overlap=T,overwrite=F)
    } else {
      
      combined_data = NULL
      combined_data2 = NULL
      for (chrom in unique(data$chrom)){
        subset = data[data$chrom==chrom,]
        new_data = windowsRecalculator(subset,recomb_file,window_size,overlap,reset_overlap=F,overwrite=T,chrom)
        new_data2 = windowsRecalculator(subset,recomb_file,window_size,overlap,reset_overlap=T,overwrite=T,chrom)
        
        new_data$chrom = rep(chrom,nrow(new_data))
        new_data2$chrom = rep(chrom,nrow(new_data2))
        
        if(is.null(combined_data)) {
          combined_data = new_data
          combined_data2 = new_data2
        } else {
          combined_data = rbind(combined_data,new_data)
          combined_data2 = rbind(combined_data2,new_data2)
        }
        
        
      }
      
      outfilename_noext = paste(recomb_file,"_w",window_size,"_o",overlap,"_genome",sep="")
      outfilename_noext2 = paste(recomb_file,"_w",window_size,"_o",overlap,"_genome_changeoverlaps",sep="")
      
      write.csv(combined_data,paste(outfilename_noext,".txt",sep=""),row.names = F)
      write.csv(combined_data2,paste(outfilename_noext2,".txt",sep=""),row.names = F)
      
      plotmax = max(data$recombRate,combined_data$weighted_recomb,
                    combined_data2$weighted_recomb,na.rm=T)
      png(paste(outfilename_noext,".png",sep=""))
      par(mfrow=c(3,1))
      plot(data$recombRate,cex=0.5,ylim=c(0,plotmax),
           col=as.numeric(as.factor(data$chrom)))
      plot(combined_data$weighted_recomb,cex=0.5,ylim=c(0,plotmax),
           col=as.numeric(as.factor(combined_data$chrom)))
      plot(combined_data2$weighted_recomb,cex=0.5,ylim=c(0,plotmax),
           col=as.numeric(as.factor(combined_data2$chrom)))
      dev.off()
      
    }
    
    ## this isn't working if the window of the file is bigger than our window
  }
  
}

## merge recomb with other data -- going to need to do this
bigplot = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.fst.dxy.taj.may2020.txt",
                     sep=",",header=T,stringsAsFactors = F)

#smallplot = bigplot[bigplot$species=="bil",]


all_recom_files=list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION/",
                           pattern="_genome.*txt",
                           recursive=T,full.names = T)
all_recom_files = all_recom_files[!(grepl("BSCORR",all_recom_files))]
all_recom_files = all_recom_files[!(grepl("sorted.sorted",all_recom_files))]
all_recom_files = all_recom_files[(grepl("o10000_genome_changeoverlaps.txt",all_recom_files))]


{
df1=read.table(all_recom_files[2],sep=",",header=T)
df2=read.table(all_recom_files[3],sep=",",header=T)
df3=read.table(all_recom_files[4],sep=",",header=T)

par(mfrow=c(2,2))
plot(df1$weighted_recomb,df2$weighted_recomb)
abline(a=0,b=1,col="red")
cor(df1$weighted_recomb,df2$weighted_recomb) # 0.89
plot(df1$weighted_recomb,df3$weighted_recomb)
abline(a=0,b=1,col="red")
cor(df1$weighted_recomb,df3$weighted_recomb) # 0.70
plot(df3$weighted_recomb,df2$weighted_recomb)
abline(a=0,b=1,col="red")
cor(df3$weighted_recomb,df2$weighted_recomb) # 0.85

means=rowMeans(cbind(df1$weighted_recomb,
                     df2$weighted_recomb,
                     df3$weighted_recomb))
dfmean = as.data.frame(cbind(df1$windowstarts,df1$windowstops,means,df1$chrom))
colnames(dfmean) = colnames(df1)

par(mfrow=c(2,2))
plot(df1$weighted_recomb,cex=0.2,type="l")
plot(df2$weighted_recomb,cex=0.2,type="l")
plot(df3$weighted_recomb,cex=0.2,type="l")
plot(dfmean$weighted_recomb,cex=0.2,type="l")

par(mfrow=c(1,1))
plot(dfmean$weighted_recomb,cex=0.2,type="l")
points(df1$weighted_recomb,cex=0.2,col=rgb(1,0,0,0.5),type="l")
points(df2$weighted_recomb,cex=0.2,col=rgb(0,1,0,0.5),type="l")
points(df3$weighted_recomb,cex=0.2,col=rgb(0,0,1,0.5),type="l")
points(dfmean$weighted_recomb,cex=0.2,type="l")
}

full_recom=NULL
for(recom_filename in all_recom_files){
  print(recom_filename)
#recom_filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION/Amphispiza-bilineata-called.geno.PREDICT.txt_w100000_o20000_genome.txt"
recom_genome = read.table(recom_filename,sep=",",header=T)
recom_genome$midPos = recom_genome$windowstarts+((recom_genome$windowstops - recom_genome$windowstarts)/2)
recom_genome$chr = gsub("PseudoNC_.+_Tgut_","",recom_genome$chrom)
shortspp=substr(strsplit(strsplit(basename(recom_filename),"\\.")[[1]][1],"-")[[1]][2],1,3)
recom_genome$species = shortspp
if(is.null(full_recom)){
  full_recom = recom_genome
} else {
  full_recom = rbind(full_recom,recom_genome)
}

}
write.table(full_recom,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigrecom_duplicates.may2020.txt",
            row.names = F)


write.table(full_recom,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigrecom.may2020.txt",
            row.names = F)

agg = merge(bigplot,full_recom,all=T)

write.table(agg,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.fst.dxy.taj.recom.may2020.txt",
                     row.names = F)

aggsmall = merge(bigplot,full_recom,all=F)


for(spp in sample(unique(aggsmall$species))){
  print(spp)
  png(paste(spp,".temp.may2020.png",sep=""))
  par(mfrow=c(3,1))
  aggtiny=aggsmall[aggsmall$species==spp,]
  ordertoplot = 1:nrow(aggtiny)
  plot(ordertoplot,aggtiny$Fst,col=as.numeric(as.factor(aggtiny$chr)))
  plot(ordertoplot,aggtiny$dxymeans,col=as.numeric(as.factor(aggtiny$chr)))
  #plot(aggtiny$midPos,aggtiny$Tajima,col=as.numeric(as.factor(aggtiny$chr)))
  plot(ordertoplot,aggtiny$weighted_recomb,col=as.numeric(as.factor(aggtiny$chr)))
  dev.off()
}

for(spp in sample(unique(agg$species))){
  print(spp)
  png(paste(spp,".temp2.may2020.png",sep=""))
  par(mfrow=c(3,1))
  aggtiny=agg[agg$species==spp,]
  ordertoplot = 1:nrow(aggtiny)
  plot(ordertoplot,aggtiny$Fst,col=as.numeric(as.factor(aggtiny$chr)))
  plot(ordertoplot,aggtiny$dxymeans,col=as.numeric(as.factor(aggtiny$chr)))
  #plot(aggtiny$midPos,aggtiny$Tajima,col=as.numeric(as.factor(aggtiny$chr)))
  plot(ordertoplot,aggtiny$weighted_recomb,col=as.numeric(as.factor(aggtiny$chr)))
  dev.off()
}

for(spp in unique(agg$species)){
  small=agg$Tajima[agg$species==spp]
  small=small[complete.cases(small)]
  if(length(small)!=0){
    print(spp)
  } 
  #plot(agg$Tajima[agg$species==spp])
}


## dxy fst etc
bigplot = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.fst.dxy.taj.recom.may2020.txt",
                   sep=" ",header = T,stringsAsFactors = F)
#boxplot(bigplot$dxymeans~bigplot$species)

bigplot = bigplot[!(grepl("Pseudo",bigplot$chr)),]

for(spp in unique(bigplot$species)) {
  #print(spp)
  small = bigplot[bigplot$species == spp, c("chr","Fst")]
  small = small[complete.cases(small),]
  cat(
    paste(
      mean(small[,2], na.rm = T),
  "+-",
  sd(small[,2], na.rm = T),
  " (",
  length(unique(small$chr)),
  ")", sep = ""),sep="\n")
}

meanagg = aggregate(
  cbind(Fst, dxymeans, weighted_recomb) ~ chr + species, data = bigplot,
  FUN=function(x){mean(x,na.rm=T)}
)
sdagg = aggregate(
  cbind(Fst, dxymeans, weighted_recomb) ~ chr + species, data = bigplot,
  FUN=function(x){sd(x,na.rm=T)}
)
  
par(mfrow=c(5,2),mar=c(2,2,1,0)) 
for(spp in unique(bigplot$species)){
  sppmean = meanagg[meanagg$species==spp,]
  sppsd = sdagg[sdagg$species==spp,]
  centers=barplot(sppmean$Fst,main=spp,las=2,cex.axis=0.5,
                  cex.names=0.75,ylim=c(0,0.125),
                  names=sppmean$chr)
  segments(centers, sppmean$Fst - sppsd$Fst * 2, centers,
           sppmean$Fst + sppsd$Fst * 2, lwd = 1.5)
}

par(mfrow=c(5,2),mar=c(2,2,1,0)) 
for(spp in unique(bigplot$species)){
  sppmean = meanagg[meanagg$species==spp,]
  sppsd = sdagg[sdagg$species==spp,]
  centers=barplot(sppmean$dxymeans,main=spp,las=2,cex.axis=0.5,
                  cex.names=0.75,ylim=c(0,0.04),
                  names=sppmean$chr)
  segments(centers, sppmean$dxymeans - sppsd$dxymeans * 2, centers,
           sppmean$dxymeans + sppsd$dxymeans * 2, lwd = 1.5)
}

par(mfrow=c(5,2),mar=c(2,3,1,0)) 
for(spp in unique(bigplot$species)){
  sppmean = meanagg[meanagg$species==spp,]
  sppsd = sdagg[sdagg$species==spp,]
  centers=barplot(sppmean$weighted_recomb,main=spp,las=2,cex.axis=0.5,
                  cex.names=0.75,ylim=c(0,1.5e-09),
                  names=sppmean$chr)
  segments(centers, sppmean$weighted_recomb - sppsd$weighted_recomb * 2, centers,
           sppmean$weighted_recomb + sppsd$weighted_recomb * 2, lwd = 1.5)
}

trimmed = bigplot[bigplot$species %in% c("bel","bru","cri","cur"),c("chr","species","Fst","dxymeans","weighted_recomb","sumquantile","plotorder")]
trimmed = trimmed[complete.cases(trimmed),]
trimmed = unique(trimmed)
cor(trimmed[,3:5])

cor(trimmed[trimmed$species=="bel",3:5])
cor(trimmed[trimmed$species=="bru",3:5])
cor(trimmed[trimmed$species=="cri",3:5])
cor(trimmed[trimmed$species=="cur",3:5])

png("fourspecies.corr.may2020.png")
par(mfrow=c(2,4))
plot(trimmed$weighted_recomb[trimmed$species=="bel"],trimmed$dxymeans[trimmed$species=="bel"])
plot(trimmed$weighted_recomb[trimmed$species=="bel"],trimmed$Fst[trimmed$species=="bel"])
plot(trimmed$weighted_recomb[trimmed$species=="bru"],trimmed$dxymeans[trimmed$species=="bru"])
plot(trimmed$weighted_recomb[trimmed$species=="bru"],trimmed$Fst[trimmed$species=="bru"])
plot(trimmed$weighted_recomb[trimmed$species=="cri"],trimmed$dxymeans[trimmed$species=="cri"])
plot(trimmed$weighted_recomb[trimmed$species=="cri"],trimmed$Fst[trimmed$species=="cri"])
plot(trimmed$weighted_recomb[trimmed$species=="cur"],trimmed$dxymeans[trimmed$species=="cur"])
plot(trimmed$weighted_recomb[trimmed$species=="cur"],trimmed$Fst[trimmed$species=="cur"])
dev.off()

png("temp.may2020.png",height=600,width=600)
par(mfrow=c(2,2))
boxplot(trimmed$weighted_recomb[trimmed$species=="bel"]~trimmed$sumquantile[trimmed$species=="bel"],
        col=c("white","white","darkred",
              "white","white","cyan",
              "white","white","magenta"))
boxplot(trimmed$weighted_recomb[trimmed$species=="bru"]~trimmed$sumquantile[trimmed$species=="bru"],
        col=c("white","white","darkred",
              "white","white","cyan",
              "white","white","magenta"))
boxplot(trimmed$weighted_recomb[trimmed$species=="cri"]~trimmed$sumquantile[trimmed$species=="cri"],
        col=c("white","white","darkred",
              "white","white","cyan",
              "white","white","magenta"))
boxplot(trimmed$weighted_recomb[trimmed$species=="cur"]~trimmed$sumquantile[trimmed$species=="cur"],
        col=c("white","white","darkred",
              "white","white","cyan",
              "white","white","magenta"))
dev.off()
## 9  = 1+8  = fst not outlier, dxy not outlier
## 10 = 2+8  = fst low,         dxy not outlier
## 12 = 4+8  = fst high,        dxy not outlier -- AMBIGUOUS
## 17 = 1+16 = fst not outlier, dxy low
## 18 = 2+16 = fst low,         dxy low
## 20 = 4+16 = fst high,        dxy low -- SWEEP
## 33 = 1+32 = fst not outlier, dxy high
## 34 = 2+32 = fst low,         dxy high
## 36 = 4+32 = fst high,        dxy high -- ISLAND

lowagg = aggregate(bigplot$plotorder~bigplot$chr,FUN=function(x){min(x,na.rm=T)})
highagg = aggregate(bigplot$plotorder~bigplot$chr,FUN=function(x){max(x,na.rm=T)})

together=cbind(lowagg,highagg)

macrochroms=c(1:10,"Z")

png("macro.may2020.png",width=900,height=900)
palette(c("lightgrey","grey"))
par(mfrow=c(4,1))
for(spp in unique(trimmed$species)){
  smalltrim=trimmed[trimmed$species==spp,]
  #smalltrim = smalltrim[smalltrim$chr %in% macrochroms,]
plot(smalltrim$plotorder,
     smalltrim$weighted_recomb,col=as.numeric(as.factor(smalltrim$chr)),cex=0.25,
     type="n")
rect(lowagg[,2], rep(-0.1,nrow(lowagg)), highagg[,2],
     rep(0.1,nrow(lowagg)), lwd = 1,col=c("lightgrey","grey"),xpd=F)
points(smalltrim$plotorder,
     smalltrim$weighted_recomb,col="black",cex=0.25)
points(smalltrim$plotorder[smalltrim$sumquantile==20],
       smalltrim$weighted_recomb[smalltrim$sumquantile==20],
       col="cyan",pch=16)
points(smalltrim$plotorder[smalltrim$sumquantile==36],
       smalltrim$weighted_recomb[smalltrim$sumquantile==36],
       col="magenta",pch=16)
}
dev.off()



for(spp in unique(trimmed$species)){
  spptrim=trimmed[trimmed$species==spp,]
  num=length(unique(spptrim$chr))
  png(paste(spp, ".may2020.png", sep = ""),
      width = 150 * ceiling(sqrt(num)),
      height = 150 * ceiling(sqrt(num)))
  par(mfrow=c(ceiling(sqrt(num)),ceiling(sqrt(num))))
  for(chrom in unique(spptrim$chr)){
    smalltrim=spptrim[spptrim$chr==chrom,]
    cat(spp,chrom,sep=" ")
    plot(smalltrim$plotorder,
         smalltrim$weighted_recomb,cex=0.25)
    points(smalltrim$plotorder[smalltrim$sumquantile==20],
           smalltrim$weighted_recomb[smalltrim$sumquantile==20],
           col="cyan",pch=16)
    points(smalltrim$plotorder[smalltrim$sumquantile==36],
           smalltrim$weighted_recomb[smalltrim$sumquantile==36],
           col="magenta",pch=16,cex=1.5)
  }
  dev.off()
}


yes_sweep = trimmed[trimmed$sumquantile==20,]
not_sweep = trimmed[trimmed$sumquantile!=20,]
yes_island = trimmed[trimmed$sumquantile==36,]
not_island = trimmed[trimmed$sumquantile!=36,]

t.test(yes_sweep$weighted_recomb,not_sweep$weighted_recomb) ## yes difference -- pval ~0, yes sweep 10.41-e9, no sweep 9.90 e-9
t.test(yes_island$weighted_recomb,not_island$weighted_recomb) ## no difference pval 0.15
t.test(yes_sweep$weighted_recomb,yes_island$weighted_recomb) ## yes difference pval 0.00039 yes sweep 10.41e-9 yes island 9.7 e-9

spp="bel"
t.test(yes_sweep$weighted_recomb[yes_sweep$species==spp],not_sweep$weighted_recomb[yes_sweep$species==spp]) ## yes difference -- pval ~0, yes sweep 10.41-e9, no sweep 9.90 e-9
t.test(yes_island$weighted_recomb[yes_sweep$species==spp],not_island$weighted_recomb[yes_sweep$species==spp]) ## no difference pval 0.15
t.test(yes_sweep$weighted_recomb[yes_sweep$species==spp],yes_island$weighted_recomb[yes_sweep$species==spp]) ## yes difference pval 0.00039 yes sweep 10.41e-9 yes island 9.7 e-9


#png("rectangles.png")

bigplot=unique(bigplot)

allsweeps = bigplot[bigplot$sumquantile==20,c("chr","species","plotorder","sumquantile","midPos")]
allsweeps=unique(allsweeps)
table(allsweeps$chr,allsweeps$species)

wheresweeps=table(allsweeps$plotorder)
plot(as.numeric(names(wheresweeps)),as.numeric(wheresweeps),
     type="p")

multiple_species=names(wheresweeps[wheresweeps>1])
multall=allsweeps[allsweeps$plotorder %in% multiple_species,]
unique(multall$chr)

allislands = bigplot[bigplot$sumquantile==36,c("chr","species","plotorder","sumquantile",
                                               "midPos")]
allislands=unique(allislands)
table(allislands$chr,allislands$species)

whereislands=table(allislands$plotorder)
plot(as.numeric(names(whereislands)),as.numeric(whereislands),
     type="p")

multiple_species=names(whereislands[whereislands>0])
multall=allislands[allislands$plotorder %in% multiple_species,]
unique(multall$chr)



trim_for_chrom = bigplot[,c("Fst","dxymeans","weighted_recomb","chr","species")]
chromagg = aggregate(cbind(Fst,dxymeans,weighted_recomb) ~ chr,data=trim_for_chrom,
                     FUN=function(x){mean(x,na.rm=T)})
chromsd = aggregate(cbind(Fst,dxymeans,weighted_recomb) ~ chr,data=trim_for_chrom,
                     FUN=function(x){sd(x,na.rm=T)})

cbind(chromagg,chromsd)
par(mfrow=c(3,1))
centers=barplot(chromagg$Fst,main="fst",las=2,cex.axis=0.5,
                  cex.names=0.75,
                  names=chromagg$chr)
segments(centers, chromagg$Fst- chromsd$Fst, centers,
         chromagg$Fst + chromsd$Fst, lwd = 1.5)
centers=barplot(chromagg$dxymeans,main="dxymeans",las=2,cex.axis=0.5,
                cex.names=0.75,
                names=chromagg$chr)
segments(centers, chromagg$dxymeans- chromsd$dxymeans, centers,
         chromagg$dxymeans + chromsd$dxymeans, lwd = 1.5)
centers=barplot(chromagg$weighted_recomb,main="weighted_recomb",las=2,cex.axis=0.5,
                cex.names=0.75,
                names=chromagg$chr)
segments(centers, chromagg$weighted_recomb- chromsd$weighted_recomb, centers,
         chromagg$weighted_recomb + chromsd$weighted_recomb, lwd = 1.5)
