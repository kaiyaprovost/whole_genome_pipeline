## YOU ARE GOING TO NEED TO REORGANIZE THIS 
skip=T
setwd("~")
if(skip == F){
merged=read.table("~/rec_taj_dxy_fst_islswp_missing.temp",fill=T,header=T,stringsAsFactors = F)
small=merged[,c("chr","midPos","species")]
small=small[duplicated(small),]
dim(small)
dim(merged)
merged=merged[order(merged$chr,merged$midPos,merged$species),]
}

## sweeps and islands, compiled, without missing chromosomes
library(plyr)
library(readr)

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

if(skip == F){

fstfiles=list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/FST/3.slidingWindowJobs/",
                    pattern="fst$",full.names = T)

fstdf=NULL
for(fstfile in fstfiles) {
  spp=substr(strsplit(basename(fstfile),"_")[[1]][4],1,3)
  if(is.null(fstdf)){
    fstdf = read.table(fstfile,header=T,sep="\t",stringsAsFactors = F,fill = T)
    fstdf$species = spp
  } else {
    toadd = read.table(fstfile,header=T,sep="\t",stringsAsFactors = F,fill = T)
    toadd$species = spp
    fstdf = rbind(fstdf,toadd)
  }
}

fst_agg = aggregate(fstdf$Fst~fstdf$chr+fstdf$species,FUN=function(x){mean(x,na.rm=T)})
summary(fst_agg)
#boxplot(fst_agg$`fstdf$Fst`)
fst_agg = aggregate(fstdf$Fst~fstdf$chr,FUN=function(x){mean(x,na.rm=T)})
summary(fst_agg)
#boxplot(fstdf$Fst~fstdf$chr)
write.table(fstdf,"fst.temp")



dxyfiles = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY",
                      pattern="_4_SON_Dxy_WINDOWS_chrfix_1-ALL.txt$",full.names = T,recursive = T)

dxydf=NULL
for(dxyfile in dxyfiles) {
  spp=strsplit(basename(dxyfile),"_")[[1]][1]
  if(is.null(dxydf)){
    dxydf = read.table(dxyfile,header=T,sep=" ",stringsAsFactors = F,fill=T)
    dxydf$species = spp
  } else {
    toadd = read.table(dxyfile,header=T,sep=" ",stringsAsFactors = F,fill=T)
    toadd$species = spp
    dxydf = rbind(dxydf,toadd)
  }
}
dxydf$midPos = as.numeric(dxydf$starts) + 59999
dxydf$chr = dxydf$scafs

write.table(dxydf,"dxy.temp")

#dxydf=dxydf[dxydf$species=="bel",]

#merged$dxymeans=NA

dxy_agg = aggregate(dxydf$means~dxydf$scafs+dxydf$species,FUN=function(x){mean(x,na.rm=T)})
summary(dxy_agg)

dxyfst = merge(fstdf[,c("chr","midPos","species","Fst")],
               dxydf[,c("chr","midPos","species","means")],
               all=T,sort=T)
write.table(dxyfst,"fst_dxy.temp",row.names = F)
dxyfst = read.table("fst_dxy.temp",header=T,stringsAsFactors = F)

chroms_with_both = dxyfst[complete.cases(dxyfst[,c("Fst","means")]),]

# z = (i - mean) / sd

dxyfst$rankFST = dplyr::percent_rank(dxyfst$Fst)
dxyfst$rankDXY = dplyr::percent_rank(dxyfst$means)
#dxyfst$zFST = (dxyfst$Fst - mean(dxyfst$Fst,na.rm=T)) / sd(dxyfst$Fst,na.rm=T)
#dxyfst$zDXY = (dxyfst$means - mean(dxyfst$means,na.rm=T)) / sd(dxyfst$means,na.rm=T)

#png("rank_z.png")
#par(mfrow=c(1,2))
#plot(dxyfst$rankFST,dxyfst$zFST)
#plot(dxyfst$rankDXY,dxyfst$zDXY)
#dev.off()

dxyfst = dxyfst[complete.cases(dxyfst),]

dxyfst$rankISLAND = as.numeric(dxyfst$rankFST >= 0.95 & dxyfst$rankDXY >= 0.95)
dxyfst$rankSWEEP = as.numeric(dxyfst$rankFST >= 0.95 & dxyfst$rankDXY <= 0.05)
#dxyfst$zISLAND = as.numeric(dxyfst$zFST >= 5 & dxyfst$zDXY >= 5)
#dxyfst$zSWEEP = as.numeric(dxyfst$zFST >= 5 & dxyfst$zDXY <= -5)
table(dxyfst[,c("rankISLAND","rankSWEEP")])
#table(dxyfst[,c("zISLAND","zSWEEP")])

dxyfst$ranksppFST = dplyr::percent_rank(dxyfst$Fst)
dxyfst$ranksppDXY = dplyr::percent_rank(dxyfst$means)
#dxyfst$zsppFST = (dxyfst$Fst - mean(dxyfst$Fst,na.rm=T)) / sd(dxyfst$Fst,na.rm=T)
#dxyfst$zsppDXY = (dxyfst$means - mean(dxyfst$means,na.rm=T)) / sd(dxyfst$means,na.rm=T)

for(spp in unique(dxyfst$species)){
  dxyfst$ranksppFST[dxyfst$species==spp] = dplyr::percent_rank(dxyfst$Fst[dxyfst$species==spp])
  dxyfst$ranksppDXY[dxyfst$species==spp] = dplyr::percent_rank(dxyfst$means[dxyfst$species==spp])
  #dxyfst$zsppFST[dxyfst$species==spp] = (dxyfst$Fst[dxyfst$species==spp] - mean(dxyfst$Fst[dxyfst$species==spp],na.rm=T)) / sd(dxyfst$Fst[dxyfst$species==spp],na.rm=T)
  #dxyfst$zsppDXY[dxyfst$species==spp] = (dxyfst$means[dxyfst$species==spp] - mean(dxyfst$means[dxyfst$species==spp],na.rm=T)) / sd(dxyfst$means[dxyfst$species==spp],na.rm=T)
}

#png("rank_z_spp.png")
#par(mfrow=c(2,2))
#plot(dxyfst$ranksppFST,dxyfst$zsppFST)
#plot(dxyfst$ranksppDXY,dxyfst$zsppDXY)
#plot(dxyfst$rankFST,dxyfst$ranksppFST)
#plot(dxyfst$zDXY,dxyfst$zsppDXY)
#dev.off()

dxyfst$ranksppISLAND = as.numeric(dxyfst$ranksppFST >= 0.95 & dxyfst$ranksppDXY >= 0.95)
dxyfst$ranksppSWEEP = as.numeric(dxyfst$ranksppFST >= 0.95 & dxyfst$ranksppDXY <= 0.05)
#dxyfst$zsppISLAND = as.numeric(dxyfst$zsppFST >= 5 & dxyfst$zsppDXY >= 5)
#dxyfst$zsppSWEEP = as.numeric(dxyfst$zsppFST >= 5 & dxyfst$zsppDXY <= -5)
table(dxyfst[,c("ranksppISLAND","ranksppSWEEP")])
#table(dxyfst[,c("zsppISLAND","zSWEEP")])
table(dxyfst[,c("ranksppSWEEP","rankSWEEP")])
table(dxyfst[,c("ranksppISLAND","rankISLAND")])

dxyfst$sppchr = paste(dxyfst$chr,dxyfst$species,sep="-")
table(dxyfst[dxyfst$ranksppSWEEP==1,c("species","ranksppSWEEP")])
table(dxyfst[dxyfst$ranksppSWEEP==1,c("chr","ranksppSWEEP")])
table(dxyfst[dxyfst$ranksppSWEEP==1,c("sppchr","ranksppSWEEP")])
table(dxyfst[dxyfst$ranksppSWEEP==1,c("chr","species")])

table(dxyfst[dxyfst$ranksppISLAND==1,c("chr","species")])

table(rbind(dxyfst$ranksppISLAND | dxyfst$ranksppSWEEP,dxyfst$species))

table(dxyfst$species)

recdf = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigrecom.may2020.txt",
                   header=T,sep=" ",stringsAsFactors = T,fill=T)

## per chrom per spp:
min(recdf$weighted_recomb,na.rm=T)

rec_agg=aggregate(recdf$weighted_recomb~recdf$chrom+recdf$species,FUN=function(x){mean(x,na.rm=T)})
summary(rec_agg)
rec_agg=aggregate(recdf$weighted_recomb~recdf$species,FUN=function(x){mean(x,na.rm=T)})
summary(rec_agg)

boxplot(recdf$weighted_recomb ~ recdf$chrom,outline=F)
mod=aov(recdf$weighted_recomb ~ recdf$chrom)
summary(mod)
tukey=TukeyHSD(mod)$`recdf$chrom`
sigs = rownames(tukey)[tukey[,4]< 0.00000001]
tukey[order(-tukey[,4]),]


tajdf = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigtaj.new.june2020.txt",
                   header=T,sep=",",stringsAsFactors = T,fill=T)
tajdf$chr = tajdf$Chr
tajdf$midPos = tajdf$WinCenter
unique(tajdf$species)
summary(tajdf)

## get new tajd data into old merged file
#intersect(names(tajdf),names(merged))

## first get the unique chr-window values
## fst is midpos
## dxy is starts, 1 indexed
## rec is starts 0 indexed, stops, midpos
## taj is wincenter

## extract one species
#spp = "bel"
#tajdf_spp = tajdf[tajdf$species==spp,]
#recdf_spp = recdf[recdf$species==spp,]
#recdf_spp <- recdf_spp[order(recdf_spp$chr, recdf_spp$midPos),]





## need to fix the overlaps, weighted averages
# windowsRecalculator(data,recomb_file,window_size,overlap,reset_overlap=F,overwrite=F)

windows = unique(rbind(fstdf[,c("chr","midPos")],
                       dxydf[,c("chr","midPos")],
                       tajdf[,c("chr","midPos")],
                       recdf[,c("chr","midPos")]))
windows <- windows[order(windows$chr, windows$midPos),]

dim(recdf)
dim(tajdf)
#tajdf = tajdf[,c("chr","midPos","species","Tajima")]
#tajdf = aggregate(tajdf$Tajima ~ tajdf$chr + tajdf$midPos + tajdf$species,
#                      FUN=function(x){mean(x,na.rm=T)})
#colnames(tajdf) = c("chr","midPos","species","Tajima")
tajdf <- tajdf[order(tajdf$chr, tajdf$midPos),]
recdf <- recdf[order(recdf$chr, recdf$midPos),]


merged = merge(recdf[,c("chr","midPos","species","weighted_recomb")],
               tajdf[,c("chr","midPos","species","Tajima")],
               all=T,sort=T)
merged <- merged[order(merged$chr, merged$midPos),]
merged = merged[merged$midPos %% 10000 == 0,]

merged$windowstarts = merged$midPos-50000
merged$windowstops = merged$midPos+50000
write.table(merged,"rec_taj.temp",row.names = F)


recdf = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigrecom.may2020.txt",
                   header=T,sep=" ",stringsAsFactors = T,fill=T)



for(row_index in 1:nrow(merged)){
  
  if(row_index %% 1000 == 0) {cat(paste(row_index,"/",nrow(merged)),sep="\n")}
  
  if(is.na(merged$weighted_recomb[row_index])){
    this_start = merged$windowstarts[row_index]
    this_stop = merged$windowstops[row_index]
    this_mid = merged$midPos[row_index]
    this_chr = merged$chr[row_index]
    this_species = merged$species[row_index]
    
    to_average=recdf[(recdf$chr == this_chr 
                      & recdf$windowstarts <= this_mid 
                      & recdf$windowstops >= this_mid
                      & recdf$species == this_species),]
    
    to_average$overlap[to_average$windowstarts < this_start] = 
      to_average$windowstops[to_average$windowstarts < this_start] - this_start
    to_average$overlap[to_average$windowstarts >= this_start] = 
      this_stop - to_average$windowstarts[to_average$windowstarts >= this_start]
    
    this_rec = sum(to_average$weighted_recomb * to_average$overlap) / sum(to_average$overlap)
    merged$weighted_recomb[row_index] = this_rec
    
  }
}

write.table(merged,"rec_taj.temp",row.names = F)
summary(merged$weighted_recomb)
summary(merged$Tajima)
merged <- merged[order(merged$species, merged$chr, merged$midPos),]


merged = read.table("rec_taj.temp",fill=T,header=T,stringsAsFactors = F)
dxydf = read.table("dxy.temp",fill=T,header=T,stringsAsFactors = F)
fstdf = read.table("fst.temp",fill=T,header=T,stringsAsFactors = F)


colnames(dxydf)[3] = "dxymeans"

merged = merge(merged,
               dxydf[,c("chr","midPos","species","dxymeans")],
               all=T,sort=T)
summary(merged$dxymeans)

write.table(merged,"rec_taj_dxy.temp",row.names = F)
}

merged = read.table("rec_taj_dxy.temp",fill=T,header=T,stringsAsFactors = F)
fstdf = read.table("fst.temp",fill=T,header=T,stringsAsFactors = F)
recdf = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigrecom.may2020.txt",
                   header=T,sep=" ",stringsAsFactors = T,fill=T)

merged = merge(merged,
               fstdf[,c("chr","midPos","species","Fst")],
               all=T,sort=T)
summary(merged$Fst)

merged$windowstarts = merged$midPos-50000
merged$windowstops = merged$midPos+50000

for(row_index in 1235000:nrow(merged)){
  
  if(row_index %% 100 == 0) {cat(paste(row_index,"/",nrow(merged)),sep="\n")}
  
  if(is.na(merged$weighted_recomb[row_index])){
    this_start = merged$windowstarts[row_index]
    this_stop = merged$windowstops[row_index]
    this_mid = merged$midPos[row_index]
    this_chr = merged$chr[row_index]
    this_species = merged$species[row_index]
    
    to_average=recdf[(recdf$chr == this_chr 
                      & recdf$windowstarts <= this_mid 
                      & recdf$windowstops >= this_mid
                      & recdf$species == this_species),]
    
    to_average$overlap[to_average$windowstarts < this_start] = 
      to_average$windowstops[to_average$windowstarts < this_start] - this_start
    to_average$overlap[to_average$windowstarts >= this_start] = 
      this_stop - to_average$windowstarts[to_average$windowstarts >= this_start]
    
    this_rec = sum(to_average$weighted_recomb * to_average$overlap) / sum(to_average$overlap)
    merged$weighted_recomb[row_index] = this_rec
    
  }
}

merged$plotorder = 1
merged = merged[order(merged$chr,merged$midPos,merged$species),]
totrows=nrow(merged)
for(rowindex in 2:totrows){
  if(rowindex %% 1000 == 0){print(paste(rowindex,"/",totrows))}
  this_chr = merged$chr[rowindex]
  prev_chr = merged$chr[rowindex-1]
  this_pos = merged$midPos[rowindex]
  prev_pos = merged$midPos[rowindex-1]
  prev_ord = merged$plotorder[rowindex-1]
  
  if(this_chr==prev_chr && this_pos==prev_pos){
    merged$plotorder[rowindex] = prev_ord
  } else {
    merged$plotorder[rowindex] = prev_ord+1
  }
  
}


write.table(merged,"rec_taj_dxy_fst.temp",row.names = F)
write.table(merged[complete.cases(merged),],"rec_taj_dxy_fst.complete.temp",row.names = F)


missingfiles = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/MISSING",
                          pattern="geno.vcf.lmiss.lmiss.windows.txt",full.names = T,
                          recursive = T)

missingdf=NULL
for(missingfile in missingfiles) {
  spp=substr(strsplit(basename(missingfile),"-")[[1]][2],1,3)
  print(spp)
  if(is.null(missingdf)){
    missingdf = read.table(missingfile,header=T,sep="\t",stringsAsFactors = F,fill = T)
    missingdf$species = spp
  } else {
    toadd = read.table(missingfile,header=T,sep="\t",stringsAsFactors = F,fill = T)
    toadd$species = spp
    missingdf = merge(missingdf,toadd,all=T)
  }
}
head(missingdf)
colnames(missingdf)[1:2] = c("chr","midPos")
missingdf = unique(missingdf)
dim(missingdf)
missingdf=missingdf[,c(1:11)]
missingdf=missingdf[complete.cases(missingdf),]
missingdf$midPos=as.numeric(missingdf$midPos)
missingdf=missingdf[missingdf$midPos %% 10000 == 0,]
missingdf$chr=missingdf
missingdf=missingdf[order(missingdf$chr,missingdf$midPos,missingdf$species),]
write.table(missingdf,"~/missing.temp",row.names=F,quote=F,sep="\t")


mergx=merge(merged,missingdf,by=c("chr","midPos","species"),all.x=T,
            all.y=F)
mergy=merge(merged,missingdf,by=c("chr","midPos","species"),all.x=F,
            all.y=T)
merg = unique(rbind(mergx,mergy))
merg=unique(merg)
write.table(merg,"~/rec_taj_dxy_fst_islswp_missing.temp",row.names=F,quote=F,sep="\t")
merged=read.table("~/rec_taj_dxy_fst_islswp_missing.temp",fill=T,header=T,stringsAsFactors = F)


## adding DXY and TAJ to existing data
smalldxydf = unique(dxydf[,c("species","chr","midPos","means")])
colnames(smalldxydf) = c("species","chr","midPos","dxymeans")
smalldxydf=unique(smalldxydf)
merged=unique(merged)
for(this_spp in sort(unique(smalldxydf$species))){
  print(sum(is.na(merged$dxymeans))) # 1246927
  print(this_spp)
  if(length(is.na(merged$dxymeans[(merged$species==this_spp)])) > 0) {
    subset1 = smalldxydf[smalldxydf$species==this_spp,]
    subset1 = subset1[complete.cases(subset1),]
    #for(this_chr in rev(sort(unique(subset1$chr)))){
    for(this_chr in c("14","Z","4")){
      print(sum(is.na(merged$dxymeans))) # 1246927
      print(this_chr)
      if(length(is.na(merged$dxymeans[(merged$species==this_spp & merged$chr==this_chr)])) > 0) {
        subset2 = subset1[subset1$chr==this_chr,]
        subset2 = subset2[complete.cases(subset2),]
        print(nrow(subset2))
        if(nrow(subset2)>0){
          for(pos_index in 1:length(unique(subset2$midPos))){
            if(pos_index %% 100==0){print(paste(pos_index,"of",length(unique(subset2$midPos))))}
            this_pos = unique(subset2$midPos)[pos_index]
            #tocheck=merged$dxymeans[(merged$species==this_spp & merged$chr==this_chr & merged$midPos==this_pos)]
            
            #if(is.na(merged$dxymeans[(merged$species==this_spp & merged$chr==this_chr & merged$midPos==this_pos)])){
            #  print("CHANGE")
            #}
            
            #if(is.na(tocheck)) {
            this_dxy = subset2$dxymeans[subset2$midPos==this_pos]
            #tocheck=merged$dxymeans[(merged$species==this_spp & merged$chr==this_chr & merged$midPos==this_pos)]
            merged$dxymeans[(merged$species==this_spp & merged$chr==this_chr & merged$midPos==this_pos)] <- this_dxy
          }
        }
      }
    }
    write.table(merged,"~/rec_taj_dxy_fst_islswp_missing.temp")
  }
  write.table(merged,"~/rec_taj_dxy_fst_islswp_missing.temp")
}
sum(is.na(merged$dxymeans)) # 1246927
merged=unique(merged)
write.table(merged,"~/rec_taj_dxy_fst_islswp_missing.temp")
smalltajdf = unique(tajdf[,c("species","chr","midPos","Tajima")])
for(this_spp in sort(unique(smalltajdf$species))){
  print(sum(is.na(merged$Tajima))) # 1246927
  print(this_spp)
  if((this_spp %in% c("sin")) && length(is.na(merged$Tajima[(merged$species==this_spp)])) > 0) {
    subset1 = smalltajdf[smalltajdf$species==this_spp,]
    for(this_chr in sort(unique(subset1$chr))){
      print(sum(is.na(merged$Tajima))) # 1246927
      print(this_chr)
      if(length(is.na(merged$Tajima[(merged$species==this_spp & merged$chr==this_chr)])) > 0) {
        subset2 = subset1[subset1$chr==this_chr,]
        print(nrow(subset2))
        for(pos_index in 1:length(unique(subset2$midPos))){
          if(pos_index %% 100==0){print(paste(pos_index,"of",length(unique(subset2$midPos))))}
          this_pos = unique(subset2$midPos)[pos_index]
          #tocheck=merged$Tajima[(merged$species==this_spp & merged$chr==this_chr & merged$midPos==this_pos)]
          
          #if(is.na(merged$Tajima[(merged$species==this_spp & merged$chr==this_chr & merged$midPos==this_pos)])){
          #  print("CHANGE")
          #}
          
          #if(is.na(tocheck)) {
          this_taj = subset2$Tajima[subset2$midPos==this_pos]
          #tocheck=merged$Tajima[(merged$species==this_spp & merged$chr==this_chr & merged$midPos==this_pos)]
          merged$Tajima[(merged$species==this_spp & merged$chr==this_chr & merged$midPos==this_pos)] <- this_taj
          #}
        }
      }
    }
    write.table(merged,"~/rec_taj_dxy_fst_islandssweeps.temp")
  }
  write.table(merged,"~/rec_taj_dxy_fst_islandssweeps.temp")
}
sum(is.na(merged$Tajima)) # 1246927
write.table(merged,"~/rec_taj_dxy_fst_islandssweeps.temp")

merged = read.table("rec_taj_dxy_fst.temp",header=T,stringsAsFactors = F,fill=T)
dim(merged)
merged=unique(merged)

dxyfst = merged[complete.cases(merged[,c("Fst","dxymeans")]),]

merged$rankFST = dplyr::percent_rank(merged$Fst)
merged$rankDXY = dplyr::percent_rank(merged$dxymeans)
merged$rankISLAND = as.numeric(merged$rankFST >= 0.95 & merged$rankDXY >= 0.95)
merged$rankSWEEP = as.numeric(merged$rankFST >= 0.95 & merged$rankDXY <= 0.05)
table(merged[,c("rankISLAND","rankSWEEP")])
merged$ranksppFST = dplyr::percent_rank(merged$Fst)
merged$ranksppDXY = dplyr::percent_rank(merged$dxymeans)

for(spp in unique(merged$species)){
  merged$ranksppFST[merged$species==spp] = dplyr::percent_rank(merged$Fst[merged$species==spp])
  merged$ranksppDXY[merged$species==spp] = dplyr::percent_rank(merged$dxymeans[merged$species==spp])
}


merged$ranksppISLAND = as.numeric(merged$ranksppFST >= 0.95 & merged$ranksppDXY >= 0.95)
merged$ranksppSWEEP = as.numeric(merged$ranksppFST >= 0.95 & merged$ranksppDXY <= 0.05)
table(merged[,c("ranksppISLAND","ranksppSWEEP")])
table(merged[,c("ranksppSWEEP","rankSWEEP")])
table(merged[,c("ranksppISLAND","rankISLAND")])

merged$sppchr = paste(merged$chr,merged$species,sep="-")
table(merged[merged$ranksppSWEEP==1,c("species","ranksppSWEEP")])
table(merged[merged$ranksppSWEEP==1,c("chr","ranksppSWEEP")])
table(merged[merged$ranksppSWEEP==1,c("sppchr","ranksppSWEEP")])
table(merged[merged$ranksppSWEEP==1,c("chr","species")])

table(merged[merged$ranksppISLAND==1,c("chr","species")])

table(rbind(merged$ranksppISLAND | merged$ranksppSWEEP,merged$species))

table(merged$species)

# merged$plotorder = paste(merged$chr,merged$midPos)
# uniques=unique(merged$plotorder)
# numuniques = uniques
# 
# lapply(unique(merged$plotorder),FUN=function(x){
#   print(x)
#   merged$numorder[merged$plotorder==x] = which(uniques %in% x)
# })


write.table(unique(merged),"rec_taj_dxy_fst_islandssweeps.temp")
merged=read.table("~/rec_taj_dxy_fst_islandssweeps.temp",header=T,fill=T,
                  stringsAsFactors = F)
dim(merged)








## plot dxy for each species
micro=T
Z1=F
smallspp=F
{
  dxymerged = merged[merged$chr %in% unique(merged$chr[!(is.na(merged$dxymeans))]),]
  if(micro==T){
    dxymerged = dxymerged[dxymerged$chr %in% c("1B",as.character(10:28),"mtDNA","LG2","LG5","LGE22"),] ## micro
    suffix="micro"
  } else {
    if(Z1==F){
      dxymerged = dxymerged[dxymerged$chr %in% c("1","1A",as.character(2:4),"4A",as.character(5:9),"Z"),] ## macro
      suffix="macro"
    } else {
      dxymerged = dxymerged[dxymerged$chr %in% c("1","Z","LGE22"),] ## macro
      suffix="Z1LGE22"
      
    }
  }
  dxymerged$color=as.numeric(as.factor(dxymerged$chr))
  dxymerged = unique(dxymerged[!(is.na(dxymerged$dxymeans)),])
  min_plot_per_chrom = aggregate(dxymerged$plotorder~dxymerged$chr,FUN=function(x){min(x,na.rm=T)}); colnames(min_plot_per_chrom)=c("chr","min")
  max_plot_per_chrom = aggregate(dxymerged$plotorder~dxymerged$chr,FUN=function(x){max(x,na.rm=T)}); colnames(max_plot_per_chrom)=c("chr","max")
  tot_plot_per_chrom = merge(min_plot_per_chrom,max_plot_per_chrom,all=F)
  dxymerged$plotindex=dxymerged$plotorder
  #print(tot_plot_per_chrom)
  for(rowindex in 2:nrow(tot_plot_per_chrom)){
    this_chr = tot_plot_per_chrom$chr[rowindex]
    #print(this_chr)
    num_to_subtract = (tot_plot_per_chrom$min[rowindex] - tot_plot_per_chrom$max[rowindex-1])-100
    #print(num_to_subtract)
    dxymerged$plotindex[dxymerged$chr==this_chr] = (dxymerged$plotorder[dxymerged$chr==this_chr])-num_to_subtract
    tot_plot_per_chrom[rowindex,2:3] = tot_plot_per_chrom[rowindex,2:3]-num_to_subtract
    
  }
  dxymerged$plotorder=dxymerged$plotindex
  min_plot_per_chrom = tot_plot_per_chrom[,c(1,2)]
  #print(tot_plot_per_chrom)
  specieslist = sort(unique(dxymerged$species[!(is.na(dxymerged$dxymeans))]))
  if(smallspp==T){
    specieslist = specieslist[specieslist %in% c("cur","bru","bil")]
    suffix=paste(suffix,"-small",sep="")
  }
  
  print(specieslist)
  png(paste("dxy_all_species_",suffix,".png",sep=""),width=800,height=(150*length(specieslist)))
  palette(c("black","grey"))
  par(mfrow=c(length(specieslist),1),
      mar=c(2,4,0,0))
  for(i in 1:length(specieslist)){
    spp=specieslist[i]
    print(spp)
    small = dxymerged[dxymerged$species==spp,c("chr","species","dxymeans","midPos","plotorder","plotindex","color")]
    small = unique(small)
    print(unique(small$chr))
    plot(dxymerged$plotorder,dxymerged$dxymeans,col=dxymerged$color,xaxt="n",ylab=spp,type="n")
    abline(v=min_plot_per_chrom[,2],col=c("darkred","red"))
    points(small$plotorder,small$dxymeans,col=small$color,xaxt="n",xaxt="n")
    axis(1,at=min_plot_per_chrom[,2],labels=min_plot_per_chrom[,1],las=2,cex=0.7)
  }
  dev.off()
}
## same as abuve for recomb
{
  recmerged = merged[merged$chr %in% unique(merged$chr[!(is.na(merged$weighted_recomb))]),]
  if(micro==T){
    recmerged = recmerged[recmerged$chr %in% c("1B",as.character(10:28),"mtDNA","LG2","LG5","LGE22"),] ## micro
    suffix="micro"
  } else {
    if(Z1==F){
      recmerged = recmerged[recmerged$chr %in% c("1","1A",as.character(2:4),"4A",as.character(5:9),"Z"),] ## macro
      suffix="macro"
    } else {
      recmerged = recmerged[recmerged$chr %in% c("1","Z","LGE22"),] ## macro
      suffix="Z1LGE22"
      
    }
  }
  recmerged$color=as.numeric(as.factor(recmerged$chr))
  recmerged = unique(recmerged[!(is.na(recmerged$weighted_recomb)),])
  min_plot_per_chrom = aggregate(recmerged$plotorder~recmerged$chr,FUN=function(x){min(x,na.rm=T)}); colnames(min_plot_per_chrom)=c("chr","min")
  max_plot_per_chrom = aggregate(recmerged$plotorder~recmerged$chr,FUN=function(x){max(x,na.rm=T)}); colnames(max_plot_per_chrom)=c("chr","max")
  tot_plot_per_chrom = merge(min_plot_per_chrom,max_plot_per_chrom,all=F)
  recmerged$plotindex=recmerged$plotorder
  #print(tot_plot_per_chrom)
  for(rowindex in 2:nrow(tot_plot_per_chrom)){
    this_chr = tot_plot_per_chrom$chr[rowindex]
    #print(this_chr)
    num_to_subtract = (tot_plot_per_chrom$min[rowindex] - tot_plot_per_chrom$max[rowindex-1])-100
    #print(num_to_subtract)
    recmerged$plotindex[recmerged$chr==this_chr] = (recmerged$plotorder[recmerged$chr==this_chr])-num_to_subtract
    tot_plot_per_chrom[rowindex,2:3] = tot_plot_per_chrom[rowindex,2:3]-num_to_subtract
    
  }
  recmerged$plotorder=recmerged$plotindex
  min_plot_per_chrom = tot_plot_per_chrom[,c(1,2)]
  #print(tot_plot_per_chrom)
  specieslist = sort(unique(recmerged$species[!(is.na(recmerged$weighted_recomb))]))
  if(smallspp==T){
    specieslist = specieslist[specieslist %in% c("cur","bru","nit")]
    suffix=paste(suffix,"-small",sep="")
  }
  
  print(specieslist)
  png(paste("rec_all_species_",suffix,".png",sep=""),width=800,height=(150*length(specieslist)))
  palette(c("black","grey"))
  par(mfrow=c(length(specieslist),1),
      mar=c(2,4,0,0))
  for(i in 1:length(specieslist)){
    spp=specieslist[i]
    print(spp)
    small = recmerged[recmerged$species==spp,c("chr","species","weighted_recomb","midPos","plotorder","plotindex","color")]
    small = unique(small)
    print(unique(small$chr))
    plot(recmerged$plotorder,recmerged$weighted_recomb,col=recmerged$color,xaxt="n",ylab=spp,type="n")
    abline(v=min_plot_per_chrom[,2],col=c("darkred","red"))
    points(small$plotorder,small$weighted_recomb,col=small$color,xaxt="n",xaxt="n")
    axis(1,at=min_plot_per_chrom[,2],labels=min_plot_per_chrom[,1],las=2,cex=0.7)
  }
  dev.off()
}
## same as above for fst
{
  fstmerged = merged[merged$chr %in% unique(merged$chr[!(is.na(merged$Fst))]),]
  if(micro==T){
    fstmerged = fstmerged[fstmerged$chr %in% c("1B",as.character(10:28),"mtDNA","LG2","LG5","LGE22"),] ## micro
    suffix="micro"
  } else {
    if(Z1==F){
      fstmerged = fstmerged[fstmerged$chr %in% c("1","1A",as.character(2:4),"4A",as.character(5:9),"Z"),] ## macro
      suffix="macro"
    } else {
      fstmerged = fstmerged[fstmerged$chr %in% c("1","Z","LGE22"),] ## macro
      suffix="Z1LGE22"
      
    }
  }
  fstmerged$color=as.numeric(as.factor(fstmerged$chr))
  fstmerged = unique(fstmerged[!(is.na(fstmerged$Fst)),])
  min_plot_per_chrom = aggregate(fstmerged$plotorder~fstmerged$chr,FUN=function(x){min(x,na.rm=T)}); colnames(min_plot_per_chrom)=c("chr","min")
  max_plot_per_chrom = aggregate(fstmerged$plotorder~fstmerged$chr,FUN=function(x){max(x,na.rm=T)}); colnames(max_plot_per_chrom)=c("chr","max")
  tot_plot_per_chrom = merge(min_plot_per_chrom,max_plot_per_chrom,all=F)
  fstmerged$plotindex=fstmerged$plotorder
  #print(tot_plot_per_chrom)
  for(rowindex in 2:nrow(tot_plot_per_chrom)){
    this_chr = tot_plot_per_chrom$chr[rowindex]
    #print(this_chr)
    num_to_subtract = (tot_plot_per_chrom$min[rowindex] - tot_plot_per_chrom$max[rowindex-1])-100
    #print(num_to_subtract)
    fstmerged$plotindex[fstmerged$chr==this_chr] = (fstmerged$plotorder[fstmerged$chr==this_chr])-num_to_subtract
    tot_plot_per_chrom[rowindex,2:3] = tot_plot_per_chrom[rowindex,2:3]-num_to_subtract
    
  }
  fstmerged$plotorder=fstmerged$plotindex
  min_plot_per_chrom = tot_plot_per_chrom[,c(1,2)]
  #print(tot_plot_per_chrom)
  specieslist = sort(unique(fstmerged$species[!(is.na(fstmerged$Fst))]))
  if(smallspp==T){
    specieslist = specieslist[specieslist %in% c("cur","bru","nit")]
    suffix=paste(suffix,"-small",sep="")
  }
  
  print(specieslist)
  png(paste("fst_all_species_",suffix,".png",sep=""),width=800,height=(150*length(specieslist)))
  palette(c("black","grey"))
  par(mfrow=c(length(specieslist),1),
      mar=c(2,4,0,0))
  for(i in 1:length(specieslist)){
    spp=specieslist[i]
    print(spp)
    small = fstmerged[fstmerged$species==spp,c("chr","species","Fst","midPos","plotorder","plotindex","color")]
    small = unique(small)
    print(unique(small$chr))
    plot(fstmerged$plotorder,fstmerged$Fst,col=fstmerged$color,xaxt="n",ylab=spp,type="n")
    abline(v=min_plot_per_chrom[,2],col=c("darkred","red"))
    points(small$plotorder,small$Fst,col=small$color,xaxt="n",xaxt="n")
    axis(1,at=min_plot_per_chrom[,2],labels=min_plot_per_chrom[,1],las=2,cex=0.7)
  }
  dev.off()
}
## same as above for sweeps and islands
{
  fstmerged = merged[merged$chr %in% unique(merged$chr[!(is.na(merged$Fst))]),]
  sweepmerged = fstmerged[fstmerged$chr %in% unique(fstmerged$chr[!(is.na(fstmerged$dxymeans))]),]
  
  if(micro==T){
    sweepmerged = sweepmerged[sweepmerged$chr %in% c("1B",as.character(10:28),"mtDNA","LG2","LG5","LGE22"),] ## micro
    suffix="micro"
  } else {
    if(Z1==F){
      sweepmerged = sweepmerged[sweepmerged$chr %in% c("1","1A",as.character(2:4),"4A",as.character(5:9),"Z"),] ## macro
      suffix="macro"
    } else {
      sweepmerged = sweepmerged[sweepmerged$chr %in% c("1","Z","LGE22"),] ## macro
      suffix="Z1LGE22"
      
    }
  }
  sweepmerged$color=as.numeric(as.factor(sweepmerged$chr))
  sweepmerged = unique(sweepmerged[!(is.na(sweepmerged$Fst)),])
  min_plot_per_chrom = aggregate(sweepmerged$plotorder~sweepmerged$chr,FUN=function(x){min(x,na.rm=T)}); colnames(min_plot_per_chrom)=c("chr","min")
  max_plot_per_chrom = aggregate(sweepmerged$plotorder~sweepmerged$chr,FUN=function(x){max(x,na.rm=T)}); colnames(max_plot_per_chrom)=c("chr","max")
  tot_plot_per_chrom = merge(min_plot_per_chrom,max_plot_per_chrom,all=F)
  sweepmerged$plotindex=sweepmerged$plotorder
  
  #print(tot_plot_per_chrom)
  for(rowindex in 2:nrow(tot_plot_per_chrom)){
    this_chr = tot_plot_per_chrom$chr[rowindex]
    #print(this_chr)
    num_to_subtract = (tot_plot_per_chrom$min[rowindex] - tot_plot_per_chrom$max[rowindex-1])-100
    #print(num_to_subtract)
    sweepmerged$plotindex[sweepmerged$chr==this_chr] = (sweepmerged$plotorder[sweepmerged$chr==this_chr])-num_to_subtract
    tot_plot_per_chrom[rowindex,2:3] = tot_plot_per_chrom[rowindex,2:3]-num_to_subtract
    
  }
  sweepmerged$plotorder=sweepmerged$plotindex
  min_plot_per_chrom = tot_plot_per_chrom[,c(1,2)]
  #print(tot_plot_per_chrom)
  specieslist = sort(unique(sweepmerged$species[!(is.na(sweepmerged$Fst))]))
  if(smallspp==T){
    specieslist = specieslist[specieslist %in% c("cur","bru","nit")]
    suffix=paste(suffix,"-small",sep="")
  }
  abline_locations =sweepmerged[which(sweepmerged$ranksppISLAND + sweepmerged$ranksppSWEEP > 0),"plotorder"]
  print(specieslist)
  png(paste("islswp_all_species_",suffix,"-2.png",sep=""),width=800,height=(150*length(specieslist)))
  palette(c("black","grey"))
  par(mfrow=c(length(specieslist),1),
      mar=c(2,4,0,0))
  for(i in 1:length(specieslist)){
    spp=specieslist[i]
    print(spp)
    small = sweepmerged[sweepmerged$species==spp,c("chr","species","Fst","midPos","plotorder","plotindex","color",
                                                   "ranksppISLAND","ranksppSWEEP")]
    small = unique(small)
    print(unique(small$chr))
    plot(sweepmerged$plotorder,sweepmerged$Fst,col=sweepmerged$color,xaxt="n",ylab=spp,type="n")
    abline(v=min_plot_per_chrom[,2],col=c("darkred","red"))
    abline(v=abline_locations,col="goldenrod")
    points(small$plotorder,small$Fst,col=small$color,xaxt="n",xaxt="n")
    points(small$plotorder[small$ranksppSWEEP==1],small$Fst[small$ranksppSWEEP==1],col="magenta",xaxt="n",xaxt="n",
           pch=16,cex=2)
    points(small$plotorder[small$ranksppISLAND==1],small$Fst[small$ranksppISLAND==1],col="cyan",xaxt="n",xaxt="n",
           pch=16,cex=3)
    axis(1,at=min_plot_per_chrom[,2],labels=min_plot_per_chrom[,1],las=2,cex=0.7)
  }
  dev.off()
}
## tajima's D same as above
{
  tajmerged = merged[merged$chr %in% unique(merged$chr[!(is.na(merged$Tajima))]),]
  if(micro==T){
    tajmerged = tajmerged[tajmerged$chr %in% c("1B",as.character(10:28),"mtDNA","LG2","LG5","LGE22"),] ## micro
    suffix="micro"
  } else {
    if(Z1==F){
      tajmerged = tajmerged[tajmerged$chr %in% c("1","1A",as.character(2:4),"4A",as.character(5:9),"Z"),] ## macro
      suffix="macro"
    } else {
      tajmerged = tajmerged[tajmerged$chr %in% c("1","Z","LGE22"),] ## macro
      suffix="Z1LGE22"
      
    }
  }
  tajmerged$color=as.numeric(as.factor(tajmerged$chr))
  tajmerged = unique(tajmerged[!(is.na(tajmerged$Tajima)),])
  min_plot_per_chrom = aggregate(tajmerged$plotorder~tajmerged$chr,FUN=function(x){min(x,na.rm=T)}); colnames(min_plot_per_chrom)=c("chr","min")
  max_plot_per_chrom = aggregate(tajmerged$plotorder~tajmerged$chr,FUN=function(x){max(x,na.rm=T)}); colnames(max_plot_per_chrom)=c("chr","max")
  tot_plot_per_chrom = merge(min_plot_per_chrom,max_plot_per_chrom,all=F)
  tajmerged$plotindex=tajmerged$plotorder
  #print(tot_plot_per_chrom)
  for(rowindex in 2:nrow(tot_plot_per_chrom)){
    this_chr = tot_plot_per_chrom$chr[rowindex]
    #print(this_chr)
    num_to_subtract = (tot_plot_per_chrom$min[rowindex] - tot_plot_per_chrom$max[rowindex-1])-100
    #print(num_to_subtract)
    tajmerged$plotindex[tajmerged$chr==this_chr] = (tajmerged$plotorder[tajmerged$chr==this_chr])-num_to_subtract
    tot_plot_per_chrom[rowindex,2:3] = tot_plot_per_chrom[rowindex,2:3]-num_to_subtract
    
  }
  tajmerged$plotorder=tajmerged$plotindex
  min_plot_per_chrom = tot_plot_per_chrom[,c(1,2)]
  #print(tot_plot_per_chrom)
  specieslist = sort(unique(tajmerged$species[!(is.na(tajmerged$Tajima))]))
  if(smallspp==T){
    specieslist = specieslist[specieslist %in% c("cur","bru","nit")]
    suffix=paste(suffix,"-small",sep="")
  }
  
  print(specieslist)
  png(paste("Tajima_all_species_",suffix,".png",sep=""),width=800,height=(150*length(specieslist)))
  palette(c("black","grey"))
  par(mfrow=c(length(specieslist),1),
      mar=c(2,4,0,0))
  for(i in 1:length(specieslist)){
    spp=specieslist[i]
    print(spp)
    small = tajmerged[tajmerged$species==spp,c("chr","species","Tajima","midPos","plotorder","plotindex","color")]
    small = unique(small)
    print(unique(small$chr))
    plot(tajmerged$plotorder,tajmerged$Tajima,col=tajmerged$color,xaxt="n",ylab=spp,type="n")
    abline(v=min_plot_per_chrom[,2],col=c("darkred","red"))
    points(small$plotorder,small$Tajima,col=small$color,xaxt="n",xaxt="n")
    axis(1,at=min_plot_per_chrom[,2],labels=min_plot_per_chrom[,1],las=2,cex=0.7)
  }
  dev.off()
}

## plot FST and DXY for cur/bil for presentation, black bg
{
  fstmerged = merged[merged$chr %in% unique(merged$chr[!(is.na(merged$Fst))]),]
  fstmerged$color=as.numeric(as.factor(fstmerged$chr))
  fstmerged = unique(fstmerged[!(is.na(fstmerged$Fst)),])
  min_plot_per_chrom = aggregate(fstmerged$plotorder~fstmerged$chr,FUN=function(x){min(x,na.rm=T)}); colnames(min_plot_per_chrom)=c("chr","min")
  max_plot_per_chrom = aggregate(fstmerged$plotorder~fstmerged$chr,FUN=function(x){max(x,na.rm=T)}); colnames(max_plot_per_chrom)=c("chr","max")
  tot_plot_per_chrom = merge(min_plot_per_chrom,max_plot_per_chrom,all=F)
  fstmerged$plotindex=fstmerged$plotorder
  for(rowindex in 2:nrow(tot_plot_per_chrom)){
    this_chr = tot_plot_per_chrom$chr[rowindex]
    num_to_subtract = (tot_plot_per_chrom$min[rowindex] - tot_plot_per_chrom$max[rowindex-1])-100
    fstmerged$plotindex[fstmerged$chr==this_chr] = (fstmerged$plotorder[fstmerged$chr==this_chr])-num_to_subtract
    tot_plot_per_chrom[rowindex,2:3] = tot_plot_per_chrom[rowindex,2:3]-num_to_subtract
  }
  fstmerged$plotorder=fstmerged$plotindex
  min_plot_per_chrom = tot_plot_per_chrom[,c(1,2)]
  specieslist = sort(unique(fstmerged$species[!(is.na(fstmerged$Fst))]))
  specieslist = specieslist[specieslist %in% c("cur","bil")]
  print(specieslist)
  fstmerged=fstmerged[fstmerged$species %in% specieslist,]
  
  png(paste("fst_dxy_cur_bil_presentation.png",sep=""),width=8,height=5,units="in",res=300)
  palette(c("white","grey"))
  par(mfrow=c(length(specieslist),2),
      mar=c(2,4,0,0))
  par(bg="black",col="white",col.axis="white",col.lab="white")
  
  m <- rbind(c(1,3), c(2,4))
  #print(m)
  layout(m)
  
  for(i in 1:length(specieslist)){
    spp=specieslist[i]
    print(spp)
    small = fstmerged[fstmerged$species==spp,c("chr","species","Fst","dxymeans","midPos","plotorder","plotindex","color")]
    small = unique(small)
    print(unique(small$chr))
    small = small[complete.cases(small[,c("Fst","dxymeans")]),]
    plot(fstmerged$Fst,col=fstmerged$color,xaxt="n",ylab=spp,type="n",xlim=c(0,nrow(small)))
    points(x=1:nrow(small),y=small$Fst,col=small$color,xaxt="n",xaxt="n")
    plot(fstmerged$dxymeans,col=fstmerged$color,xaxt="n",ylab=spp,type="n",xlim=c(0,nrow(small)))
    points(x=1:nrow(small),y=small$dxymeans,col=small$color,xaxt="n",xaxt="n")
  }
  
  dev.off()
}
{
  fstmerged = merged[merged$chr %in% unique(merged$chr[!(is.na(merged$Fst))]),]
  sweepmerged = fstmerged[fstmerged$chr %in% unique(fstmerged$chr[!(is.na(fstmerged$dxymeans))]),]
  sweepmerged$color=as.numeric(as.factor(sweepmerged$chr))
  sweepmerged = unique(sweepmerged[!(is.na(sweepmerged$Fst)),])
  min_plot_per_chrom = aggregate(sweepmerged$plotorder~sweepmerged$chr,FUN=function(x){min(x,na.rm=T)}); colnames(min_plot_per_chrom)=c("chr","min")
  max_plot_per_chrom = aggregate(sweepmerged$plotorder~sweepmerged$chr,FUN=function(x){max(x,na.rm=T)}); colnames(max_plot_per_chrom)=c("chr","max")
  tot_plot_per_chrom = merge(min_plot_per_chrom,max_plot_per_chrom,all=F)
  sweepmerged$plotindex=sweepmerged$plotorder
  for(rowindex in 2:nrow(tot_plot_per_chrom)){
    this_chr = tot_plot_per_chrom$chr[rowindex]
    #print(this_chr)
    num_to_subtract = (tot_plot_per_chrom$min[rowindex] - tot_plot_per_chrom$max[rowindex-1])-100
    #print(num_to_subtract)
    sweepmerged$plotindex[sweepmerged$chr==this_chr] = (sweepmerged$plotorder[sweepmerged$chr==this_chr])-num_to_subtract
    tot_plot_per_chrom[rowindex,2:3] = tot_plot_per_chrom[rowindex,2:3]-num_to_subtract
    
  }
  sweepmerged$plotorder=sweepmerged$plotindex
  min_plot_per_chrom = tot_plot_per_chrom[,c(1,2)]
  #print(tot_plot_per_chrom)
  specieslist = sort(unique(sweepmerged$species[!(is.na(sweepmerged$Fst))]))
  specieslist = specieslist[specieslist %in% c("cur")]
  abline_locations =sweepmerged[which(sweepmerged$ranksppISLAND + sweepmerged$ranksppSWEEP > 0),"plotorder"]
  print(specieslist)
  png(paste("islswp_cur_plotting.png",sep=""),width=10,height=4,units="in",res=300)
  palette(c("white","grey"))
  par(mfrow=c(length(specieslist),1),
      mar=c(2,4,2,0))
  par(bg="black",col="white",col.axis="white",col.lab="white")
  for(i in 1:length(specieslist)){
    spp=specieslist[i]
    print(spp)
    small = sweepmerged[sweepmerged$species==spp,c("chr","species","Fst","midPos","plotorder","plotindex","color",
                                                   "ranksppISLAND","ranksppSWEEP")]
    small = unique(small)
    print(unique(small$chr))
    plot(sweepmerged$plotorder,sweepmerged$Fst,col=sweepmerged$color,xaxt="n",ylab=spp,type="n")
    points(small$plotorder,small$Fst,col=small$color,xaxt="n",xaxt="n")
    points(small$plotorder[small$ranksppSWEEP==1],small$Fst[small$ranksppSWEEP==1],col="magenta",xaxt="n",xaxt="n",
           pch=16,cex=2)
    points(small$plotorder[small$ranksppISLAND==1],small$Fst[small$ranksppISLAND==1],col="cyan",xaxt="n",xaxt="n",
           pch=16,cex=3)
  }
  dev.off()
}


png("cur_dxy_fst_square.png")
par(bg="black",col="white",col.axis="white",col.lab="white")
plot(merged$dxymeans[merged$species=="cur"],
     merged$Fst[merged$species=="cur"])
abline(v=quantile(merged$dxymeans[merged$species=="cur"],c(0.05,0.95),na.rm=T),col="red")
abline(h=quantile(merged$Fst[merged$species=="cur"],c(0.05,0.95),na.rm=T),col="red")
points(merged$dxymeans[merged$species=="cur" & merged$ranksppSWEEP==1],
       merged$Fst[merged$species=="cur" & merged$ranksppSWEEP==1],
       col="magenta",pch=16,cex=2)
points(merged$dxymeans[merged$species=="cur" & merged$ranksppISLAND==1],
       merged$Fst[merged$species=="cur" & merged$ranksppISLAND==1],
       col="cyan",pch=16,cex=2)
dev.off()



## if you have
png("islands_vs_sweeps_relative.png",width=10,height=4,units="in",res=300)
par(bg="black",col="white",col.axis="white",col.lab="white")
x=t(table(merged$ranksppISLAND,merged$species))
islands=x[,2]
islands_rel = islands/sum(islands)
y=t(table(merged$ranksppSWEEP,merged$species))
sweeps=y[,2]
sweeps_rel = sweeps/sum(sweeps)
relative=cbind(islands_rel,sweeps_rel)
relative=relative[c("cur","bel","cri","fus",
                    "fla","bru","sin",
                    "mel","nit","bil"),]
par(mar=c(2,2,0,0))
b=barplot(t(relative),beside=T,col=c("cyan","magenta"),
          ylim=c(0,0.4),border="white")
box(col="white")
dev.off()


#merged$plotorder = factor(merged$plotorder, 
#                          levels = unique(as.character(merged$plotorder)))
#merged$plotorder = as.numeric(as.factor(merged$plotorder))
#write.table(merged,"rec_taj_dxy_fst_islandssweeps.temp")

#subset = merged[merged$chr %in% c(1,2,4,"4A","Z"),c("chr","midPos","species","Fst","dxymeans",
#                                                    "ranksppISLAND","ranksppSWEEP",
#                                                    "ranksppFST","plotorder")]
subset = merged[merged$chr %in% c(1,"1A",2,4,"4A",5,10,12:13,15,17:20,
                                  22:27,"Z"),
                c("chr","midPos","species","Fst","dxymeans",
                  "ranksppISLAND","ranksppSWEEP",
                  "ranksppFST","plotorder")]

subset=subset[complete.cases(subset),]
subset=unique(subset)
#only_islands = subset[subset$ranksppISLAND==1,]
#subset$plotorder = 1
subset = subset[order(subset$chr,subset$midPos,subset$species),]
dim(subset)
# totrows=nrow(subset)
# for(rowindex in 2:totrows){
#   if(rowindex %% 1000 == 0){print(paste(rowindex,"/",totrows))}
#   this_chr = subset$chr[rowindex]
#   prev_chr = subset$chr[rowindex-1]
#   this_pos = subset$midPos[rowindex]
#   prev_pos = subset$midPos[rowindex-1]
#   prev_ord = subset$plotorder[rowindex-1]
#   
#   if(this_chr==prev_chr && this_pos==prev_pos){
#     subset$plotorder[rowindex] = prev_ord
#   } else {
#     subset$plotorder[rowindex] = prev_ord+1
#   }
#   
# }

#only_sweeps = subset[subset$ranksppSWEEP==1,]
#names(merged)

# png("test_island.png",width=800,height=1500)
# palette(c("black","darkgrey","darkgreen","darkred","darkblue"))
# par(mfrow=c(10,1),
#     mar=c(4,4,0,0))
# for(spp in sort(unique(subset$species))){
#   plot(subset$Fst[subset$species==spp],col=as.numeric(as.factor(subset$chr[subset$species==spp])),
#        main=spp,cex=0.5)
#   points(subset$Fst[subset$species==spp],col=c(rgb(0,0,0,0),"cyan")[(subset$ranksppISLAND[subset$species==spp]+1)])
# }
# dev.off()

abline_locations =subset[which(subset$ranksppISLAND + subset$ranksppSWEEP > 0),c("chr","midPos","species","ranksppISLAND","ranksppSWEEP")]
dim(abline_locations)
dim(unique(abline_locations))
abline_locations=unique(abline_locations)
abline_locations$chrpos = paste(abline_locations$chr,abline_locations$midPos)

tab=table(abline_locations$chrpos)

overlapping=abline_locations[abline_locations$chrpos %in% names(tab[which(as.numeric(tab)>1)]),]
#write.table(overlapping,"overlaps.txt",row.names = F,sep="\t")

png("test_both_final.png",width=800,height=1500)
palette(c("black","darkblue"))
par(mfrow=c(10,1),
    mar=c(2,2,0,0))
abline_locations =subset$plotorder[which(subset$ranksppISLAND + subset$ranksppSWEEP > 0)]
for(spp in sort(unique(subset$species))){
  print(spp)
  plot(subset$plotorder[subset$species==spp],subset$Fst[subset$species==spp],
       col=as.numeric(as.factor(subset$chr[subset$species==spp])),
       main=spp,cex=0.25,type="n",xlim=c(min(subset$plotorder,na.rm=T),max(subset$plotorder,na.rm=T)))
  abline(v=abline_locations,col="grey")
  points(subset$plotorder[subset$species==spp],subset$Fst[subset$species==spp],
         col=as.numeric(as.factor(subset$chr[subset$species==spp])),
         main=spp,cex=0.25)
  points(subset$plotorder[subset$species==spp],subset$Fst[subset$species==spp],
         col=c(rgb(0,0,0,0),"magenta")[(subset$ranksppSWEEP[subset$species==spp]+1)],
         pch=16)
  points(subset$plotorder[subset$species==spp],subset$Fst[subset$species==spp],
         col=c(rgb(0,0,0,0),"cyan")[(subset$ranksppISLAND[subset$species==spp]+1)],
         pch=16)
  
}
dev.off()

pdf("test_both_final_eachchrom.pdf",width=8,height=15)
palette(c("black","darkblue"))
abline_locations =subset$plotorder[which(subset$ranksppISLAND + subset$ranksppSWEEP > 0)]
#png("chrom1_isl_sweeps_onlyfocal.png",width=6,height=4,units = "in",res=300)
png("chrom1_isl_sweeps.png",width=6,height=8,units = "in",res=300)
for(chrom in sort(unique(subset$chr))){
  print(chrom)
  thischrom = subset[subset$chr==chrom,]
  par(mfrow=c(length(unique(thischrom$species)),1),
      #par(mfrow=c(5,1),
      mar=c(2,3,0,0))
  for(spp in sort(unique(thischrom$species))){
    #for(spp in c("bel","bil","bru","cri","cur")){
    print(spp)
    plot(thischrom$plotorder[thischrom$species==spp],thischrom$Fst[thischrom$species==spp],
         col=as.numeric(as.factor(thischrom$chr[thischrom$species==spp])),yaxt="n",
         main="",cex=0.25,type="n",xlim=c(min(thischrom$plotorder,na.rm=T),max(thischrom$plotorder,na.rm=T)))
    axis(2,las=2)
    
    abline(v=abline_locations,col="grey")
    points(thischrom$plotorder[thischrom$species==spp],thischrom$Fst[thischrom$species==spp],
           col=as.numeric(as.factor(thischrom$chr[thischrom$species==spp])),
           main=spp,cex=0.25)
    points(thischrom$plotorder[thischrom$species==spp],thischrom$Fst[thischrom$species==spp],
           col=c(rgb(0,0,0,0),"magenta")[(thischrom$ranksppSWEEP[thischrom$species==spp]+1)],
           pch=16)
    points(thischrom$plotorder[thischrom$species==spp],thischrom$Fst[thischrom$species==spp],
           col=c(rgb(0,0,0,0),"cyan")[(thischrom$ranksppISLAND[thischrom$species==spp]+1)],
           pch=16)
    #legend("topleft",legend=paste(chrom,spp))
  }
}
dev.off()

## non-overlaps

non_overlap = data.frame()
for(spp in sort(unique(merged$species))){
  thisspp = merged[merged$species==spp,]
  #for(chrom in sort(unique(thisspp$chr))){
  thischrom=thisspp[thisspp$chr==chrom,]
  for(rowindex in seq(1,nrow(thischrom),10)){
    subset = thischrom[rowindex:(rowindex+9),]
    means=colMeans(subset[,c("plotorder","Fst","dxymeans","Tajima","weighted_recomb")],
                   na.rm = T)
    means=c(means,spp,chrom)
    names(means) = c(names(means)[1:5],"species","chrom")
    non_overlap = rbind(non_overlap,means)
  }
}
}




subset=subset[subset$chr %in% c(1,2,4,"4A","Z"),
              c("chr","midPos","species","Fst","dxymeans",
                "ranksppISLAND","ranksppSWEEP","ranksppFST","plotorder")]
png("test_both_islandzoom_final.png",width=800,height=1500)
palette(c("black","darkblue"))
par(mfrow=c(10,1),
    mar=c(2,2,0,0))
abline_locations =subset$plotorder[which(subset$ranksppISLAND + subset$ranksppSWEEP > 0)]
for(spp in sort(unique(subset$species))){
  print(spp)
  plot(subset$plotorder[subset$species==spp],subset$Fst[subset$species==spp],
       col=as.numeric(as.factor(subset$chr[subset$species==spp])),
       main=spp,cex=0.25,type="n",xlim=c(min(subset$plotorder,na.rm=T),max(subset$plotorder,na.rm=T)))
  abline(v=abline_locations,col="grey")
  points(subset$plotorder[subset$species==spp],subset$Fst[subset$species==spp],
         col=as.numeric(as.factor(subset$chr[subset$species==spp])),
         main=spp,cex=0.25)
  points(subset$plotorder[subset$species==spp],subset$Fst[subset$species==spp],
         col=c(rgb(0,0,0,0),"magenta")[(subset$ranksppSWEEP[subset$species==spp]+1)],
         pch=16)
  points(subset$plotorder[subset$species==spp],subset$Fst[subset$species==spp],
         col=c(rgb(0,0,0,0),"cyan")[(subset$ranksppISLAND[subset$species==spp]+1)],
         pch=16)
  
}
dev.off()






# png("test_both_bigger.png",width=800,height=1500)
# palette(c("black","darkgrey","darkgreen","darkred","darkblue"))
# par(mfrow=c(10,1),
#     mar=c(4,4,0,0))
# for(spp in sort(unique(subset2$species))){
#   plot(subset2$Fst[subset2$species==spp],col=as.numeric(as.factor(subset2$chr[subset2$species==spp])),
#        main=spp,cex=0.5)
#   points(subset2$Fst[subset2$species==spp],
#          col=c(rgb(0,0,0,0),"cyan")[(subset2$ranksppISLAND[subset2$species==spp]+1)],
#          pch=16)
#   points(subset2$Fst[subset2$species==spp],
#          col=c(rgb(0,0,0,0),"magenta")[(subset2$ranksppSWEEP[subset2$species==spp]+1)],
#          pch=16)
#   
# }
# dev.off()


# png("test_sweep.png",width=800,height=1500)
# palette(c("black","darkgrey","darkgreen","darkred","darkblue"))
# par(mfrow=c(10,1),
#     mar=c(4,4,0,0))
# for(spp in sort(unique(subset$species))){
#   plot(subset$Fst[subset$species==spp],col=as.numeric(as.factor(subset$chr[subset$species==spp])),
#        main=spp)
#   points(subset$Fst[subset$species==spp],col=c(rgb(0,0,0,0),"magenta")[(subset$ranksppSWEEP[subset$species==spp]+1)])
# }
# dev.off()

agg = aggregate(merged$windowstops ~ merged$chr,FUN=function(x){max(x,na.rm=T)})


cor(merged[,c("Tajima","weighted_recomb","Fst","dxymeans")],use="pairwise.complete.obs")

png("all_species_fst_dxy.png")
xval=merged$dxymeans[complete.cases(merged[,c("dxymeans","Fst")])]
yval=merged$Fst[complete.cases(merged[,c("dxymeans","Fst")])]
mod=lm(xval~yval)
mod_2 = lm(xval~poly(yval,2,raw=TRUE))
mod_3 = lm(xval~poly(yval,3,raw=TRUE))
mod_4 = lm(xval~poly(yval,4,raw=TRUE))
mod_e = lm(xval ~ exp(yval))
mod_l = lm(xval ~ log(yval))
mod_e2 = lm(exp(xval) ~ (yval))
mod_l2 = lm(log(xval) ~ (yval))
plot(yval,xval,type="p",col=rgb(0,0,0,0.05)#,
     #main=signif(summary(mod)$adj.r.squared,3)
)
xx <- seq(min(xval,na.rm=T),max(xval,na.rm=T), length=50)
yy <- seq(min(yval,na.rm=T),max(yval,na.rm=T), length=50)
#abline(mod,col="red")
#abline(mod_2,col="cyan")
lines(yy, predict(mod, list(yval=yy)), col="red")
lines(yy, predict(mod_2, list(yval=yy)), col="cyan")
lines(yy, predict(mod_3, list(yval=yy)), col="green")
lines(yy, predict(mod_4, list(yval=yy)), col="grey")
lines(yy, predict(mod_e, (list(yval=yy))), col="magenta")
lines(yy, predict(mod_l, (list(yval=yy))), col="lightgrey")
lines(yy, predict(mod_e2, (list(yval=yy))), col="blue") ## not vis
lines(yy, predict(mod_l2, (list(yval=yy))), col="darkgreen") ## not vis
dev.off()

png("all_species_taj_dxy.png")
xval=merged$Tajima[complete.cases(merged[,c("Tajima","dxymeans")])]
yval=merged$dxymeans[complete.cases(merged[,c("Tajima","dxymeans")])]
#mod=lm(xval~yval)
#mod_2 = lm(xval~poly(yval,2,raw=TRUE))
#mod_3 = lm(xval~poly(yval,3,raw=TRUE))
#mod_4 = lm(xval~poly(yval,4,raw=TRUE))
#mod_e = lm(xval ~ exp(yval))
mod_l = lm(xval ~ log(yval))
#mod_e2 = lm(exp(xval) ~ (yval))
#mod_l2 = lm(log(xval) ~ (yval))
plot(yval,xval,type="p",col=rgb(0,0,0,1),
     main=signif(summary(mod_l)$adj.r.squared,3),
     ylab="Tajima's D",xlab="DXY"
)
xx <- seq(min(xval,na.rm=T),max(xval,na.rm=T), length=50)
yy <- seq(min(yval,na.rm=T),max(yval,na.rm=T), length=50)
#abline(mod,col="red")
#abline(mod_2,col="cyan")
#lines(yy, predict(mod, list(yval=yy)), col="red",lwd=2,lty=2)
#lines(yy, predict(mod_2, list(yval=yy)), col="cyan",lwd=2,lty=2)
#lines(yy, predict(mod_3, list(yval=yy)), col="green",lwd=2,lty=2)
#lines(yy, predict(mod_4, list(yval=yy)), col="grey",lwd=2,lty=2)
#lines(yy, predict(mod_e, (list(yval=yy))), col="magenta",lwd=2,lty=2)
lines(yy, predict(mod_l, (list(yval=yy))), col="red",lwd=2,lty=2) ###
#lines(yy, predict(mod_e2, (list(yval=yy))), col="blue",lwd=2,lty=2)
#lines(yy, predict(mod_l2, (list(yval=yy))), col="darkgreen",lwd=2,lty=2)

dev.off()


only_chrs = unique(merged[,c("chr","midPos")])
t(table(only_chrs$chr))

# duprecomdf = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigrecom_duplicates.may2020.txt",
#                         sep=" ",header=T,stringsAsFactors = F,fill=T)
# duprecdf_spp = duprecomdf[duprecomdf$species==spp,]
# 
# windows = unique(rbind(windows,duprecdf_spp[,c("chr","midPos")]))
# windows <- windows[order(windows$chr, windows$midPos),]
# bothrecdf_spp = unique(rbind(recdf_spp,duprecdf_spp))
# 
# aggrecdf_spp = aggregate(bothrecdf_spp$weighted_recomb ~ bothrecdf_spp$chr + bothrecdf_spp$midPos + bothrecdf_spp$species,
#                          FUN=function(x){mean(x,na.rm=T)})

## looking at by species correlations
#merged=read.table("rec_taj_dxy_fst.temp",header=T,stringsAsFactors = F,fill=T)

mergx=merged[,c("Fst","plotorder","species","weighted_recomb","dxymeans","Tajima")]
mergx=unique(mergx)

fsttable=as.data.frame(as.matrix(with(merged, tapply(Fst, list(plotorder, species), FUN=function(x){mean(x,na.rm=T)}))))
rectable=as.data.frame(as.matrix(with(merged, tapply(weighted_recomb, list(plotorder, species), FUN=function(x){mean(x,na.rm=T)}))))
dxytable=as.data.frame(as.matrix(with(merged, tapply(dxymeans, list(plotorder, species), FUN=function(x){mean(x,na.rm=T)}))))
tajtable=as.data.frame(as.matrix(with(merged, tapply(Tajima, list(plotorder, species), FUN=function(x){mean(x,na.rm=T)}))))

neworder=c("bel","fla","nit","mel","bru","cri","cur","bil","fus","sin")
#numorder=c(1,6,9,8,3,4,5,2,7,10)
fsttable=fsttable[,neworder]
rectable=rectable[,neworder]
dxytable=dxytable[,neworder]
tajtable=tajtable[,neworder]


#dxysnptable=as.data.frame(as.matrix(with(dxydf, tapply(snps, list(chr, midPos, species), 
#                                                       FUN=function(x){mean(x,na.rm=T)}))))
#dxysnptable=dxysnptable[,neworder]


#pdf("corrplots_fst_rec_dxy_taj.pdf",width=6,height=6)
png("corrplots_fst_rec_dxy_taj.png",width=800,height=300)
#png("corrplots_fst_rec_dxy_color.png",width=600,height=300)
#par(mfrow=c(1,3),mar=c(0,0,0,0),cex=1)
par(mfrow=c(1,4),mar=c(0,0,0,0),cex=1)
corrplot::corrplot(cor(fsttable,use="pairwise.complete.obs"),
                   method="color",main="fst",
                   col=colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(200))
corrplot::corrplot(cor(rectable,use="pairwise.complete.obs"),
                   method="color",main="recomb rate",
                   col=colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(200))
corrplot::corrplot(cor(dxytable,use="pairwise.complete.obs"),
                   method="color",main="dxy",
                   col=colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(200))
corrplot::corrplot(cor(tajtable,use="pairwise.complete.obs"),
                   method="color",main="taj d",#order="hclust",
                   col=colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(200))
dev.off()

# merged = merged[merged$chr %in% c(1,2,4,5),]
# merged$plotorder = 1
# merged = merged[order(merged$chr,merged$midPos,merged$species),]
# totrows=nrow(merged)
# for(rowindex in 2:totrows){
#   if(rowindex %% 10 == 0){print(paste(rowindex,"/",totrows))}
#   this_chr = merged$chr[rowindex]
#   prev_chr = merged$chr[rowindex-1]
#   this_pos = merged$midPos[rowindex]
#   prev_pos = merged$midPos[rowindex-1]
#   prev_ord = merged$plotorder[rowindex-1]
#   
#   if(this_chr==prev_chr && this_pos==prev_pos){
#     merged$plotorder[rowindex] = prev_ord
#   } else {
#     merged$plotorder[rowindex] = prev_ord+1
#   }
#   
# }

png("dxttest.png",height=600,width=600)
palette(c("black","darkblue"))
par(mfrow=c(10,1),mar=c(2,2,0,0))
min_plot_per_chrom = aggregate(merged$plotorder~merged$chr,FUN=function(x){min(x,na.rm=T)})
max_plot_per_chrom = aggregate(merged$plotorder~merged$chr,FUN=function(x){max(x,na.rm=T)})

for(spp in sort(unique(merged$species))){
  plot(merged$plotorder[merged$species==spp],
       merged$dxymeans[merged$species==spp],
       col=as.numeric(as.factor(merged$chr[merged$species==spp])),
       xaxt="n",type="n")
  abline(v=min_plot_per_chrom[,2],col=c("darkgrey","lightgrey"))
  points(merged$plotorder[merged$species==spp],
         merged$dxymeans[merged$species==spp],
         col=as.numeric(as.factor(merged$chr[merged$species==spp])),
         xaxt="n")
  axis(1,at=min_plot_per_chrom[,2],labels=min_plot_per_chrom[,1])
}
dev.off()


merged$island50 = as.numeric(merged$ranksppDXY >= 0.5 & merged$ranksppFST >= 0.95)
merged$sweep50 = as.numeric(merged$ranksppDXY < 0.5 & merged$ranksppFST >= 0.95)
sum(merged$island50,na.rm=T)
sum(merged$sweep50,na.rm=T)
merged$color = (merged$island50-merged$sweep50)+2
merged$color[is.na(merged$color)] = 4

#smallmerged = unique(merged[!(is.na(merged$color)),])

png("fifty_sweep_test.png",width=1500,height=800)
par(mfrow=c(5,2),mar=c(0.1,0.1,0.1,0.1))
palette(c("magenta","black","cyan","grey"))
for(spp in sort(unique(merged$species))){
  plot(merged$plotorder[merged$species==spp],
       merged$Fst[merged$species==spp],
       col="grey",
       cex=0.25)
  points(merged$plotorder[merged$species==spp],
         merged$Fst[merged$species==spp],
         col=merged$color[merged$species==spp],
         cex=0.25)
}
dev.off()

write.table(merged,"islands_sweeps_fiftypercent.temp")
# merged=unique(merged)
# merged$rankFST = dplyr::percent_rank(merged$Fst)
# merged$rankDXY = dplyr::percent_rank(merged$dxymeans)
# merged$rankISLAND = as.numeric(merged$rankFST >= 0.95 & merged$rankDXY >= 0.95)
# merged$rankSWEEP = as.numeric(merged$rankFST >= 0.95 & merged$rankDXY <= 0.05)
# merged$ranksppFST = dplyr::percent_rank(merged$Fst)
# merged$ranksppDXY = dplyr::percent_rank(merged$dxymeans)
# for(spp in unique(merged$species)){
#   merged$ranksppFST[merged$species==spp] = dplyr::percent_rank(merged$Fst[merged$species==spp])
#   merged$ranksppDXY[merged$species==spp] = dplyr::percent_rank(merged$dxymeans[merged$species==spp])
# }
# merged$ranksppISLAND = as.numeric(merged$ranksppFST >= 0.95 & merged$ranksppDXY >= 0.95)
# merged$ranksppSWEEP = as.numeric(merged$ranksppFST >= 0.95 & merged$ranksppDXY <= 0.05)
# merged$sppchr = paste(merged$chr,merged$species,sep="-")
# png("test_both.png",width=800,height=1500)
# palette(c("black","darkblue"))
# par(mfrow=c(10,1),
#     mar=c(2,2,0,0))
# abline_locations =merged$plotorder[which(merged$ranksppISLAND + merged$ranksppSWEEP > 0)]
# for(spp in sort(unique(merged$species))){
#   plot(merged$plotorder[merged$species==spp],merged$Fst[merged$species==spp],col=as.numeric(as.factor(merged$chr[merged$species==spp])),
#        main=spp,cex=0.25,type="n")
#   abline(v=abline_locations,col="grey")
#   points(merged$plotorder[merged$species==spp],merged$Fst[merged$species==spp],col=as.numeric(as.factor(merged$chr[merged$species==spp])),
#          main=spp,cex=0.25)
#   points(merged$plotorder[merged$species==spp],merged$Fst[merged$species==spp],
#          col=c(rgb(0,0,0,0),"cyan")[(merged$ranksppISLAND[merged$species==spp]+1)],
#          pch=16)
#   points(merged$plotorder[merged$species==spp],merged$Fst[merged$species==spp],
#          col=c(rgb(0,0,0,0),"magenta")[(merged$ranksppSWEEP[merged$species==spp]+1)],
#          pch=16)
#   
# }
# dev.off()
# 
# write.table(merged,"merged_merged.txt",append=T)

dxyfst = read.table("fst_dxy.temp",header=T,stringsAsFactors = F)
for(spp in unique(dxyfst$species)){
  dxyfst$ranksppFST[dxyfst$species==spp] = dplyr::percent_rank(dxyfst$Fst[dxyfst$species==spp])
  dxyfst$ranksppDXY[dxyfst$species==spp] = dplyr::percent_rank(dxyfst$means[dxyfst$species==spp])
}
dxyfst$ranksppSWEEP = as.numeric(dxyfst$ranksppDXY < 0.05 & dxyfst$ranksppFST >= 0.95)
dxyfst$ranksppISLAND = as.numeric(dxyfst$ranksppDXY >= 0.95 & dxyfst$ranksppFST >= 0.95)

dxyfst=unique(dxyfst)
test2= dxyfst[dxyfst$ranksppISLAND==1,]
test2=test2[,c("chr","midPos","species")]
test2=unique(test2)
table(test2$species)

test = merged[merged$ranksppISLAND==1,c("chr","midPos","species")]
test=unique(test)
table(test$species)

mergx=merged[,c("chr","midPos","species","ranksppDXY","ranksppFST")]
mergx=unique(mergx)
mergx$ranksppSWEEP = as.numeric(mergx$ranksppDXY < 0.05 & mergx$ranksppFST >= 0.95)
mergx$ranksppISLAND = as.numeric(mergx$ranksppDXY >= 0.95 & mergx$ranksppFST >= 0.95)
mergx=unique(mergx)
mergx=unique(mergx[(mergx$ranksppSWEEP==1 | mergx$ranksppISLAND == 1),])
dim(mergx)

#onlyswp = mergx[mergx$ranksppSWEEP==1,c("chr","midPos","species")]
#onlyswp=unique(onlyswp)

fst0.95_dxysweep=c()
fst0.95_dxyisland=c()

island_cutoff = matrix(nrow=21,ncol=21)
for(i in 1:nrow(island_cutoff)){
  for(j in 1:ncol(island_cutoff)){
    dxy_i = seq(0,1,0.05)[i]
    fst_j = seq(0,1,0.05)[j]
    island_cutoff[i,j] = (sum(mergx$ranksppDXY >= dxy_i & mergx$ranksppFST >= fst_j,na.rm=T))
    print(paste(i,"dxy:",dxy_i,"/",j,"fst:",fst_j,"/ cutoff:",island_cutoff[i,j]))
  }
}
#island_cutoff = island_cutoff/max(island_cutoff,na.rm=T)
rownames(island_cutoff)=paste("DXY",seq(0,1,0.05),sep="-")
colnames(island_cutoff)=paste("FST",seq(0,1,0.05),sep="-")
cutoff_raster = raster(island_cutoff)

## not varying dxy_i for this
sweep_cutoff = matrix(nrow=21,ncol=21)
for(i in 1:nrow(sweep_cutoff)){
  for(j in 1:ncol(sweep_cutoff)){
    #dxy_i = seq(0,1,0.05)[i]
    dxy_i = 0.05
    fst_j = seq(0,1,0.05)[j]
    sweep_cutoff[i,j] = (sum(mergx$ranksppDXY < dxy_i & mergx$ranksppFST >= fst_j,na.rm=T))
    print(paste(i,"dxy:",dxy_i,"/",j,"fst:",fst_j,"/ cutoff:",sweep_cutoff[i,j]))
  }
}
#sweep_cutoff = sweep_cutoff/max(sweep_cutoff,na.rm=T)
#rownames(sweep_cutoff)=paste("DXY",seq(0,1,0.05),sep="-")
rownames(sweep_cutoff)=paste("DXY",rep(0.05,21),sep="-")
colnames(sweep_cutoff)=paste("FST",seq(0,1,0.05),sep="-")
cutoff_raster_sweep = raster(sweep_cutoff)
#plot(cutoff_raster_sweep,col=pal(200),las=2,axes=F,interpolate=F,
#     xlab="FST percentile",ylab="DXY percentile")

fst_upper = 0.05
for(fst_upper in seq(0,0.95,0.05)){
  print(fst_upper)
  dxy_upper = (seq(0,1,0.05))
  dxy_lower = rev(seq(0,1,0.05))
  dxy_dxy_matrix = (matrix(nrow=21,ncol=21))
  for (i in 1:length(dxy_lower)){
    for(j in 1:length(dxy_upper)){
      
      this_dxy_lower = dxy_lower[i]
      this_dxy_upper = dxy_upper[j]
      if(this_dxy_lower<=this_dxy_upper){
        sweeps_expected = sum((merged$ranksppDXY <= this_dxy_lower) & (merged$ranksppFST >= fst_upper),na.rm=T)
        islands_expected = sum((merged$ranksppDXY >= this_dxy_upper) & (merged$ranksppFST >= fst_upper),na.rm=T)
        #expected_ratio = ((sweeps_expected-islands_expected)/(sweeps_expected+islands_expected))
        
        if(sweeps_expected > islands_expected){
          expected_ratio=log10(sweeps_expected/islands_expected)
          if(islands_expected==0){expected_ratio=log10(sweeps_expected)}
        } else if(islands_expected > sweeps_expected) {
          expected_ratio=-log10(islands_expected/sweeps_expected)
          if(sweeps_expected==0){expected_ratio=-log10(islands_expected)}
        } else {
          if(islands_expected!=0){expected_ratio=0} else {expected_ratio=NA}
        }
        
        # expected_ratio = (sweeps_expected-islands_expected)/islands_expected
        # if(islands_expected==0){
        #   expected_ratio = 1000000
        #   if(sweeps_expected==0){
        #     expected_ratio=NA
        #   }
        # }
        
        dxy_dxy_matrix[i,j] = expected_ratio
        #print(paste("lower",this_dxy_lower,"upper",this_dxy_upper,"sweeps",
        #            sweeps_expected,"islands",islands_expected,"ratio",signif(expected_ratio,2)))
        #rownames(dxy_dxy_matrix)[i] = as.character(this_dxy_lower)
        #colnames(dxy_dxy_matrix)[j] = as.character(this_dxy_upper)
      } else {
        dxy_dxy_matrix[i,j] = NA
      }
    }
  }
  dxy_dxy_raster = raster(dxy_dxy_matrix)
  print(summary(dxy_dxy_raster))
  #log_raster = log(abs(dxy_dxy_raster))
  #log_raster = log_raster[dxy_dxy_raster]
  png(paste("expected_sweep_to_island_ratio_fst",fst_upper,".png",sep=""),width=400,height=400)
  pal = colorRampPalette(c("white","cyan"))
  pal2 = colorRampPalette(c("white","magenta"))
  pal3 = colorRampPalette(c("darkcyan","cyan","white","magenta","darkmagenta"))
  par(mar=c(4,4,0,0))
  #cuts=c(-1,-0.5,-0.01,0.01,0.5,1,10,100,1000,10000,100000)
  #cuts=seq(-10.5,10.5,1)
  cuts=c(-6,-2,-1,-0.0001,0.0001,1,2,6)
  plot(dxy_dxy_raster,xlab="Island Cutoff",ylab="Sweep Cutoff",
       breaks=cuts,
       col=pal3(7),
       #col=c(rev(pal(8))[5:8],pal2(9)),asp=NA,
       #col=c(rev(pal(3)),pal2(8)),
       las=2,axes=F,interpolate=F,colNA = "grey",
       main="",legend=F,zlim=c(-6,6))
  legend(x=0,y=1,
         #fill=c(rev(pal(8))[5:8],pal2(8)),
         fill=pal3(7),
         title="Islands:Sweeps",
         legend=c("Over 100x I","100x I","10x I","Equal",
                  "10x S","100x S","Over 100x S"),
         bty="n")
  #text(dxy_dxy_raster,cex=0.5,digits=2)
  at_points = seq(0,1,length.out=43)[seq(2,43,2)]
  axis(side = 1, at = at_points,las=2,labels=seq(0,1,0.05)*100)
  axis(side = 2, at = at_points,las=2,labels=seq(0,1,0.05)*100)
  polygon(x=c(0.95,0.9,0.9,0.95),
          y=c(0.05,0.05,0.1,0.1))
  text(x=at_points[9],y=at_points[20],labels=paste("FST=",fst_upper),
       adj=c(0,1),cex=2)
  dev.off()
}

png("fst_vs_dxy_sweep_proportion.png",width=400,height=400)
par(mar=c(4,4,0,0))
#par(mfrow=c(1,3),mar=c(4,4,4,0))
# plot(cutoff_raster,col=pal(200),las=2,axes=F,interpolate=F,
#      xlab="",ylab="DXY percentile",main="Percentage of Genome as Islands")
# #polygon(x=c(0,1,1,0),y=c(0.05,0.05,0.10,0.10))
# polygon(x=c(0,1,1,0),y=c(0.9,0.9,0.95,0.95))
# polygon(y=c(0,1,1,0),x=c(0.9,0.9,0.95,0.95))
# at_points = seq(0,1,length.out=43)[seq(2,43,2)]
# axis(side = 1, at = at_points,las=2,labels=seq(0,1,0.05)*100)
# axis(side = 2, at = at_points,las=2,labels=seq(0,1,0.05)*100)
# plot(cutoff_raster_sweep,col=pal2(200),las=2,axes=F,interpolate=F,
#      xlab="FST percentile",ylab="",main="Percentage of Genome as Sweeps")
# polygon(x=c(0,1,1,0),y=c(0.05,0.05,0.10,0.10))
# #polygon(x=c(0,1,1,0),y=c(0.9,0.9,0.95,0.95))
# polygon(y=c(0,1,1,0),x=c(0.9,0.9,0.95,0.95))
# at_points = seq(0,1,length.out=43)[seq(2,43,2)]
# axis(side = 1, at = at_points,las=2,labels=seq(0,1,0.05)*100)
# axis(side = 2, at = at_points,las=2,labels=seq(0,1,0.05)*100)
plot(((cutoff_raster_sweep - cutoff_raster) / (cutoff_raster_sweep + cutoff_raster)),
     col=pal3(200),las=2,axes=F,interpolate=F,
     xlab="Upper DXY cutoff",ylab="FST cutoff",
     main="")
polygon(x=c(0,1,1,0),y=c(0.05,0.05,0.10,0.10))
polygon(x=c(0,1,1,0),y=c(0.9,0.9,0.95,0.95))
polygon(y=c(0,1,1,0),x=c(0.9,0.9,0.95,0.95))
dev.off()


island_cutoff_empirical = island_cutoff[c(20),]
sweep_cutoff_empirical = sweep_cutoff[c(2),]
barplot(rbind(island_cutoff_empirical,sweep_cutoff_empirical))

barplot(rbind(island_cutoff_empirical,sweep_cutoff_empirical),las=2,
        col=c("cyan","magenta"),yaxt="n",
        names=seq(0,1,0.05))
axis(2,at=seq(0,13938,length.out = 5),labels=seq(0,1,length.out=5),las=2)



for(value in seq(0,1,0.05)){
  fst0.95_dxyisland =  c(fst0.95_dxyisland,(sum(mergx$ranksppDXY >= value & mergx$ranksppFST >= 0.95,na.rm=T)))
  fst0.95_dxysweep =  c(fst0.95_dxysweep,(sum(mergx$ranksppDXY < value & mergx$ranksppFST >= 0.95,na.rm=T)))
  
}

plot(cbind(seq(0,1,0.05),fst0.95_dxysweep),type="b")
points(cbind(seq(0,1,0.05),fst0.95_dxyisland),type="b",col="red",lty=2,pch=16)

barplot(rbind(fst0.95_dxysweep,fst0.95_dxyisland),col=c("magenta","cyan"),
        names=seq(0,1,0.05),las=2,ylab="Total Number",xlab=("DXY cutoff"))
box()


##
agg1 = aggregate(merged$ranksppISLAND~merged$chr+merged$species,FUN=function(x){sum(x==1,na.rm=T)})
agg2 = aggregate(merged$ranksppISLAND~merged$chr+merged$species,FUN=function(x){sum(x==0,na.rm=T)})
agg3 = aggregate(merged$Fst~merged$chr+merged$species,FUN=function(x){mean(x,na.rm=T)})
colnames(agg1)=c("chr","species","yesfst")
colnames(agg2)=c("chr","species","nofst")
colnames(agg3)=c("chr","species","fst")
aggm = merge(agg1,agg2,all=T)
aggm = merge(aggm,agg3,all=T)
#aggm=aggm[aggm$yesfst<=10,]
plot(aggm$fst,aggm$yesfst)
plot(aggm$fst,aggm$nofst)
plot(aggm$fst,aggm$yesfst/(aggm$nofst+aggm$yesfst))
#plot(aggm$yesfst+aggm$nofst,aggm$yesfst)
summary(lm((aggm$yesfst/(aggm$nofst+aggm$yesfst))~aggm$fst))
aggm$species[aggm$species=="cur"]=1
aggm$species[aggm$species=="bel"]=2
aggm$species[aggm$species=="cri"]=3
aggm$species[aggm$species=="fus"]=4
aggm$species[aggm$species=="fla"]=5
aggm$species[aggm$species=="bru"]=6
aggm$species[aggm$species=="sin"]=7
aggm$species[aggm$species=="mel"]=8
aggm$species[aggm$species=="nit"]=9
aggm$species[aggm$species=="bil"]=10
plot(aggm$species,aggm$fst)
plot(aggm$species,aggm$yesfst)
plot(aggm$species,aggm$yesfst/(aggm$yesfst+aggm$nofst))
summary(lm(aggm$yesfst~aggm$species))



png("FST_DXY_peaks_troughs_comparison_hists.png")
mergx=merged[merged$ranksppFST>=0.95,]
#png("FST_peaks_by_dxy_values.png",height=300)
#par(mfrow=c(1,2),mar=c(4,4,1,0.1))
par(mfrow=c(3,2),mar=c(4,4,1,0.1))
h = hist(mergx$dxymeans,plot=F,breaks=20) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE,ylab="Percentage of FST Peaks",xlab="DXY value",main="",
     ylim=c(0,100))
cumulative = sapply(1:length(h$density),FUN=function(x){sum(h$density[1:x])})
points(c(0,h$mids),c(0,cumulative),col="red",type="l")

h = hist(mergx$ranksppDXY,plot=F,breaks=20) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE,ylab="Percentage of FST Peaks",xlab="DXY rank",main="",
     ylim=c(0,100),col=c("white","grey",rep("white",17),"grey"))
cumulative = sapply(1:length(h$density),FUN=function(x){sum(h$density[1:x])})
points(c(0,h$mids),c(0,cumulative),col="red",type="l")
#dev.off()

mergx=merged[merged$ranksppDXY>=0.95,]
#png("DXY_peaks_by_fst_values.png",height=300)
#par(mfrow=c(1,2),mar=c(4,4,1,0.1))
h = hist(mergx$Fst,plot=F,breaks=20) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE,ylab="Percentage of DXY Peaks",xlab="FST value",main="",
     ylim=c(0,100))
cumulative = sapply(1:length(h$density),FUN=function(x){sum(h$density[1:x])})
points(c(0,h$mids),c(0,cumulative),col="red",type="l")

h = hist(mergx$ranksppFST,plot=F,breaks=20) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE,ylab="Percentage of DXY Peaks",xlab="FST rank",main="",
     ylim=c(0,100),col=c(rep("white",19),"grey"))
cumulative = sapply(1:length(h$density),FUN=function(x){sum(h$density[1:x])})
points(c(0,h$mids),c(0,cumulative),col="red",type="l")
#dev.off()

mergx=merged[merged$ranksppDXY<0.05,]
#png("DXY_troughs_by_fst_values.png",height=300)
#par(mfrow=c(1,2),mar=c(4,4,1,0.1))
h = hist(mergx$Fst,plot=F,breaks=20) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE,ylab="Percentage of DXY Troughs",xlab="FST value",main="",
     ylim=c(0,100))
cumulative = sapply(1:length(h$density),FUN=function(x){sum(h$density[1:x])})
points(c(0,h$mids),c(0,cumulative),col="red",type="l")

h = hist(mergx$ranksppFST,plot=F,breaks=20) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)*100
plot(h,freq=FALSE,ylab="Percentage of DXY Troughs",xlab="FST rank",main="",
     ylim=c(0,100),col=c(rep("white",19),"grey"))
cumulative = sapply(1:length(h$density),FUN=function(x){sum(h$density[1:x])})
points(c(0,h$mids),c(0,cumulative),col="red",type="l")
dev.off()

















# MISSING DATA


only_miss = merged[!(is.na(merged$mean_F_MISS)),c("chr","midPos","species","weighted_recomb",
                                                  "Tajima","Fst","dxymeans","ranksppISLAND",
                                                  "ranksppSWEEP","plotorder","mean_F_MISS",
                                                  "ranksppDXY")]
only_miss_fst = only_miss[!(is.na(only_miss$Fst)),]
only_miss_dxy = only_miss[!(is.na(only_miss$dxymeans)),]
only_miss_taj = only_miss[!(is.na(only_miss$Tajima)),]
only_miss_rec = only_miss[!(is.na(only_miss$weighted_recomb)),]
only_miss_swp = only_miss[!(is.na(only_miss$ranksppSWEEP)),]
only_miss_isl = only_miss[!(is.na(only_miss$ranksppISLAND)),]

only_miss_dxyrank = only_miss[!(is.na(only_miss$ranksppDXY)),]
only_miss_dxyrank$ranksppMISS = NA
for(spp in unique(only_miss_dxyrank$species)){
  only_miss_dxyrank$ranksppMISS[only_miss_dxyrank$species==spp] = dplyr::percent_rank(as.numeric(only_miss_dxyrank$mean_F_MISS[only_miss_dxyrank$species==spp]))
}

only_miss_dxyrank$lowdxy = only_miss_dxyrank$ranksppDXY<0.05
only_miss_dxyrank$highdxy = only_miss_dxyrank$ranksppDXY>=0.95

boxplot(as.numeric(only_miss_dxyrank$mean_F_MISS)~as.character(only_miss_dxyrank$lowdxy)
        +as.character(only_miss_dxyrank$highdxy),
        names=c("middle","low","high","n/a"))


plot(only_miss_dxyrank$ranksppMISS,only_miss_dxyrank$ranksppDXY)
mod=lm(as.numeric(only_miss_dxyrank$ranksppDXY)~as.numeric(only_miss_dxyrank$ranksppMISS))
abline(mod,col="red")
summary(mod)

cor(cbind(fst=as.numeric(only_miss$Fst),dxy=as.numeric(only_miss$dxymeans),
          taj=as.numeric(only_miss$Tajima),rec=as.numeric(only_miss$weighted_recomb),
          miss=as.numeric(only_miss$mean_F_MISS)),
    use="pairwise.complete.obs")

## dxy and missing are fairly correlated
plot(as.numeric(only_miss_dxy$mean_F_MISS),as.numeric(only_miss_dxy$dxymeans),
     col=rgb(0,0,0,0.1))
mod=lm(as.numeric(only_miss_dxy$dxymeans)~as.numeric(only_miss_dxy$mean_F_MISS))
abline(mod,col="red")
summary(mod)

nooutliers = only_miss_dxy[dplyr::percent_rank(as.numeric(only_miss_dxy$dxymeans))<0.95,]
plot(as.numeric(nooutliers$mean_F_MISS),as.numeric(nooutliers$dxymeans),
     col=rgb(0,0,0,0.1))
mod=lm(as.numeric(nooutliers$dxymeans)~as.numeric(nooutliers$mean_F_MISS))
abline(mod,col="red")
summary(mod)

plot(dplyr::percent_rank(as.numeric(nooutliers$dxymeans)),
     as.numeric(nooutliers$dxymeans))
abline(v=0.95)


nooutliers = nooutliers[dplyr::percent_rank(as.numeric(nooutliers$mean_F_MISS))>0.05,]
plot(as.numeric(nooutliers$mean_F_MISS),as.numeric(nooutliers$dxymeans),
     col=rgb(0,0,0,0.1))
mod=lm(as.numeric(nooutliers$dxymeans)~as.numeric(nooutliers$mean_F_MISS))
abline(mod,col="red")
summary(mod)



nooutliers = only_miss_dxy[dplyr::percent_rank(as.numeric(only_miss_dxy$mean_F_MISS))>0.05,]
plot(as.numeric(nooutliers$mean_F_MISS),as.numeric(nooutliers$dxymeans),
     col=rgb(0,0,0,0.1))
mod=lm(as.numeric(nooutliers$dxymeans)~as.numeric(nooutliers$mean_F_MISS))
abline(mod,col="red")
summary(mod)

plot(dplyr::percent_rank(as.numeric(nooutliers$mean_F_MISS)),
     as.numeric(nooutliers$mean_F_MISS))
abline(v=0.05)


png("~/testmissing.png",width=800,height=500)
par(mfrow=c(2,3))
plot(only_miss_fst$mean_F_MISS,only_miss_fst$Fst,col=as.numeric(as.factor(only_miss_fst$species)))
plot(only_miss_dxy$mean_F_MISS,only_miss_dxy$dxymeans,col=as.numeric(as.factor(only_miss_dxy$species)))
plot(only_miss_taj$mean_F_MISS,only_miss_taj$Tajima,col=as.numeric(as.factor(only_miss_taj$species)))
plot(only_miss_rec$mean_F_MISS,only_miss_rec$weighted_recomb,col=as.numeric(as.factor(only_miss_rec$species)))
boxplot(as.numeric(only_miss_swp$mean_F_MISS)~as.character(only_miss_swp$ranksppSWEEP),
        outline=T,main="no outliers")
boxplot(as.numeric(only_miss_isl$mean_F_MISS)~as.character(only_miss_isl$ranksppISLAND),
        outline=T,main="no outliers")
#mod = lm(only_miss$Fst~only_miss$mean_F_MISS)
dev.off()


#merged=unique(merged)
#merged$type="none"
#merged$type[merged$ranksppISLAND==1] = "island"
#merged$type[merged$ranksppSWEEP==1] = "sweep"

boxplot(merged$mean_F_MISS~merged$type)
