dxyfiles = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY",
                      pattern="_4_SON_Dxy_WINDOWS_chrfix_1-ALL.txt$",full.names = T,recursive = T)

sin=read.table(dxyfiles[10],header=T,stringsAsFactors = F)
bel=read.table(dxyfiles[1],header=T,stringsAsFactors = F)

bel_1m = bel[bel$scafs %in% unique(sin$scafs),]

both=rbind(sin,bel)

notdup = data.frame()
for(chrom in unique(sin$scafs)){
  print(chrom)
  for(position in unique(sin$starts[sin$scafs==chrom])){
    subset = sin[(sin$starts==position & sin$scafs==chrom),]
    if (nrow(subset) > 1) {
      tocheck = bel[(bel$starts==position & bel$scafs==chrom),]
      subset=subset[!(subset$means %in% tocheck$means),]
    }
    notdup=rbind(notdup,subset)
  }
}
notdup=unique(notdup)

write.table(notdup,paste(dxyfiles[10],"fix.txt",sep="."),
            row.names = F,sep="\t",quote = F)



merged=read.table("~/rec_taj_dxy_fst_islandssweeps.temp",header=T,fill=T,
                  stringsAsFactors = F)
merged$dxymeans[merged$species=="sin"] = NA
dxyfstsin = dxyfst[dxyfst$species=="sin",]
colnames(dxyfstsin)[5] = "dxymeans"
dxyfstsin = dxyfstsin[!(is.na(dxyfstsin$dxymeans)),]
head(dxyfstsin)
for(rowindex in 1:nrow(dxyfstsin)){
  if(rowindex %% 100 == 0){print(rowindex)}
  thisrow=dxyfstsin[rowindex,]
  chrom = thisrow$chr
  midpos = thisrow$midPos
  fst = thisrow$Fst
  merged$dxymeans[(merged$species=="sin" & merged$chr==chrom & merged$midPos==midpos)] = thisrow$dxymeans
}

png("belsin.png")
plot(merged$dxymeans[merged$species=="bel"],
     merged$dxymeans[merged$species=="sin"])
dev.off()

merged$ranksppDXY[merged$species=="sin"] = NA
merged$ranksppDXY[merged$species=="sin"] = dplyr::percent_rank(merged$dxymeans[merged$species=="sin"])
merged$ranksppISLAND[merged$species=="sin"] = as.numeric(merged$ranksppFST[merged$species=="sin"] >= 0.95 
                                                         & merged$ranksppDXY[merged$species=="sin"] >= 0.95)
merged$ranksppSWEEP[merged$species=="sin"] = as.numeric(merged$ranksppFST[merged$species=="sin"] >= 0.95 
                                                        & merged$ranksppDXY[merged$species=="sin"] < 0.05)

merged$rankDXY = dplyr::percent_rank(merged$dxymeans)
merged$rankISLAND = as.numeric(merged$rankFST >= 0.95 & merged$rankDXY >= 0.95)
merged$rankSWEEP = as.numeric(merged$rankFST >= 0.95 & merged$rankDXY < 0.05)

write.table(merged,"rec_taj_dxy_fst_islandssweeps.temp")
