merged=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore_2.temp",
                  fill=T,header=T,stringsAsFactors = F,strip.white = T,comment.char = "")

fix_dups=F
if(fix_dups==T) {

dups = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore_2.FIXED.txt",header=T,fill=T,stringsAsFactors = F,strip.white = T,comment.char = "")
not_dups = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore_2.NOTDUPS.txt",header=T,fill=T,stringsAsFactors = F,strip.white = T,comment.char = "")

merged$color[merged$color=="empty"] = NA
merged$sppchr = paste(merged$chr,merged$species,sep="-")
merged$sppchrpos = paste(merged$chr,merged$species,as.character(merged$midPos),sep="-")
## MUST fix the duplicates issue

#original = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore.temp",
#                      fill=T,header=T,stringsAsFactors = F,strip.white = T)
#kept = unique(merged$sppchrpos)
#original$sppchrpos = paste(original$chr,original$species,as.character(original$midPos),sep="-")
#notkept = original[!(original$sppchrpos %in% kept),]
#dim(notkept)
#merged = rbind(merged,notkept)

posnums=rev(table(merged$sppchrpos))
names = unique(names(posnums)[posnums>1])

dups = merged[merged$sppchrpos %in% names,]
dups = unique(dups)
not_dups = merged[!(merged$sppchrpos %in% names),]
not_dups = unique(not_dups)

write.table(not_dups,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore_2.NOTDUPS.temp",row.names = F,quote = F)
write.table(dups,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore_2.FIXED.temp",row.names = F,quote = F)

chromlist=c("LG2","1B","LGE22",25,22,28:26,8,24,23,21,
            13:11,20,10,19,17,18,15,14,"4A",9,6,7,4,5,"Z",3,1,2)
spplist=c("mel","fus","cur","bru","cri","nit","sin","bil","bel")

## write a function that does what I want

# numeric = sapply(merged,varhandle::check.numeric)
# colSums(numeric)==nrow(merged)
# summary(as.numeric(x$zscoresppFST[!(is.na(x$Fst))]))
# x=merged[order(merged$species,merged$chr,merged$midPos),]
# match=sapply(2:nrow(x),FUN=function(i){
#   if(x$chr[i]==x$chr[i-1] & x$species[i]==x$species[i-1] & x$outlier[i]==x$outlier[i-1]) {
#     return(0)
#   } else {
#     return(1)
#   }
# })
# x$match = c(1,match)
# write.table(x,"fstoutlier.small.csv",quote = F,sep=",")

fix_column_dups = function(df,columns=1:ncol(df)) {
  classes = sapply(df,class)
  numeric = sapply(df,varhandle::check.numeric)
  classes = unique(union(which(colSums(numeric==FALSE)==0),which(classes=="numeric")))
  
  tocolsum=intersect(columns,intersect(which(colSums(!(is.na(df)))==1),classes))
  df[,tocolsum] <- apply(df[,tocolsum], 2, as.numeric)
  colsums = rbind(colSums(df[,tocolsum],na.rm=T))
  toreplace = do.call("rbind", replicate(nrow(df), colsums, simplify = FALSE))
  df[,tocolsum] = toreplace
  
  uniquevalcols = (intersect(which(colSums(!(is.na(df)))==1),which(!(columns %in% classes))))
  uniques=sapply(uniquevalcols,FUN=function(x){
    return(df[!(is.na(df[,x])),x])
  })
  toreplace2 = do.call("rbind", replicate(nrow(df), uniques, simplify = FALSE))
  df[,uniquevalcols] = toreplace2
  
  are_blank = which(colSums((is.na(df)))>1)
  
  leftover = which(!(columns %in% union(classes,union(uniquevalcols,are_blank))))
  
  pasted_values=sapply(leftover,FUN=function(x){
    u = unique(df[,x])
    u = u[!(is.na(u))]
    return(paste(u,collapse=";"))
  })
  toreplace3 = do.call("rbind", replicate(nrow(df), pasted_values, simplify = FALSE))
  df[,leftover] = toreplace3
  
  return(df)
}

for (chr in chromlist) {
  print(chr)
  df_c = dups[dups$chr==chr,]
  for(spp in spplist) {
    print(spp)
    #df_s = dups[dups$species==spp,]
    df_cs = df_c[df_c$species==spp,]
    #df_cs = df_s[df_s$chr==chr,]
    
    posnums = (table(df_cs$midPos))
    names = unique(names(posnums)[posnums>1])
    
    if(length(names)>1) {
      print(length(names))
      print(names)
      for(midPos in names) {
        print(paste(spp,chr,midPos))
        df = fix_column_dups(df_cs[df_cs$midPos==midPos,])
        dups[dups$species==spp & dups$chr==chr & dups$midPos==midPos,] = df
      }
      
      #dups = unique(dups)
      #write.table(dups,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore_2.FIXED",row.names = F,quote = F)
    } }
  #dups = unique(dups)
  #write.table(dups,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore_2.FIXED",row.names = F,quote = F)
}
dups = unique(dups)
not_dups = unique(not_dups)
#write.table(dups,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore_2.FIXED",row.names = F,quote = F)
#dim(dups)
merged = rbind(not_dups,dups)
merged = unique(merged)
posnums=rev(table(merged$sppchrpos))
names = unique(names(posnums)[posnums>1])
dups = merged[merged$sppchrpos %in% names,]

#d_a = dups[dups$mean_N_DATA==40,]
#write.table(d_a,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore_2.FIXED.bad.txt",row.names = F,quote = F)
#dups = dups[dups$mean_N_DATA!=40,]

not_dups = merged[!(merged$sppchrpos %in% names),]
write.table(not_dups,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore_2.NOTDUPS.txt",row.names = F,quote = F)
write.table(dups,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore_2.FIXED.txt",row.names = F,quote = F)
write.table(merged,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_taj_dxy_fst_islswp_miss_lostruct_zscore_2.temp",row.names = F,quote = F)

# plot(not_dups$rankDXY,not_dups$dxymeans,ylim=c(0.008,0.013),xlim=c(0.64,0.81))
# points(dups$rankDXY,dups$dxymeans,col="red")
# View(unique(not_dups[not_dups$rankDXY>0.7628 & not_dups$rankDXY<0.7629,c("dxymeans","rankDXY")]))



weird_missing = merged[merged$mean_N_MISS<100 & !(is.na(merged$mean_N_MISS)),]
for(spp in sort(unique(merged$species))) {
  temp=merged[merged$species==spp,]
  
  try({
    agg1 = aggregate(temp$weighted_recomb~temp$chr+temp$species,
                     FUN=function(x){mean(x,na.rm=T)},drop=F)
    print(paste(spp," recomb x 10^-10: ",signif(mean(agg1[,3],na.rm=T)/10^-10,3)," ± ",signif(sd(agg1[,3],na.rm=T)/10^-10,3)," (",sum(!(is.na(agg1[,3]))),")",sep=""))
    
  })
  
  try({
    agg2 = aggregate(temp$Fst~temp$chr+temp$species,
                     FUN=function(x){mean(x,na.rm=T)},drop=F)
    print(paste(spp," fst: ",signif(mean(agg2[,3],na.rm=T),3)," ± ",signif(sd(agg2[,3],na.rm=T),3)," (",sum(!(is.na(agg2[,3]))),")",sep=""))
    
  })
  
  try({
    agg3 = aggregate(temp$dxymeans~temp$chr+temp$species,
                     FUN=function(x){mean(x,na.rm=T)},drop=F)
    print(paste(spp," dxy: ",signif(mean(agg3[,3],na.rm=T),3)," ± ",signif(sd(agg3[,3],na.rm=T),3)," (",sum(!(is.na(agg3[,3]))),")",sep=""))
    
  })
  
  try({
    agg4 = aggregate(temp$Tajima~temp$chr+temp$species,
                     FUN=function(x){mean(x,na.rm=T)},drop=F)
    print(paste(spp," taj: ",signif(mean(agg4[,3],na.rm=T),3)," ± ",signif(sd(agg4[,3],na.rm=T),3)," (",sum(!(is.na(agg4[,3]))),")",sep=""))
    
    
  })
  
  try({
    agg5 = aggregate(temp$mean_F_MISS~temp$chr+temp$species,
                     FUN=function(x){mean(x,na.rm=T)},drop=F)
    
    print(paste(spp," missing: ",signif(mean(agg5[,3],na.rm=T),3)," ± ",signif(sd(agg5[,3],na.rm=T),3)," (",sum(!(is.na(agg5[,3]))),")",sep=""))
    
  })
  

}
}

merged$chr = as.character(merged$chr)
merged$species = as.character(merged$species)
merged$midPos = as.numeric(merged$midPos)

head(merged)


## need to add the three tajd -- remove the noweirds 

listfiles = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas/",
                       pattern="CHRFIX.pestPG",recursive=F,full.names = T)
listfiles = listfiles[!(grepl("NOWEIRD",listfiles))]

all_files = listfiles[!(grepl("CHI",listfiles)) & !(grepl("SON",listfiles))]
chi_files = listfiles[grepl("CHI",listfiles)]
son_files = listfiles[grepl("SON",listfiles)]

df_list=lapply(sort(unique(merged$species)),FUN=function(spp){
  print(spp)
  all = all_files[grepl(spp,all_files)]
  chi = chi_files[grepl(spp,chi_files)]
  son = son_files[grepl(spp,son_files)]
  df_all = read.table(all,header=T,strip.white = T,comment.char = "")
  df_chi = read.table(chi,header=T,strip.white = T,comment.char = "")
  df_son = read.table(son,header=T,strip.white = T,comment.char = "")
  df_all$species = spp
  df_chi$species = spp
  df_son$species = spp
  colnames(df_all)[colnames(df_all)=="WinCenter"] = "midPos"
  colnames(df_chi)[colnames(df_chi)=="WinCenter"] = "midPos"
  colnames(df_son)[colnames(df_son)=="WinCenter"] = "midPos"
  colnames(df_all)[colnames(df_all)=="Chr"] = "chr"
  colnames(df_chi)[colnames(df_chi)=="Chr"] = "chr"
  colnames(df_son)[colnames(df_son)=="Chr"] = "chr"
  df_all = df_all[,c("chr","midPos","species","Tajima")]
  df_son = df_son[,c("chr","midPos","species","Tajima")]
  df_chi = df_chi[,c("chr","midPos","species","Tajima")]
  colnames(df_all)[colnames(df_all)=="Tajima"] = "Tajima_all"
  colnames(df_chi)[colnames(df_chi)=="Tajima"] = "Tajima_chi"
  colnames(df_son)[colnames(df_son)=="Tajima"] = "Tajima_son"
  df = merge(df_all,merge(df_son,df_chi,all=T),all=T)
  
  df$chr = as.character(df$chr)
  df$species = as.character(df$species)
  df$midPos = as.numeric(df$midPos)
  
  return(df)
})

df = do.call(rbind,df_list)

merged2 = merge(merged,df,all=T)

png("test.png")
plot(merged2$Tajima,merged2$Tajima_all)
abline(a=0,b=1,col="red")
dev.off()

merged2$Tajima[is.na(merged2$Tajima) & !(is.na(merged2$Tajima_all))] = merged2$Tajima_all[is.na(merged2$Tajima) & !(is.na(merged2$Tajima_all))]

if(sum(merged2$Tajima!=merged2$Tajima_all,na.rm=T)==0){
  ## removed merged2 tajima_all
  merged2=merged2[ , ! names(merged2) %in% c("Tajima_all")]
}
write.table(merged2,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/rec_tajimafix_dxy_fst_islswp_miss_lostruct_zscore.temp",row.names = F,quote = F)

## now for recombination


listfiles = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/RECOMBINATION/",
                       pattern="PREDICT.txt_w100000_o10000_genome.txt$",recursive=T,full.names = T)
listfiles = listfiles[!(grepl("NOWEIRD",listfiles))]
listfiles = listfiles[!(grepl("Tgut",listfiles))]


all_files = listfiles[!(grepl("CHI",listfiles)) & !(grepl("SON",listfiles))]
chi_files = listfiles[grepl("CHI",listfiles)]
son_files = listfiles[grepl("SON",listfiles)]

## get the unique values of merged
merged_toadd = merged[,c("chr","midPos","windowstarts","windowstops")]
merged_toadd = unique(merged_toadd)

df_list=lapply(sort(unique(merged$species)),FUN=function(spp){
  print(spp)
  all = all_files[grepl(spp,all_files)]
  chi = chi_files[grepl(spp,chi_files)]
  son = son_files[grepl(spp,son_files)]
  df_all = read.table(all,header=T,strip.white = T,comment.char = "")
  df_chi = read.table(chi,header=T,strip.white = T,comment.char = "")
  df_son = read.table(son,header=T,strip.white = T,comment.char = "")
  df_all$species = spp
  df_chi$species = spp
  df_son$species = spp
  
  df_all$midPos = round((df_all$start+df_all$end) / 2)
  df_son$midPos = round((df_son$start+df_son$end) / 2)
  df_chi$midPos = round((df_chi$start+df_chi$end) / 2)
  
  colnames(df_all)[colnames(df_all)=="chrom"] = "chr"
  colnames(df_chi)[colnames(df_chi)=="chrom"] = "chr"
  colnames(df_son)[colnames(df_son)=="chrom"] = "chr"
  
  df_all$chr = sub("PseudoNC_","",df_all$chr)
  df_son$chr = sub("PseudoNC_","",df_son$chr)
  df_chi$chr = sub("PseudoNC_","",df_chi$chr)
  
  df_all$chr=as.character(sapply(df_all$chr,FUN=function(x){
    y=strsplit(x,"_")[[1]]
    return(y[length(y)])
    }))
  df_son$chr=as.character(sapply(df_son$chr,FUN=function(x){
    y=strsplit(x,"_")[[1]]
    return(y[length(y)])
  }))
  df_chi$chr=as.character(sapply(df_chi$chr,FUN=function(x){
    y=strsplit(x,"_")[[1]]
    return(y[length(y)])
  }))
  
  df_all = df_all[,c("chr","midPos","species","recombRate")]
  df_son = df_son[,c("chr","midPos","species","recombRate")]
  df_chi = df_chi[,c("chr","midPos","species","recombRate")]
  colnames(df_all)[colnames(df_all)=="recombRate"] = "recombRate_all"
  colnames(df_chi)[colnames(df_chi)=="recombRate"] = "recombRate_chi"
  colnames(df_son)[colnames(df_son)=="recombRate"] = "recombRate_son"
  df = merge(df_all,merge(df_son,df_chi,all=T),all=T)
  
  df$chr = as.character(df$chr)
  df$species = as.character(df$species)
  df$midPos = as.numeric(df$midPos)
  
  return(df)
})

