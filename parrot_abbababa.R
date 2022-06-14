with_melop = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/abbababa/Parrots_abbababa_jackknife.txt"
abbadf = read.csv(with_melop,sep="\t")
names(abbadf)

minagg = aggregate(abbadf$Dstat~abbadf$H1+abbadf$H2+abbadf$H3,FUN=function(x) {min(x,na.rm = T)})
maxagg = aggregate(abbadf$Dstat~abbadf$H1+abbadf$H2+abbadf$H3,FUN=function(x) {max(x,na.rm = T)})
meanagg = aggregate(abbadf$Dstat~abbadf$H1+abbadf$H2+abbadf$H3,FUN=function(x) {mean(x,na.rm = T)})
sdagg = aggregate(abbadf$Dstat~abbadf$H1+abbadf$H2+abbadf$H3,FUN=function(x) {sd(x,na.rm = T)})

colnames(minagg) = c("H1","H2","H3","MIN")
colnames(maxagg) = c("H1","H2","H3","MAX")
colnames(meanagg) = c("H1","H2","H3","MEAN")
colnames(sdagg) = c("H1","H2","H3","SD")

aggmerge = merge(minagg,maxagg)
aggmerge = merge(aggmerge,meanagg)
aggmerge = merge(aggmerge,sdagg)

aggmerge$H1 = substr(aggmerge$H1,1,3)
aggmerge$H2 = substr(aggmerge$H2,1,3)
aggmerge$H3 = substr(aggmerge$H3,1,3)

aggmerge$ALLH = paste(aggmerge$H1,aggmerge$H2,aggmerge$H3,sep="-")

aggmerge = aggmerge[order(aggmerge$MEAN),]

aggmerge=aggmerge[aggmerge$H1!="Cac",]
aggmerge=aggmerge[aggmerge$H1!="Ama",]
aggmerge=aggmerge[aggmerge$H1!="Ali",]
aggmerge=aggmerge[aggmerge$H2!="Cac",]
aggmerge=aggmerge[aggmerge$H2!="Ama",]
aggmerge=aggmerge[aggmerge$H2!="Ali",]
aggmerge=aggmerge[aggmerge$H3!="Cac",]
aggmerge=aggmerge[aggmerge$H3!="Ama",]
aggmerge=aggmerge[aggmerge$H3!="Ali",]

aggmerge=aggmerge[aggmerge$H1!="Apr",]
aggmerge=aggmerge[aggmerge$H2!="Apr",]


aggmerge=aggmerge[aggmerge$H1!=aggmerge$H2,]
aggmerge=aggmerge[aggmerge$H1!=aggmerge$H3,]
aggmerge=aggmerge[aggmerge$H2!=aggmerge$H3,]


indexes = 1:length(aggmerge$MEAN)

lowers = aggmerge$MEAN-aggmerge$SD
uppers = aggmerge$MEAN+aggmerge$SD

overlap1 = lowers < 0
overlap1.1 = lowers > 0
overlap2 = uppers > 0
overlap2.1 = uppers < 0 

overlapzero = overlap1 & overlap2
abovezero = overlap1.1 & overlap2
belowzero = overlap1 & overlap2.1

cols = as.numeric(overlapzero)+1
cols[overlapzero==T] = "black"
cols[abovezero==T] = "red"
cols[belowzero==T] = "cyan"

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/abbababa/Dstatistics_withmelop.png")
par(mar=c(4,7,0,0))
plot(aggmerge$MEAN,indexes,yaxt="n",xlab="Mean D-stat",ylab="",pch=16,
     col=cols,cex=1.5)
abline(v=0,col="black")
segments((aggmerge$MEAN-aggmerge$SD),indexes,aggmerge$MEAN+aggmerge$SD,indexes,col=cols)
#points(aggmerge$MEAN,indexes,pch=1,
#     col="black")
for(j in indexes) {
  axis(side=2, at=j, col.axis=as.character(cols)[j], labels=aggmerge$ALLH[j], las=2,cex=1.5) # Add habitat as labels, each with corresponding color, on the left margin
}
text(y=1,x=0.5,"H2<->H3",col="red")
text(y=1,x=-0.5,"H1<->H3",col="cyan")
dev.off()




## manual window abba baba with individuals
# D = [sum(ABBA) – sum(BABA)] / [sum(ABBA) + sum(BABA)]
abba="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/abbababa/Parrots_abbababa2_REHEAD_spp.abbababa"
babadf = read.csv(abba,sep="\t") ## big, takes a while
colnames(babadf) = substr(colnames(babadf),1,11)
unique(colnames(babadf))

headers=1:3
data_abba = seq(4,ncol(babadf)-1,2)
data_baba = seq(5,ncol(babadf),2)

abba_only = babadf[,c(headers,data_abba)]
baba_only = babadf[,c(headers,data_baba)]

## first, consolidate abba and baba
abba_trim = abba_only[,1:3]
colnames(abba_only) = substr(colnames(abba_only),1,11)
uniquecols = unique(colnames(abba_only))
uniquecols = uniquecols[c(-1,-2,-3)]
for(name in sort(uniquecols)){
 match = which(colnames(abba_only) %in% name)
 compiled = cbind(rowSums(abba_only[,match]))
 colnames(compiled) = name
 abba_trim=cbind(abba_trim,compiled)
}

baba_trim = baba_only[,1:3]
colnames(baba_only) = substr(colnames(baba_only),1,11)
uniquecols = unique(colnames(baba_only))
uniquecols = uniquecols[c(-1,-2,-3)]
for(name in sort(uniquecols)){
  match = which(colnames(baba_only) %in% name)
  compiled = cbind(rowSums(baba_only[,match]))
  colnames(compiled) = name
  baba_trim=cbind(baba_trim,compiled)
}

## exclude ones we don't need
#keep=c(1,2,3,14,15,16,17,18,23,24,25,26,27,38,39,40,41,42,43,44,45,46,47,58,59,60,61,62,68,69,70,71,72)
#denom=c(1,2,3,10,23,29,41,64,71)
#numer=c(15,24,26,39,43,44,46,59,61,69)
keepspp = c("ale.Apr.Apr","ale.ant.Apr","ale.ant.ant","ale.ant.swa","ale.swa.Apr","ale.swa.ant","ale.swa.swa","ant.Apr.Apr",
            "ant.ale.Apr","ant.ale.ale","ant.ale.swa","ant.swa.Apr","ant.swa.ale","ant.swa.swa","swa.Apr.Apr","swa.ale.Apr",
            "swa.ale.ale","swa.ale.ant","swa.ant.Apr","swa.ant.ale","swa.ant.ant")
denomspp =c("ale.Apr.Apr","ale.ant.ant","ale.swa.swa","ant.Apr.Apr","ant.ale.ale","ant.swa.swa","swa.Apr.Apr",
            "swa.ale.ale","swa.ant.ant")
numerspp =c("ale.ant.Apr","ale.ant.swa","ale.swa.Apr","ale.swa.ant","ant.ale.Apr","ant.ale.swa","ant.swa.Apr","ant.swa.ale","swa.ale.Apr",
            "swa.ale.ant","swa.ant.Apr","swa.ant.ale")

keep = which(colnames(abba_trim) %in% keepspp)
denom = which(colnames(abba_trim) %in% denomspp)
numer = which(colnames(abba_trim) %in% numerspp)

## F STATISTICS PASS 1 
abba_denom = abba_trim[,denom]
baba_denom = baba_trim[,denom]
abba_numer = abba_trim[,numer]
baba_numer = baba_trim[,numer]

abbababa_denom = abba_denom
abbababa_numer = abba_numer
for(colnum in 1:ncol(abba_denom)) {
  if(colnum %in% 1:3) {
    abbababa_denom[,colnum] = abba_denom[,colnum]
  } else {
    abba_colname = substr(colnames(abba_denom)[colnum],1,11)
    baba_colname = substr(colnames(baba_denom)[colnum],1,11)
    if(abba_colname==baba_colname) {
      abbababa_denom[,colnum] = (abba_denom[,colnum]-baba_denom[,colnum])
      colnames(abbababa_denom)[colnum] = abba_colname
    } else {
      print("ERRORRRRRRRR")
      print(colnum)
    }
  }
}
for(colnum in 1:ncol(abba_numer)) {
  if(colnum %in% 1:3) {
    abbababa_numer[,colnum] = abba_numer[,colnum]
  } else {
    abba_colname = substr(colnames(abba_numer)[colnum],1,11)
    baba_colname = substr(colnames(baba_numer)[colnum],1,11)
    if(abba_colname==baba_colname) {
      abbababa_numer[,colnum] = (abba_numer[,colnum]-baba_numer[,colnum])
      colnames(abbababa_numer)[colnum] = abba_colname
    } else {
      #print("ERRORRRRRRRR")
      print(paste(colnum,abba_colname,baba_colname))
    }
  }
}

## compare the first and third

colnames(abbababa_denom) = substr(colnames(abbababa_denom),1,11)
colnames(abbababa_numer) = substr(colnames(abbababa_numer),1,11)

spp_denom = colnames(abbababa_denom)
spp_numer = colnames(abbababa_numer)

spp_denom = cbind(substr((spp_denom),1,3),
                  substr((spp_denom),5,7),
                  substr((spp_denom),9,11))

spp_denom_13 = paste(spp_denom[,1],spp_denom[,3])
spp_denom_23 = paste(spp_denom[,2],spp_denom[,3])

## find the ones where 1 and 3 are the same 

abbababa_fstat = abbababa_numer
for (colnum in 4:ncol(abbababa_numer)) {
  print(colnum)
  data = abbababa_numer[,colnum]
  column_name = colnames(abbababa_numer)[colnum]
  temp = cbind(substr((column_name),1,3),
        substr((column_name),5,7),
        substr((column_name),9,11))
  
  temp_13 = paste(temp[,1],temp[,3])
  temp_23 = paste(temp[,2],temp[,3])
  
  match13 = which(spp_denom_13 %in% temp_13)
  #match23 = which(spp_denom_23 %in% temp_23)
  print(match13)
  #print(match23)
  print("xxxxxx")
  average_denom_13 = (abbababa_denom[,match13])
  #average_denom_23 = (abbababa_denom[,match23])
  #average_denom = rowSums(cbind(average_denom_13,average_denom_23))
  
  abbababa_fstat[,colnum] = data/average_denom_13
  
}

colnames(abbababa_fstat) = colnames(abbababa_numer)
abbababa_fstat = cbind(babadf[,1:3],abbababa_fstat)
colnames(abbababa_fstat)[4:ncol(abbababa_fstat)] = substr(colnames(abbababa_fstat)[4:ncol(abbababa_fstat)],1,11)
  
write.csv(abbababa_fstat,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/abbababa/F-statistics_calculated.abbababa2.txt",row.names = F)

fstat = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/abbababa/F-statistics_calculated.abbababa2.txt")
plot(1:nrow(fstat),as.numeric(as.character(fstat[,4])),
     col=as.numeric(as.factor(fstat[,1])))


chromtype = substr(as.character(fstat[,1]),1,8)
smallagg = fstat[chromtype %in% c("PseudoNC"),]
#smallagg = smallagg[,-13]

smallagg_bounded = smallagg[,4:ncol(smallagg)]
smallagg_bounded[smallagg_bounded >= 0.5] = 0.6
smallagg_bounded[smallagg_bounded <= -0.5] = -0.6
smallagg_bounded = cbind(smallagg[,1:3],smallagg_bounded)

colors = as.numeric(as.factor(as.character(smallagg[,1]))) %% 10
colors = c("red","orange","goldenrod","green","blue","purple",
           "cyan","black","brown","magenta")[colors]

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/abbababa/F-statstics_aggregated_new.png",
    height=450,width=700)
par(mfrow=c(2,3))
for(i in 4:ncol(smallagg_bounded)) {
  print(i)
  plot(1:nrow(smallagg_bounded),(as.numeric(as.character(smallagg_bounded[,i]))),col=colors,type="p",pch=16,main=colnames(smallagg)[i],
       xlab="Scaffold",ylab="F-statistic",ylim=c(-0.7,0.7))
  abline(h=0.51)
  abline(h=-0.51)
}
dev.off()



#ABBA_1_2_3 = abba(freq_table[,P1], freq_table[,P2], freq_table[,P3])
#BABA_1_2_3 = baba(freq_table[,P1], freq_table[,P2], freq_table[,P3])
#ABBA_1_3_3 = abba(freq_table[,P1], freq_table[,P3], freq_table[,P3])
#BABA_1_3_3 = baba(freq_table[,P1], freq_table[,P3], freq_table[,P3])
#f = (sum(ABBA_1_2_3) - sum(BABA_1_2_3))/  (sum(ABBA_1_3_3) - sum(BABA_1_3_3))




D_df = babadf[,1:3]

lengthtodo = ncol(babadf)

for(i in data_abba) {
  print(paste(i,i+1,"/",lengthtodo))
  comparison = colnames(babadf)[i]
  print("xxxxxx")
  #D = [sum(ABBA) – sum(BABA)] / [sum(ABBA) + sum(BABA)]
  Dstat = (babadf[,i]-babadf[,i+1]) / ((babadf[,i]+babadf[,i+1]))
  Dstat = as.data.frame(Dstat)
  colnames(Dstat) = comparison
  D_df = cbind(D_df,Dstat)
}

write.csv(D_df,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/abbababa/D-statistics_calculated.abbababa2.txt",row.names = F)

dstatrehead = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/abbababa/D-statistics_calculated.abbababa2_REHEAD.txt",header = F)

transposed = t(dstatrehead)

headers = transposed[1:3,]

test = transposed[4:nrow(transposed),]
test = as.data.frame(test)
test = taRifx::remove.factors(test)

for (colnum in 2:ncol(test)) {
  if(colnum %% 100 == 0) { 
  print(paste(colnum,"/",ncol(test)))
  }
  test[,colnum] = as.numeric(test[,colnum])
}




aggdf = c()
for(species in sort(unique(test$V1))) {
  species = as.character(species)
  print(species)
  subtest = test[test$V1==species,]
  print(subtest[1:3,1:min(5,nrow(subtest))])
  agg = colMeans(subtest[,2:ncol(subtest)],na.rm = T)
  agg = c(species,agg)
  
  if(is.null(aggdf)) {
    aggdf = agg
  } else {
    aggdf =  rbind(aggdf, agg)
  }
}

combined_agg = rbind(headers,aggdf)
combined_agg = t(combined_agg)
colnames(combined_agg) = combined_agg[1,]
combined_agg = combined_agg[-1,]

combined_agg = combined_agg[,c(1:3,15,24,26,27,39,43,44,46,59:61,69)]

write.csv(combined_agg,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/abbababa/D-statistics_aggregated_trimmed.abbababa2.txt",row.names = F)

chromtype = substr(as.character(combined_agg[,1]),1,8)
smallagg = combined_agg[chromtype %in% c("PseudoNC"),]

colors = as.numeric(as.factor(as.character(smallagg[,1]))) %% 10
colors = c("red","orange","goldenrod","green","blue","purple",
           "cyan","black","brown","magenta")[colors]


png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/abbababa/D-statstics_aggregated.png",
    height=600,width=800)
par(mfrow=c(3,4))
for(i in 4:ncol(smallagg)) {
  plot(as.numeric(smallagg[,i]),col=colors,type="p",pch=16,ylim=c(-1,1),main=colnames(smallagg)[i],
       xlab="Scaffold",ylab="D-statistic")
abline(h=0,col="black")
  }
dev.off()


write.csv(table(smallagg[,1],smallagg[,3]),"test.txt")

