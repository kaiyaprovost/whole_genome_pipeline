#filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/merged_empirical_stats_same_cols_TRIMMED.txt"
#filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/all_merged_8_april_2020_TRIMMED_max_empirical.txt"
filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/MERGED_COMPLETE_SUMSTATS_MLPIPELINE_29may2020.txt" ## old: - final .csv
df = read.table(filename,sep="\t",stringsAsFactors = F,header=T,fill=T)
file2="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/STATS/MERGED_all_merged_29_may_2020.txt" ## old: 8 april
df2 = read.table(file2,sep="\t",stringsAsFactors = F,header=T,fill=T)

tail(df[df$DEMOG=="EMPIRICAL",goodcols])
df$YEAR[df$DEMOG=="EMPIRICAL"] = "EMPIRICAL"
df$MIGRATION.RATE[df$DEMOG=="EMPIRICAL"] = "EMPIRICAL"
df$MUTATIONPOWER[df$DEMOG=="EMPIRICAL"] = "EMPIRICAL"
df$NUMRUN[df$DEMOG=="EMPIRICAL"] = "EMPIRICAL"
df$SCALE[df$DEMOG=="EMPIRICAL"] = "EMPIRICAL"
df$SECCON[df$DEMOG=="EMPIRICAL"] = "EMPIRICAL"
df$IBD[df$DEMOG=="EMPIRICAL"] = "EMPIRICAL"
df$EFFPOPSIZE[df$DEMOG=="EMPIRICAL"] = "UNK"
df$RUN[df$DEMOG=="EMPIRICAL"] = "EMPIRICAL"
df$WINDOWNUM[is.na(df$WINDOWNUM)] = df$WINDOW[is.na(df$WINDOWNUM)]
df$WINDOW[is.na(df$WINDOW)] = df$WINDOWNUM[is.na(df$WINDOW)]
df$WINDOWSIZE[is.na(df$WINDOWSIZE)] = df$WINDOW.SIZE[is.na(df$WINDOWSIZE)]
df$WINDOW.SIZE[is.na(df$WINDOW.SIZE)] = df$WINDOWSIZE[is.na(df$WINDOW.SIZE)]
df$SPECIES[df$FILE==""] = "SIMULATION"
df$SPECIES[df$SPECIES==""] = "SIMULATION"
df$SPECIES[grepl("RECAP",df$FILE)] = "SIMULATION"
df=unique(df)
write.table(df,"/Users/kprovost/Dropbox (AMNH)/Dissertation/MERGED_COMPLETE_SUMSTATS_MLPIPELINE_29may2020.txt",
            sep="\t",row.names = F)

df2 = df2[df2$DEMOG!="EMPIRICAL",]
df2$WINDOW = 0
df2$WINDOW.SIZE = 100000
df2$WINDOWNUM = 0
df2$WINDOWSIZE = 100000
df2$SPECIES = "SIMULATION"
df2$SPECIES[grepl("TGUT",df2$FILE)] = "EMPIRICAL"
df2$DEMOG[grepl("TGUT",df2$FILE)] = "EMPIRICAL"
df2$SPECIES[grepl("-CALLED",df2$FILE)] = "EMPIRICAL"
df2$DEMOG[grepl("-CALLED",df2$FILE)] = "EMPIRICAL"
df2$FILE[grepl("-CALLED",df2$FILE)] = "EMPIRICAL"
df2 = df2[df2$DEMOG!="EMPIRICAL",]
df2 = df2[df2$SPECIES!="EMPIRICAL",]
df2 = df2[df2$YEAR!="EMPIRICAL",]
df2 = df2[df2$FILE!="EMPIRICAL",]
df2$FILE[df2$FILE==""] = "SIMULATION"
df2$CHROMOSOME[df2$CHROMOSOME==""] = "SIMULATION"
df2$EFFPOPSIZE[df2$EFFPOPSIZE==""] = 400000
df2$RECOMPOWER[df2$RECOMPOWER==""] = -8
df2$MUTATIONPOWER[df2$MUTATIONPOWER==""] = -9
df2$SECCON[df2$DEMOG=="SECCON"] = 1
df2$SECCON[df2$DEMOG!="SECCON"] = 0
df2=df2[df2$DEMOG!="SCALE1",]
df2$CHROM=df2$CHROMOSOME
df2$RUN[df2$RUN==""] = "SIMULATION"
df2$OVERLAP=100000
df2$OVERLAP.SIZE=100000
df2$NUMRUN[is.na(df2$NUMRUN)] = paste(df2$DEMOG[is.na(df2$NUMRUN)],df2$IBD[is.na(df2$NUMRUN)])
df2$NUMRUN[df2$NUMRUN=="PANMIXIA FALSE"]=1
df2$NUMRUN[df2$NUMRUN=="PANMIXIA TRUE"]=2
df2$NUMRUN[df2$NUMRUN=="ISOLATION FALSE"]=3
df2$NUMRUN[df2$NUMRUN=="ISOLATION TRUE"]=4
df2$NUMRUN[df2$NUMRUN=="GENEFLOW FALSE"]=5
df2$NUMRUN[df2$NUMRUN=="GENEFLOW TRUE"]=6
df2$NUMRUN[df2$NUMRUN=="SECCON FALSE"]=7
df2$NUMRUN[df2$NUMRUN=="SECCON TRUE"]=8
df2$MIGRATION.RATE[df2$NUMRUN>=5] = 0.1
df2$MIGRATION.RATE[df2$NUMRUN<5] = 0
df2 = unique(df2)
write.table(df2,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/STATS/MERGED_all_merged_29_may_2020.txt",
            row.names = F,sep="\t")

goodcols = sort(unique(c("CHROM","DEMOG","SPECIES","FILE","EFFPOPSIZE","IBD","MIGRATION.RATE",
             "MUTATIONPOWER","NUMRUN","OVERLAP.SIZE","SCALE",
             "RUN","WINDOW","WINDOW.SIZE","YEAR","OVERLAP","SECCON",
             "WINDOWSIZE","WINDOWNUM","CHROMOSOME")))

agg = merge(df2,df,all=T)
#agg=df

agg=unique(agg)

badcols = c()
for (colnumber in 1:ncol(agg)){
  thiscol = unique(agg[,colnumber])
  thiscol[thiscol==""] = NA
  thiscol[thiscol=="NA"] = NA
  thiscol = thiscol[!(is.na(thiscol))]
  thiscol = thiscol[!(is.null(thiscol))]
  thiscol=unique(thiscol)
  print(thiscol[1:min(5,length(thiscol))])
  if(length(thiscol)<=1){
    badcols=c(badcols,colnames(agg)[colnumber])
  }
}
badcols = badcols[!(badcols %in% goodcols)]
agg = agg[,-c(which(colnames(agg) %in% badcols))]

badcols = c()
for(colnumber in 1:ncol(agg)) {
  nonblank=sum(!(is.na(agg[,colnumber])) & (agg[,colnumber] != ""))
  cat(nonblank," ")
  if(nonblank <= 1000){
    badcols=c(badcols,colnames(agg)[colnumber])
  }
}
badcols = badcols[!(badcols %in% goodcols)]

agg = agg[,-c(which(colnames(agg) %in% badcols))]
agg = unique(agg)

agg_emp = agg[agg$DEMOG=="EMPIRICAL",]
agg_noemp = agg[agg$DEMOG!="EMPIRICAL",]

badcols = c()
for(colnumber in 1:ncol(agg_emp)) {
  nonblank=sum(!(is.na(agg_emp[,colnumber])) & (agg_emp[,colnumber] != ""))
  print(nonblank)
  if(nonblank <= 1000){
    badcols=c(badcols,colnames(agg_emp)[colnumber])
  }
}
badcols = badcols[!(badcols %in% goodcols)]
agg_emp = agg_emp[,-c(which(colnames(agg_emp) %in% badcols))]

badcols = c()
for(colnumber in 1:ncol(agg_noemp)) {
  nonblank=sum(!(is.na(agg_noemp[,colnumber])) & (agg_noemp[,colnumber] != ""))
  cat(nonblank, " ")
  if(nonblank <= 10000){
    badcols=c(badcols,colnames(agg_noemp)[colnumber])
  }
}
badcols = badcols[!(badcols %in% goodcols)]
if(length(badcols)>0){
agg_noemp = agg_noemp[,-c(which(colnames(agg_noemp) %in% badcols))]
}

colskeep = intersect(names(agg_emp),names(agg_noemp))

agg = agg[,colskeep]
agg = unique(agg)
head(agg)

## calculate pca 8:12,15:80

agg_meta = agg[,goodcols]
agg_numb = agg[,-which(colnames(agg) %in% goodcols)]

agg_meta=agg_meta[complete.cases(agg_numb),]
agg_numb=agg_numb[complete.cases(agg_numb),]

pca = prcomp(agg_numb,center=T,scale. = T)

summary(pca)
#pcadata = cbind(good_notnum,pca$x)
pcadata = cbind(agg_meta,pca$x)

palette(c(RColorBrewer::brewer.pal(8,"Dark2"),
          RColorBrewer::brewer.pal(12,"Paired"))
        )
#plot(pcadata$PC1,pcadata$PC2,col=as.numeric(as.factor(pcadata$DEMOG))-1,
#     pch=as.numeric(as.factor(pcadata$SPECIES)))
plot(pcadata$PC1,pcadata$PC2,col=as.numeric(as.factor(pcadata$SPECIES)),
     pch=as.numeric(as.factor(pcadata$SPECIES)))

pcadata$ISSIM = as.numeric(!(pcadata$DEMOG %in% c("EMPIRICAL")))
plot(pcadata$PC1,pcadata$PC2,col="black",cex=0.75,pch=22,
     #bg=as.numeric(as.factor(pcadata$SPECIES)),
     bg=pcadata$ISSIM,
     lwd=0.5#,xlim=c(-5,5),ylim=c(-5,5)
     )


pca_emp = pcadata[pcadata$ISSIM==0,] 
## pc1 -11.3183 to 5.3304
## pc2 -10.4567 to 7.8451
pca_sim = pcadata[pcadata$ISSIM==1,]
## pc1 -4.989 to 11.180
## pc2 -1.207 to 6.469

## overlap zone PC1: -4.989 - 5.3304
## overlap zone PC2: -1.207 - 6.469

subset = pcadata[pcadata$PC1 >= min(pca_sim$PC1),]
subset = subset[subset$PC1 <= max(pca_emp$PC1),]
subset = subset[subset$PC2 <= max(pca_sim$PC2),]
subset = subset[subset$PC2 >= min(pca_sim$PC2),]
subset = subset[subset$PC2 <= 4,]
subset = subset[(subset$PC1+subset$PC2) <= 4,]
subset = subset[(subset$PC1+subset$PC2) >= -4,]
subset = subset[(subset$PC1*subset$PC2) >= -10,]

plot(subset$PC1,subset$PC2,col="black",cex=0.75,pch=22,
     #bg=as.numeric(as.factor(pcadata$SPECIES)),
     bg=subset$ISSIM,
     lwd=0.5#,xlim=c(-5,5),ylim=c(-5,5)
)
points(pcadata$PC1[pcadata$SPECIES=="SIMULATION"],
       pcadata$PC2[pcadata$SPECIES=="SIMULATION"],
       cex=1,pch=16,
       col=as.numeric(as.factor(pcadata$DEMOG[pcadata$SPECIES=="SIMULATION"])))
legend(x="bottomleft",cex=1,
       legend=unique(as.factor(
         pcadata$DEMOG[pcadata$SPECIES=="SIMULATION"])),
       col=unique(as.numeric(as.factor(
         pcadata$DEMOG[pcadata$SPECIES=="SIMULATION"]))),
       pch=16)


subset[subset$SPECIES=="SIMULATION",]



boxplot(pcadata$PC1~pcadata$ISSIM)
boxplot(pcadata$PC2~pcadata$ISSIM)




for(i in 1:length(levels(as.factor(pcadata$SPECIES)))) {
  spp=unique(levels(as.factor(pcadata$SPECIES)))[i]
  
car::ellipse(center = colMeans(pcadata[pcadata$SPECIES==spp,
                                       c("PC1","PC2")]), 
        shape = cov(pcadata[pcadata$SPECIES==spp,
                            c("PC2","PC1")]),
        center.pch=1,center.cex=2,lwd=3,lty=1,
        radius = sqrt(qchisq(.5, df=2)),col = i)
}
legend("bottomright",
       legend=levels(as.factor(pcadata$SPECIES)),
       col=1:length(levels(as.factor(pcadata$SPECIES))),
       fill=1:length(levels(as.factor(pcadata$SPECIES))),
       pch=22,
       cex=0.5,
       ncol=2)

cl=kmeans(pca$x,centers=11)
plot(pcadata$PC1,pcadata$PC2,col="lightgrey",cex=0.3,pch=1)
points(pcadata$PC1[pcadata$SPECIES=="SIMULATION"],
       pcadata$PC2[pcadata$SPECIES=="SIMULATION"],cex=2,col="red")
points(pcadata$PC1,pcadata$PC2,col=cl$cluster,pch=cl$cluster)
#points(cl$centers,col=1:11,pch=1:11,cex=2)


#df$GENOME.haplotype.diversity = as.numeric(as.character(df$GENOME.haplotype.diversity))


df = df[complete.cases(df),]

#numeric = df[,c(8:12,15:80)]
#notnumeric = df[,c(1:7,13:14)] ## 1:7, 13:14

notnumeric=df[,c(1,2,4,5,10,16,17,23,25,26,30,42,43)]
numeric=df[,-c(1,2,4,5,10,16,17,23,25,26,30,42,43)]


#sims=df[df$species=="simulation",]
#emps=df[df$species!="simulation",]
## sims missing: 76:80, 38:48
## emps missing: 49:80

sims=df[df$SPECIES=="SIMULATION",]
emps=df[df$SPECIES!="SIMULATION",]

#good_for_both = df[,c(1:37)]
#good_for_both = good_for_both[complete.cases(good_for_both),]
#good_num = good_for_both[,c(8:12,15:37)]
#good_notnum = good_for_both[,c(1:7,13:14)]

summary(numeric)
summary(notnumeric)

#pca = prcomp(good_num,center=T,scale. = T)
pca = prcomp(numeric,center=T,scale. = T)

summary(pca)
#pcadata = cbind(good_notnum,pca$x)
pcadata = cbind(notnumeric,pca$x)

palette(viridis::viridis(11))
#plot(pcadata$PC1,pcadata$PC2,col=as.numeric(as.factor(pcadata$DEMOG))-1,
#     pch=as.numeric(as.factor(pcadata$SPECIES)))
plot(pcadata$PC1,pcadata$PC2,col=as.numeric(as.factor(pcadata$SPECIES)),
     pch=as.numeric(as.factor(pcadata$SPECIES)))

palette(viridis::viridis(3))
plot(pcadata$PC1,pcadata$PC2,col=as.numeric(as.factor(pcadata$SIMDEFAULT.)),
     pch=as.numeric(as.factor(pcadata$SPECIES)))

par(mfrow=c(3,4))
for(spp in unique(pcadata$species)){
  toplot = pcadata[pcadata$species==spp,]
  plot(pcadata$PC1,pcadata$PC2,main=spp,col="lightgrey")
  points(toplot$PC1,toplot$PC2,col="red")
}

palette(c("white","red","orange","goldenrod","green",
          "cyan","blue","magenta","white","purple",
          "brown","black","grey","darkred"))

emponly = pcadata[pcadata$species!="simulation",]
emponly = emponly[emponly$species!="",]
plot(emponly$PC1,emponly$PC2,col=as.numeric(as.factor(emponly$species)),
     pch=as.numeric(as.factor(emponly$species)))
legend("topleft",legend=unique(as.factor(emponly$species)),col=unique(as.numeric(as.factor(emponly$species))),
       pch=unique(as.numeric(as.factor(emponly$species))),
       bty="n",cex=0.75)

par(mfrow=c(3,4))
for(spp in unique(emponly$species)){
  toplot = emponly[emponly$species==spp,]
  plot(emponly$PC1,emponly$PC2,main=spp,col="lightgrey")
  points(toplot$PC1,toplot$PC2,col="red")
}


good_for_emp = emps[,c(1:37)]
good_for_emp = good_for_emp[complete.cases(good_for_emp),]
good_emp = good_for_emp[,c(8:12,15:37)]
good_empnot = good_for_emp[,c(1:7,13:14)]
pcaemp = prcomp(good_emp,center=T,scale. = T)
summary(pcaemp)
pcadataemp = cbind(good_empnot,pcaemp$x)
plot(pcadataemp$PC1,pcadataemp$PC2,col=as.numeric(as.factor(pcadataemp$species)),
     pch=as.numeric(as.factor(pcadataemp$species)),
     ylim=c(-7,10),xlim=c(-10,10))
legend("topleft",legend=unique(as.factor(pcadataemp$species)),
       col=unique(as.numeric(as.factor(pcadataemp$species))),
       pch=unique(as.numeric(as.factor(pcadataemp$species))),
       bty="n",cex=0.75)


## analysis where for each window you look and find 10 most nearby sims












#####

trimmed = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/merged_empirical_stats_same_cols_TRIMMED_extra.txt"
scaled = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/MERGEDPOP_TOGETHER.txt"
bigdata = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS_testing/TESTING/RECAPSTATS/MERGEDPOP_FINAL_trimmed.txt"

t_df = read.csv(trimmed,sep="\t")
s_df = read.csv(scaled,sep="\t")
b_df = read.csv(bigdata,sep="\t")

all_df = smartbind(t_df,s_df,b_df)
write.csv(all_df,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/alldf.txt")

t_s_overlap = intersect(colnames(t_df),colnames(s_df))
t_b_overlap = intersect(colnames(t_df),colnames(b_df))
s_b_overlap = intersect(colnames(b_df),colnames(s_df))
t_s_b_overlap = intersect(t_s_overlap,colnames(b_df))
t_s_b_overlap = c(t_s_b_overlap,colnames(t_df)[c(1:7,13:14)])


t_trim = t_df[,which(colnames(t_df) %in% t_s_b_overlap)]
s_trim = s_df[,which(colnames(s_df) %in% t_s_b_overlap)]
b_trim = b_df[,which(colnames(b_df) %in% t_s_b_overlap)]


all_trim = smartbind(t_trim,s_trim,b_trim)
all_trim = unique(all_trim)
all_trim = all_trim[complete.cases(all_trim),]

#####

all_trim = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/alldf_trim.txt",
                     sep="\t",row.names=NULL)
all_trim = all_trim[,c(1:13,16:17,19:33,37:40)]
all_trim$WINDOW[is.na(all_trim$WINDOW)] = 0
all_trim$WINDOW.SIZE[is.na(all_trim$WINDOW.SIZE)] = 100000
all_trim$OVERLAP.SIZE[is.na(all_trim$OVERLAP.SIZE)] = 100000
all_trim = all_trim[complete.cases(all_trim),]

all_trim_num = all_trim[,c(9:ncol(all_trim))]
all_trim_notnum = all_trim[,c(1:8)]

pca = prcomp(all_trim_num,center=T,scale. = T)
summary(pca)
pcadata = cbind(all_trim_notnum,pca$x)
plot(pcadata$PC1,pcadata$PC2,col=as.numeric(as.factor(pcadata$DEMOG)),
     pch=as.numeric(as.factor(pcadata$SPECIES)))

unique(pcadata$SPECIES)

pcadata$SPECIES = as.character(pcadata$SPECIES)
pcadata$SPECIES[pcadata$SPECIES==""] = as.character(pcadata$DEMOG[pcadata$SPECIES==""])

plot(pcadata$PC1,pcadata$PC2,col=as.numeric(as.factor(pcadata$DEMOG)),
     pch=as.numeric(as.factor(pcadata$SPECIES)),type="n")

par(mfrow=c(2,2))
for (i in 2:length(unique(pcadata$DEMOG))) {
  category = unique(pcadata$DEMOG)[i]
  d = pcadata[pcadata$DEMOG==category,]
  plot(d$PC1,d$PC2,main=category)
  
  
  for (i in 1:10) {
    category = unique(pcadata$SPECIES)[i]
    d = pcadata[pcadata$SPECIES==category,]
    
    car::ellipse(center = colMeans( pcadata[pcadata$SPECIES==category,c("PC2","PC1")],na.rm = T), 
                 shape = cov( pcadata[pcadata$SPECIES==category,c("PC2","PC1")],use="pairwise.complete.obs"),
                 center.pch=i,center.cex=2,lwd=2,
                 radius = sqrt(qchisq(.95, df=2)),col = i,
                 lty=1)
  }
  
}

for (i in 1:length(unique(pcadata$DEMOG))) {
  category = unique(pcadata$DEMOG)[i]
  d = pcadata[pcadata$DEMOG==category,]
  
car::ellipse(center = colMeans( pcadata[pcadata$DEMOG==category,c("PC1","PC2")],na.rm = T), 
             shape = cov( pcadata[pcadata$DEMOG==category,c("PC1","PC2")],use="pairwise.complete.obs"),
             center.pch=i,center.cex=2,lwd=2,
             radius = sqrt(qchisq(.5, df=2)),col = i,
             lty=1)
}

plot(pcadata$PC1,pcadata$PC2,col=as.numeric(as.factor(pcadata$DEMOG)),
     pch=as.numeric(as.factor(pcadata$SPECIES)),type="n")
for (i in 1:length(unique(pcadata$YEAR))) {
  category = unique(pcadata$YEAR)[i]
  d = pcadata[pcadata$YEAR==category,]
  
  car::ellipse(center = colMeans( pcadata[pcadata$YEAR==category,c("PC1","PC2")],na.rm = T), 
               shape = cov( pcadata[pcadata$YEAR==category,c("PC1","PC2")],use="pairwise.complete.obs"),
               center.pch=i,center.cex=2,lwd=2,
               radius = sqrt(qchisq(.5, df=2)),col = i,
               lty=1)
}

plot(pcadata$PC2,pcadata$PC1,col=as.numeric(as.factor(pcadata$DEMOG)),
     pch=as.numeric(as.factor(pcadata$SPECIES)),type="n")
for (i in 1:length(unique(pcadata$SPECIES))) {
  category = unique(pcadata$SPECIES)[i]
  d = pcadata[pcadata$SPECIES==category,]
  
  car::ellipse(center = colMeans( pcadata[pcadata$SPECIES==category,c("PC2","PC1")],na.rm = T), 
               shape = cov( pcadata[pcadata$SPECIES==category,c("PC2","PC1")],use="pairwise.complete.obs"),
               center.pch=i,center.cex=2,lwd=2,
               radius = sqrt(qchisq(.5, df=2)),col = i,
               lty=1)
}




