filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/merged_empirical_stats_same_cols_TRIMMED.txt"
df = read.csv(filename,sep="\t")
## calculate pca 8:12,15:80

df$GENOME.haplotype.diversity = as.numeric(as.character(df$GENOME.haplotype.diversity))


numeric = df[,c(8:12,15:80)]
notnumeric = df[,c(1:7,13:14)] ## 1:7, 13:14

sims=df[df$species=="simulation",]
emps=df[df$species!="simulation",]
## sims missing: 76:80, 38:48
## emps missing: 49:80

good_for_both = df[,c(1:37)]
good_for_both = good_for_both[complete.cases(good_for_both),]
good_num = good_for_both[,c(8:12,15:37)]
good_notnum = good_for_both[,c(1:7,13:14)]

#summary(numeric)
#summary(notnumeric)

pca = prcomp(good_num,center=T,scale. = T)
summary(pca)
pcadata = cbind(good_notnum,pca$x)
plot(pcadata$PC1,pcadata$PC2,col=as.numeric(as.factor(pcadata$demog))-1,
     pch=as.numeric(as.factor(pcadata$species)))

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
