## summarize fst across species across chromosomes

setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/")

spplist=c("bellii",
          "bilineata",
          "brunneicapillus",
          "crissale",
          "curvirostre","flaviceps",
          "fusca","melanura",
          "nitens","sinuatus"
          )

bigord = c()

smallplot=F
doFST = F
doDXY = F
doBIG = F

if(doFST == T){
png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/FST/3.slidingWindowJobs/fst_boxplot_with_best_model.png",width=800)
par(mfrow=c(2,5),mar=c(4,4,1,0))
for (species in spplist) {

#species="crissale"
capspecies=toupper(species)


colors=RColorBrewer::brewer.pal(8,"Dark2")
palette(c("white",colors[1],colors[3],colors[7],colors[4],"black",
          colors[c(2,5,6,8)]))

path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/FST/3.slidingWindowJobs/"

setwd(path)

filename=list.files(path=path,
                    recursive = F,
                    pattern=species)

if(length(filename) != 0) {

df = read.table(filename,
                header=T)
df=df[order(df$chr, df$midPos),]

meanagg = aggregate(df$Fst~df$chr,FUN=function(x){mean(x,na.rm=T)})
colnames(meanagg) = c("chrom","meanFst")
sdagg = aggregate(df$Fst~df$chr,FUN=function(x){sd(x,na.rm=T)})
colnames(sdagg) = c("chrom","sdFst")

countsagg=as.data.frame(cbind(table(df$chr)))
colnames(countsagg) = "countFst"
countsagg$chrom = rownames(countsagg)

agg = merge(meanagg,sdagg)
agg = merge(agg,countsagg)
agg$sdFst[is.na(agg$sdFst)]=0
agg$stderrFst = agg$sdFst / sqrt(agg$countFst)

models=read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/bivariate/bivariate_gdm_results.csv",
                skip=1)
thismod=models[models$SPECIES==capspecies,]
thismod=thismod[,c("NUMERICDATASET","MOSTA.1MODEL1",
                   "MOSTA.1MODEL2","MOSTA.1MODEL3")]
colnames(thismod) = c("chrom","m1","m2","m3")
thismod$m2[is.na(thismod$m2)] = ""

agg=merge(agg,thismod)

if(smallplot==T){
par(mar=c(4,4,0,0))
boxplot(df$Fst~df$chr,las=2,col="grey",
        border=as.numeric(agg$m2),ylab="Fst",xlab="")
legend("topleft",legend=unique((agg$m2)),
       fill=unique(as.numeric(agg$m2)),bty="n",
       ncol=3,title=species)
readline("press enter: ")
}

agg=agg[order(agg$m2, agg$meanFst, agg$chrom),]

if(smallplot==T){
#par(mfrow=c(3,1))
par(mar=c(4,4,0.1,0))
# b=barplot(agg$meanFst,names=agg$chrom,las=2,
#           col=as.numeric(agg$m1))
# segments(b, agg$meanFst-agg$stderrFst,
#          b, agg$meanFst+agg$stderrFst, lwd=2)
# legend("top",legend=unique((agg$m1)),
#        fill=unique(as.numeric(agg$m1)),bty="n",ncol=6)
b=barplot(agg$meanFst,names=agg$chrom,las=2,cex.names=0.75,
          col=as.numeric(agg$m2))
segments(b, agg$meanFst-agg$stderrFst,
         b, agg$meanFst+agg$stderrFst, lwd=2)
legend("top",legend=unique((agg$m2)),
       fill=unique(as.numeric(agg$m2)),bty="n",ncol=3,title=species)
# b=barplot(agg$meanFst,names=agg$chrom,las=2,
#           col=as.numeric(agg$m3))
# segments(b, agg$meanFst-agg$stderrFst,
#          b, agg$meanFst+agg$stderrFst, lwd=2)
# legend("top",legend=unique((agg$m3)),
#        fill=unique(as.numeric(agg$m3)),bty="n",ncol=6)
readline("press enter: ")
}

agg2 = aggregate(agg$meanFst~agg$m2,FUN=function(x)
  {median(x,na.rm = T)}
  )


ord=order(-agg2[,2])

p=boxplot(agg$meanFst~agg$m2,col=1:6,
          xlab=species,ylab="Mean Fst",las=2)

ord2=order(-p$stats[3,])

print(species)
if(sum(ord2==ord)==length(ord)) {
  print("FINE")
  names(ord) = agg2[,1]
  
} else {
  print("MISMATCH")
  print(ord)
  print(ord2)
  ord=ord2
  names(ord) = p$names
  
}

ord = as.data.frame(ord)
colnames(ord) = species
ord$model = rownames(ord)

if(is.null(bigord)) {
  bigord=ord
} else {
  bigord = merge(bigord,ord,all=T,by="model")
}



}
}
dev.off()
bigord
rownames(bigord) = bigord$model
rownames(bigord)[rownames(bigord)==""] = "NOTTESTED"
bigord=bigord[,2:ncol(bigord)]
colnames(bigord)[colnames(bigord)==""] = "NOTTESTED"
bigord=t(bigord)
 corrplot::corrplot(bigord,is.corr=F,method="color")
colMeans(bigord,na.rm=T)

smallord=bigord[,c(2:5)]
mins=sapply(1:nrow(smallord),FUN=function(x){
  row=smallord[x,]
  return(colnames(smallord)[which(row==min(row,na.rm=T))])
})
maxs=sapply(1:nrow(smallord),FUN=function(x){
  row=smallord[x,]
  return(colnames(smallord)[which(row==max(row,na.rm=T))])
})
}

if(doDXY ==T){
png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/images/dxy_boxplot_with_best_model.png",width=800)
par(mfrow=c(2,5),mar=c(4,4,1,0))
for (species in spplist) {
  
  #species="crissale"
shortspp = substr(species,1,3)  
print(species)
capspecies=toupper(species)
  
palette(c("lightblue","pink","lightgrey"))
  
  path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/"
  
  setwd(path)
  
  filename=list.files(path=path,
                      recursive = F,
                      pattern=paste("scaff_",shortspp,"_bothFST_DXY_ALL.txt$",sep=""))
  
  if(length(filename) != 0) {
    
    df = read.table(filename,
                    header=T,
                    numerals="no.loss",
                    as.is=F,stringsAsFactors = F)
    
    
  
    #levels(df$scafs) = c(1,"1A","1B",2:4,"4A",5:28,
    levels(df$chr) = c(1,"1A","1B",2:4,"4A",5:28,                    
                         "mtDNA","Z","LG2","LG5","LGE22")
    df=df[order(df$chr, df$chr),]
    #df=df[order(df$scafs, df$starts),]
    
    df$order = 1:nrow(df)
    
    #scafchrs=paste(df$scafs,df$starts)
    scafchrs=paste(df$chr,df$starts)
    
    if( length(scafchrs) == length(unique(scafchrs)) ) {
      nodups=df
    } else {
      nodups=df[-(which (scafchrs%in% (scafchrs[duplicated(scafchrs)]))),]
    }
    
    if(smallplot==T){
    png(paste(species,"_dxy_2020.png",sep=""),width=600)
    #plot(nodups$order,nodups$means,col=as.numeric(as.factor(nodups$scafs)))
    plot(nodups$order,nodups$dxymeans,col=as.numeric(as.factor(nodups$chr)))
    text(x=nodups$order[nodups$starts==1],
         y=c(min(nodups$dxymeans,na.rm=T),mean(nodups$dxymeans,na.rm=T),
             min(nodups$dxymeans,na.rm=T)+mean(nodups$dxymeans,na.rm=T)),
         labels=nodups$chr[nodups$starts==1],
         cex=0.75)
    #text(x=nodups$order[nodups$starts==1],
    #     y=c(min(nodups$means),mean(nodups$means),
    #         min(nodups$means)+mean(nodups$means)),
    #     labels=nodups$scafs[nodups$starts==1],
    #     cex=0.75)
    dev.off()
    }
    agg = aggregate(nodups$dxymeans~nodups$chr,FUN=function(x){mean(x,na.rm=T)})
    #agg = aggregate(nodups$means~nodups$scafs,FUN=function(x){mean(x,na.rm=T)})
    colnames(agg) = c("chrom","meanDxy")
    
    models=read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/bivariate/bivariate_gdm_results.csv",
                    skip=1)
    thismod=models[models$SPECIES==capspecies,]
    thismod=thismod[,c("NUMERICDATASET","MOSTA.1MODEL2")]
    colnames(thismod) = c("chrom","m2")
    thismod$m2[is.na(thismod$m2)] = ""
    
    agg=merge(agg,thismod,all=T)
    
    colors=RColorBrewer::brewer.pal(8,"Dark2")
    palette(c("white",colors[1],colors[3],colors[7],colors[4],"grey",
              colors[c(2,5,6,8)]))
    p=boxplot(agg$meanDxy~agg$m2,col=1:6,
              xlab=species,ylab="Mean Dxy",las=2)
    
    
  }
}
dev.off()
}

if(doBIG==T){
models=read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/bivariate/bivariate_gdm_results.csv",
                skip=1)
bigplot=read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.txt",
                 header=T)
missing=read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/MISSING/missing_data_summarized.txt",
                 sep="\t")
colnames(missing) = c("SPECIES","DATASET","MISSING")
missing$SPECIES = substr(missing$SPECIES,1,3)

agg = aggregate(cbind(bigplot$Fst,bigplot$dxymeans),by=list(bigplot$chr,bigplot$species),FUN=function(x){mean(x,na.rm=T)})
colnames(agg) = c("DATASET","SPECIES","MEANFST","MEANDXY")

genomeagg = aggregate(cbind(bigplot$Fst,bigplot$dxymeans),by=list(bigplot$species),FUN=function(x){mean(x,na.rm=T)})
colnames(genomeagg) = c("SPECIES","MEANFST","MEANDXY")
genomeagg$DATASET = "GENOME"

agg = merge(agg,genomeagg,all=T)
agg$SPECIES = toupper(agg$SPECIES)
agg$DATASET = toupper(agg$DATASET)
agg = unique(agg)

submod = models[,c("DATASET","SPECIES","MOSTA.1MODEL2","MOSTA.1MODEL1","MOSTA.1MODEL3")]
submod$DATASET=(substr(submod$DATASET,4,10)) 
submod$DATASET[submod$DATASET=="OME"] = "GENOME"
submod$NUMDAT = as.numeric(submod$DATASET)
submod$NUMDAT[is.na(submod$NUMDAT)] = submod$DATASET[is.na(submod$NUMDAT)]
submod$DATASET = submod$NUMDAT
submod=submod[,1:length(submod)-1]
submod$SPECIES = substr(submod$SPECIES,1,3)
submod$DATASET[submod$DATASET %in% c("04A")] = "4A"
submod$DATASET[submod$DATASET %in% c("01A")] = "1A"
submod$DATASET[submod$DATASET %in% c("01B")] = "1B"
submod = unique(submod)

together = merge(agg,submod,all=T)
together = merge(together,missing,all=T)

together$MOSTA.1MODEL2 = as.character(together$MOSTA.1MODEL2)
together$MOSTA.1MODEL2[together$MOSTA.1MODEL2==""] = "NONE"
together$MOSTA.1MODEL2[is.na(together$MOSTA.1MODEL2)] = "NONE"
together$MOSTA.1MODEL1 = as.character(together$MOSTA.1MODEL1)
together$MOSTA.1MODEL1[together$MOSTA.1MODEL1==""] = "NONE"
together$MOSTA.1MODEL1[is.na(together$MOSTA.1MODEL1)] = "NONE"
together$MOSTA.1MODEL3 = as.character(together$MOSTA.1MODEL3)
together$MOSTA.1MODEL3[together$MOSTA.1MODEL3==""] = "NONE"
together$MOSTA.1MODEL3[is.na(together$MOSTA.1MODEL3)] = "NONE"

together=unique(together)

colors=RColorBrewer::brewer.pal(8,"Dark2")
boxplot(together$MEANFST~together$MOSTA.1MODEL2,ylab="Mean Fst",las=2,xlab="Best Model",
        col=c(colors[1],colors[3],colors[7],colors[4],"grey","white",
              colors[c(2,5,6,8)]))
boxplot(together$MEANFST~together$MOSTA.1MODEL2,ylab="Mean Fst",las=2,ylim=c(0,0.15),
        xlab="Best Model",varwidth=T,notch=F,
        col=c(colors[1],colors[3],colors[7],colors[4],"grey","white",
              colors[c(2,5,6,8)]))
anova <- aov(together$MEANFST~together$MOSTA.1MODEL2
             +together$SPECIES
             )
summary.aov(anova)

anova1 <- aov(together$MEANFST~together$MOSTA.1MODEL1
             +together$SPECIES
)
summary.aov(anova1)
anova3 <- aov(together$MEANFST~together$MOSTA.1MODEL3
             +together$SPECIES
)
summary.aov(anova3)

values=TukeyHSD(anova)$`together$MOSTA.1MODEL2`
values = values[order(values[,4]),]
barplot(values[,4],las=2,ylim=c(0,0.05))

boxplot(together$MEANDXY~together$MOSTA.1MODEL2,ylab="Mean Dxy",las=2,xlab="Best Model",
        col=c(colors[1],colors[3],colors[7],colors[4],"grey","white",
              colors[c(2,5,6,8)]))
anova <- aov(together$MEANDXY~together$MOSTA.1MODEL2
             +together$SPECIES
)
summary.aov(anova) ## significant


anova1 <- aov(together$MEANDXY~together$MOSTA.1MODEL1
              +together$SPECIES
)
summary.aov(anova1)
anova3 <- aov(together$MEANDXY~together$MOSTA.1MODEL3
              +together$SPECIES
)
summary.aov(anova3)

values=TukeyHSD(anova)$`together$MOSTA.1MODEL2`
values = values[order(values[,4]),]
barplot(values[,4],las=2,ylim=c(0,0.05))

values=TukeyHSD(anova1)$`together$MOSTA.1MODEL1`
values = values[order(values[,4]),]
barplot(values[,4],las=2,ylim=c(0,0.05))

values=TukeyHSD(anova3)$`together$MOSTA.1MODEL3`
values = values[order(values[,4]),]
barplot(values[,4],las=2,ylim=c(0,0.05))


boxplot(together$MISSING~together$MOSTA.1MODEL2,ylab="Percent Missing",las=2,xlab="Best Model",
        col=c(colors[1],colors[3],colors[7],colors[4],"grey","white",
              colors[c(2,5,6,8)]))
anova <- aov(together$MISSING~together$MOSTA.1MODEL2
             #+together$SPECIES
)
summary.aov(anova) ## significant before and after species
values=TukeyHSD(anova)$`together$MOSTA.1MODEL2`
values = values[order(values[,4]),]
barplot(values[,4],las=2,ylim=c(0,0.05))

anova1 <- aov(together$MISSING~together$MOSTA.1MODEL1
             #+together$SPECIES
)
summary.aov(anova1) ## significant before and after species
values=TukeyHSD(anova1)$`together$MOSTA.1MODEL1`
values = values[order(values[,4]),]
barplot(values[,4],las=2,ylim=c(0,0.05))

anova3 <- aov(together$MISSING~together$MOSTA.1MODEL3
              #+together$SPECIES
)
summary.aov(anova3) ## significant before and after species
values=TukeyHSD(anova3)$`together$MOSTA.1MODEL3`
values = values[order(values[,4]),]
barplot(values[,4],las=2,ylim=c(0,0.05))

together_small = together[together$MOSTA.1MODEL2!="NONE",]

## final plot
png("missing_fst_dxy_allspecies.png",width=7,height=2,units = "in",
    res=300)
par(mfrow=c(1,3))
par(mar=c(4,4,0.5,0))
boxplot(together_small$MEANFST~together_small$MOSTA.1MODEL2,ylab="Mean Fst",las=2,xlab="",
        col=c(colors[1],colors[3],colors[7],colors[4],"grey","white",
              colors[c(2,5,6,8)]))
boxplot(together_small$MEANDXY~together_small$MOSTA.1MODEL2,ylab="Mean Dxy",las=2,xlab="Best Model",
        col=c(colors[1],colors[3],colors[7],colors[4],"grey","white",
              colors[c(2,5,6,8)]))
boxplot(together_small$MISSING~together_small$MOSTA.1MODEL2,ylab="Percent Missing",las=2,xlab="",
        col=c(colors[1],colors[3],colors[7],colors[4],"grey","white",
              colors[c(2,5,6,8)]))
dev.off()

}



