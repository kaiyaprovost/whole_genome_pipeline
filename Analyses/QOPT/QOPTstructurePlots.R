## This script makes Figure 4 and Appendix Figure 1, the structure plots.

## import color scheme
library(RColorBrewer)
col = brewer.pal(12,"Paired")
red=col[6]
blue=col[2]
green=col[3]
yellow=col[7]
grey=col[9]

## Appendix Figure 1 Steps

## Import Data
## Note that these must be sorted by the column "SORTORDER" for this to work properly
structure_all <- read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/AllSpeciesMetadata_allK_9june2020.csv",
                          row.names=NULL)
names(structure_all)

structure_all = structure_all[!(is.na(structure_all$K2_A)),]

## remove if no K2 values


## pull out the relevant data
metadata_cols = 1:12
k2_cols = 13:14
k3_cols = 15:17
k4_cols = 18:21
k5_cols = 22:26
p2_cols = 27:28
label_cols = 5:6

pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/teststructure.pdf")

for (s in 1:length(unique(structure_all$SP))) {
  
  spp = unique(structure_all$SP)[s]
  print(spp)
  
  structure_big = structure_all[structure_all$SP==spp,]
  
  ## subset the clusters of interest with the metadata
  k2 = structure_big[,c(metadata_cols,k2_cols,label_cols = 6:7)] ## clustering for K2 with metadata
  k3 = structure_big[,c(metadata_cols,k3_cols,label_cols = 6:7)] ## clustering for K3 with metadata
  k4 = structure_big[,c(metadata_cols,k4_cols,label_cols = 6:7)] ## clustering for K4 with metadata
  k5 = structure_big[,c(metadata_cols,k5_cols,label_cols = 6:7)] ## clustering for K5 with metadata
  
  ## subset the clusters of interest without metadata, and 
  ## transpose for easier processing later 
  k2m = t(structure_big[,c(k2_cols)]) ## just K2 clusters, and transposed for ease
  k3m = t(structure_big[,c(k3_cols)]) ## just K3 clusters, and transposed for ease
  k4m = t(structure_big[,c(k4_cols)]) ## just K4 clusters, and transposed for ease
  k5m = t(structure_big[,c(k5_cols)]) ## just K5 clusters, and transposed for ease
  
  
  
  ## plot Appendix Figure 1
  #pdf("AppendixFigure1_structure.pdf",
  #    height=4,width=4,pointsize=10)
  
  par(mfrow=c(4,1),mar=c(0.5,0,0,0),lwd=0.1)
  two = barplot(k2m,col=c(green,blue),axisnames=F,axes=F,xpd=F)
  mtext(text = k2$LABEL, side = 1, at = two, line = 0.1)
  three = barplot(k3m,col=c(green,blue,red),axisnames=F,axes=F,xpd=F)
  mtext(text = k3$LABEL, side = 1, at = two, line = 0.1)
  four = barplot(k4m,col=c(green,blue,red,yellow),axisnames=F,axes=F,xpd=F)
  mtext(text = k4$LABEL, side = 1, at = two, line = 0.1)
  five = barplot(k5m,col=c(green,blue,red,yellow,grey),axisnames=F,axes=F,xpd=F)
  mtext(text = k5$LABEL, side = 1, at = two, line = 0.1)
  
}

dev.off()


pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/teststructureP1P2.pdf")

structure_all = structure_all[!(is.na(structure_all$P1)),]


for (s in 1:length(unique(structure_all$SP))) {
  
  spp = unique(structure_all$SP)[s]
  print(spp)
  
  structure_big = structure_all[structure_all$SP==spp,]
  
  ## subset the clusters of interest with the metadata
  p2 = structure_big[,c(metadata_cols,p2_cols,label_cols = 6:7)] ## clustering for K2 with metadata
  ## subset the clusters of interest without metadata, and 
  ## transpose for easier processing later 
  p2m = t(structure_big[,c(p2_cols)]) ## just K2 clusters, and transposed for ease
  
  par(mfrow=c(1,1),mar=c(0.5,0,0,0),lwd=0.1)
  two = barplot(p2m,col=c(green,blue),axisnames=F,axes=F,xpd=F)
  mtext(text = p2$LABEL, side = 1, at = two, line = 0.1)

  }

dev.off()


## Figure 5 Steps
## Note that the variable names are all reused here
## Import data
structure <- read.csv("cardinalis_bothdatasets_forplotting_ingroup.csv",
                      row.names=1)
names(structure)

## pull out the relevant data
metadata_cols = 1:6
k2_cols = 7:8
k3_cols = 9:11
k4_cols = 12:15
k5_cols = 16:20
label_cols = 21:22

## subset the clusters of interest with the metadata
k2 = structure[,c(metadata_cols,k2_cols,label_cols = 21:22)] ## clustering for K2 with metadata
k3 = structure[,c(metadata_cols,k3_cols,label_cols = 21:22)] ## clustering for K3 with metadata
k4 = structure[,c(metadata_cols,k4_cols,label_cols = 21:22)] ## clustering for K4 with metadata
k5 = structure[,c(metadata_cols,k5_cols,label_cols = 21:22)] ## clustering for K5 with metadata

## subset the clusters of interest without metadata, and 
## transpose for easier processing later 
k2m = t(structure[,c(k2_cols)]) ## just K2 clusters, and transposed for ease
k3m = t(structure[,c(k3_cols)]) ## just K3 clusters, and transposed for ease
k4m = t(structure[,c(k4_cols)]) ## just K4 clusters, and transposed for ease
k5m = t(structure[,c(k5_cols)]) ## just K5 clusters, and transposed for ease

## Plot Figure 5
## Note post-processing in another image program is needed for the final labels
pdf("Figure 4 - structure.pdf",
    height=4,width=4,pointsize=10)
#png("/Users/kprovost/Dropbox (AMNH)/Cardinalis MS/FOR_DRYAD/data/str_files/Figure 4 - structure.png",
#    height=11,width=11,res=600,units="cm",pointsize=10)
par(mfrow=c(4,1),mar=c(0.5,0,0,0),lwd=0.1)
barplot(k2m,col=c(green,blue),axisnames=F,axes=F,
        #space=c(rep(0,1),1,rep(0,3-1),1,rep(0,3-1),1,rep(0,33-1),1,
        #        rep(0,11-1),1,rep(0,6-1),1,rep(0,2-1),1,rep(0,28-1),1,
        #        rep(0,1-1)),xpd=F)
        space=c(rep(0,2),0,rep(0,33-1),0,
                rep(0,11-1),0,rep(0,6-1),1,rep(0,2-1),0,rep(0,27-1),0,
                rep(0,1-1)),xpd=F)
barplot(k3m,col=c(green,blue,red),axisnames=F,axes=F,
        space=c(rep(0,2),0,rep(0,33-1),0,
                rep(0,11-1),0,rep(0,6-1),1,rep(0,2-1),0,rep(0,27-1),0,
                rep(0,1-1)),xpd=F)
#space=c(rep(0,1),1,rep(0,3-1),1,rep(0,3-1),1,rep(0,33-1),1,
#        rep(0,11-1),1,rep(0,6-1),1,rep(0,2-1),1,rep(0,28-1),1,
#        rep(0,1-1)))
barplot(k4m,col=c(green,blue,red,yellow),axisnames=F,axes=F,
        space=c(rep(0,2),0,rep(0,33-1),0,
                rep(0,11-1),0,rep(0,6-1),1,rep(0,2-1),0,rep(0,27-1),0,
                rep(0,1-1)),xpd=F)
#space=c(rep(0,1),1,rep(0,3-1),1,rep(0,3-1),1,rep(0,33-1),1,
#        rep(0,11-1),1,rep(0,6-1),1,rep(0,2-1),1,rep(0,28-1),1,
#        rep(0,1-1)))
barplot(k5m,col=c(green,blue,red,yellow,grey),axisnames=F,axes=F,
        space=c(rep(0,2),0,rep(0,33-1),0,
                rep(0,11-1),0,rep(0,6-1),1,rep(0,2-1),0,rep(0,27-1),0,
                rep(0,1-1)),xpd=F)
#        space=c(rep(0,1),1,rep(0,3-1),1,rep(0,3-1),1,rep(0,33-1),1,
#                rep(0,11-1),1,rep(0,6-1),1,rep(0,2-1),1,rep(0,28-1),1,
#                rep(0,1-1)),xpd=F)

dev.off()

## OPTIONAL STEPS

## You may also want to generate these plots by themselves, here is the code to do that and have them look nice
## This also includes the spacing between localities

png("/Users/kprovost/Dropbox (AMNH)/Cardinalis MS/FOR_DRYAD/data/str_files/Figure4_k2only.png",
    height=3,width=11,res=600,units="cm",pointsize=10)
par(mfrow=c(1,1),mar=c(0.5,0,0,0),lwd=0.1)
barplot(k2m,col=c(green,blue),axisnames=F,axes=F,
        space=c(rep(0,2),1,rep(0,33-1),1,
                rep(0,11-1),1,rep(0,6-1),1,rep(0,2-1),1,rep(0,27-1),1,
                rep(0,1-1)),xpd=F)
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Cardinalis MS/FOR_DRYAD/data/str_files/Figure4_k3only.png",
    height=3,width=11,res=600,units="cm",pointsize=10)
par(mfrow=c(1,1),mar=c(0.5,0,0,0),lwd=0.1)
barplot(k3m,col=c(green,blue,red),axisnames=F,axes=F,
        space=c(rep(0,2),1,rep(0,33-1),1,
                rep(0,11-1),1,rep(0,6-1),1,rep(0,2-1),1,rep(0,27-1),1,
                rep(0,1-1)),xpd=F)
dev.off()

png("/Users/kprovost/Dropbox (AMNH)/Cardinalis MS/FOR_DRYAD/data/str_files/Figure4_k4only.png",
    height=3,width=11,res=600,units="cm",pointsize=10)
par(mfrow=c(1,1),mar=c(0.5,0,0,0),lwd=0.1)
barplot(k4m,col=c(green,blue,red,yellow),axisnames=F,axes=F,
        space=c(rep(0,2),1,rep(0,33-1),1,
                rep(0,11-1),1,rep(0,6-1),1,rep(0,2-1),1,rep(0,27-1),1,
                rep(0,1-1)),xpd=F)
dev.off()



png("/Users/kprovost/Dropbox (AMNH)/Cardinalis MS/FOR_DRYAD/data/str_files/Figure4_k5only.png",
    height=3,width=11,res=600,units="cm",pointsize=10)
par(mfrow=c(1,1),mar=c(0.5,0,0,0),lwd=0.1)
barplot(k5m,col=c(green,blue,red,yellow,grey),axisnames=F,axes=F,
        space=c(rep(0,2),1,rep(0,33-1),1,
                rep(0,11-1),1,rep(0,6-1),1,rep(0,2-1),1,rep(0,27-1),1,
                rep(0,1-1)),xpd=F)
dev.off()

