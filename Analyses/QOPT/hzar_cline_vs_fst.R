df = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/GDM_results/bivariate_gdm_results_useold.csv",
                sep=",",
                header=T,fill=T,
                stringsAsFactors = F,
                skip = 0)

genome = df[df$DATASET=="GENOME",]
FST=genome[,c("MEAN_FST_100","SPECIES")]
mw=aggregate(df$HZAR_width[df$DATASET!="GENOME"]~df$SPECIES[df$DATASET!="GENOME"],FUN=function(x){mean(x,na.rm=T)})
sw=aggregate(df$HZAR_width[df$DATASET!="GENOME"]~df$SPECIES[df$DATASET!="GENOME"],FUN=function(x){sd(x,na.rm=T)})
mc=aggregate(df$HZAR_center[df$DATASET!="GENOME"]~df$SPECIES[df$DATASET!="GENOME"],FUN=function(x){mean(x,na.rm=T)})
sc=aggregate(df$HZAR_center[df$DATASET!="GENOME"]~df$SPECIES[df$DATASET!="GENOME"],FUN=function(x){sd(x,na.rm=T)})

colnames(mw)=c("SPECIES","MW")
colnames(mc)=c("SPECIES","MC")
colnames(sw)=c("SPECIES","SW")
colnames(sc)=c("SPECIES","SC")

newdf = merge(merge(merge(merge(FST,mw),sw),mc),sc)

## spp is alphabetical

## NEED TO UPDATE THESE NUMS
#mw=c(15.84,7.51,15.32,13.77,7.62,13.28,15.89,10.97,6.94,10.28)
#sw=c(0.72,4.04,0.52,1.96,3.12,1.82,1.2,0.154,1.78,3.2)
#mc=c(11.63,7.95,12.7,7.45,6.55,3.58,8.91,7.96,7.61,9.45)
#sc=c(3.93,1.07,4.82,2.72,0.49,2.81,2.64,0.9,0.59,0.96)
#FST=c(0.02,0.05,0.03,0.03,0.04,0.03,0.02,0.04,0.1,0.06)

pdf("hzar_cline_vs_fst_pval_6may2022.pdf")
par(mfrow=c(2,2),mar=c(4,4,1.5,1.5))
plot(newdf$MEAN_FST_100,newdf$MW,ylab="Mean Cline Width",ylim=c(5,17),xlim=c(0.01,0.13),
     type="n",col="red",xlab="FST")
points(newdf$MEAN_FST,newdf$MW,col="black",pch=0)
#mw_mod=lm(mw~FST)
#abline(mw_mod)
mw_mod2=lm(newdf$MW~newdf$MEAN_FST_100)
abline(mw_mod2,col="red")
## BOTH ARE SIGNIFICANT 
mtext(paste("adjR2 =",round(summary(mw_mod2)$adj.r.squared,2)))
#mtext(paste("R2 =",round(summary(mw_mod)$r.squared,2)))
#mtext(paste("p =",round(summary(mw_mod)$coefficients[2,4],5)))

plot(newdf$MEAN_FST_100,newdf$SW,ylab="S.D. Cline Width",col="red",xlim=c(0.01,0.13), type="n",xlab="FST")
points(newdf$MEAN_FST_100,newdf$SW,ylab="S.D. Cline Width",col="black",pch=0)
sw_mod2=lm(newdf$SW~newdf$MEAN_FST_100)
#abline(sw_mod)
#sw_mod2=lm(newdf$SW~newdf$MEAN_FST)
## BOTH ARE NOT SIGNIFICANT 
abline(sw_mod2,col="red")
mtext(paste("adjR2 =",round(summary(sw_mod)$adj.r.squared,2)))
#mtext(paste("R2 =",round(summary(sw_mod)$r.squared,2)))
#mtext(paste("p =",round(summary(sw_mod)$coefficients[2,4],5)))

plot(newdf$MEAN_FST_100,newdf$MC,ylab="Mean Cline Center",col="red",xlim=c(0.01,0.13), type="n",xlab="FST")
points(newdf$MEAN_FST_100,newdf$MC,col="black",pch=0)
#mc_mod=lm(mc~FST)
#abline(mc_mod)
mc_mod2=lm(newdf$MC~newdf$MEAN_FST_100)
abline(mc_mod2,col="red")
## BOTH ARE NOT SIGNIFICANT 
mtext(paste("adjR2 =",round(summary(mc_mod2)$adj.r.squared,2)))
#mtext(paste("R2 =",round(summary(mc_mod)$r.squared,2)))
#mtext(paste("p =",round(summary(mc_mod)$coefficients[2,4],5)))

plot(newdf$MEAN_FST_100,newdf$SC,ylab="S.D. Cline Center",col="red",xlim=c(0.01,0.13),type="n",xlab="FST")
points(newdf$MEAN_FST_100,newdf$SC,col="black",pch=0)
#sc_mod=lm(sc~FST)
#abline(sc_mod)
## GOES FROM SIG TO NEARLY SIG
sc_mod2=lm(newdf$SC~newdf$MEAN_FST_100)
abline(sc_mod2,col="red")
mtext(paste("adjR2 =",round(summary(sc_mod2)$adj.r.squared,2)))
#mtext(paste("R2 =",round(summary(sc_mod)$r.squared,2)))
#mtext(paste("p =",round(summary(sc_mod)$coefficients[2,4],5)))

dev.off()

## add in clines and recombination vs chromosome length 

recom = df[,c("HZAR_width","CHROM_LENGTH","MEAN_RECOMB_100","SPECIES")]
#recom = recom[recom$CHROM_LENGTH!=max(df$CHROM_LENGTH,na.rm=T),]

cl=aggregate(df$CHROM_LENGTH~df$DATASET,FUN=function(x){mean(x,na.rm=T)})
mw=aggregate(df$HZAR_width~df$DATASET,FUN=function(x){mean(x,na.rm=T)})
sw=aggregate(df$HZAR_width~df$DATASET,FUN=function(x){sd(x,na.rm=T)})
mr=aggregate(df$MEAN_RECOMB_100~df$DATASET,FUN=function(x){mean(x,na.rm=T)})
sr=aggregate(df$MEAN_RECOMB_100~df$DATASET,FUN=function(x){sd(x,na.rm=T)})

colnames(cl)=c("DATASET","CL")
colnames(mw)=c("DATASET","MW")
colnames(mr)=c("DATASET","MR")
colnames(sw)=c("DATASET","SW")
colnames(sr)=c("DATASET","SR")

newdf2 = merge(merge(merge(merge(cl,mw),sw),mr),sr)

pdf("hzar_cline_recom_vs_chromlength_pval_6may2022.pdf")
par(mfrow=c(2,2),mar=c(4,4,1.5,1.5))
plot(log10(newdf2$CL),newdf2$MR,ylab="Mean Recombination Rate",col="red",type="n",xlab="Log Chromosome Length")
points(log10(newdf2$CL),newdf2$MR,col="black",pch=0)
mr_mod2=lm(newdf2$MR~log10(newdf2$CL))
abline(mr_mod2,col="red")
mtext(paste("adjR2 =",round(summary(mr_mod2)$adj.r.squared,2)))

plot(log10(newdf2$CL),newdf2$SR,ylab="S.D. Recombination Rate",col="red",type="n",xlab="Log Chromosome Length")
points(log10(newdf2$CL),newdf2$SR,col="black",pch=0)
sr_mod2=lm(newdf2$SR~log10(newdf2$CL))
abline(sr_mod2,col="red")
mtext(paste("adjR2 =",round(summary(sr_mod2)$adj.r.squared,2)))

plot(log10(newdf2$CL),newdf2$MW,ylab="Mean Cline Width",col="red",type="n",xlab="Log Chromosome Length")
points(log10(newdf2$CL),newdf2$MW,col="black",pch=0)
mw_mod2=lm(newdf2$MW~log10(newdf2$CL))
abline(mw_mod2,col="red")
mtext(paste("adjR2 =",round(summary(mw_mod2)$adj.r.squared,2)))

plot(log10(newdf2$CL),newdf2$SW,ylab="S.D. Cline Width",col="red",type="n",xlab="Log Chromosome Length")
points(log10(newdf2$CL),newdf2$SW,col="black",pch=0)
sw_mod2=lm(newdf2$SW~log10(newdf2$CL))
abline(sw_mod2,col="red")
mtext(paste("adjR2 =",round(summary(sw_mod2)$adj.r.squared,2)))
dev.off()
