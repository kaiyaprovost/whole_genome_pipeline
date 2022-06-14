x = read.csv("/Users/kprovost/Documents/Dissertation/CHAPTER2_GENOMES/momi_test_forr.txt")

png("/Users/kprovost/Documents/Dissertation/CHAPTER2_GENOMES/momi_test_asym.png",h=800,w=800,u="px")
par(mfrow=c(2,2),mar=c(2,2,2,0),cex=2)
plot(x$divTimeLow,ylim=c(450000,2500000),type="l",lwd=1,main="Div Time",xaxt="n",xlab="")
axis(1,at=1:6,labels=x$priorset,las=1)
points(x$divTimeHigh,ylim=c(450000,2500000),type="l",lwd=1)
points(x$divTimeEst,ylim=c(450000,2500000),pch=15)

plot(x$extNElow,ylim=c(0,2000000),type="l",lwd=1,main="Ne",xaxt="n",xlab="")
axis(1,at=1:6,labels=x$priorset,las=1)
points(x$extNEhigh,ylim=c(0,2000000),type="l",lwd=1)
points(x$n_s,ylim=c(0,2000000),pch=15,col="green")
points(x$n_c,ylim=c(0,2000000),pch=1,col="blue")

plot(x$migLo,ylim=c(0,1),type="l",lwd=1,main="Mig Rate",xaxt="n",xlab="")
axis(1,at=1:6,labels=x$priorset,las=1)
points(x$migHi,ylim=c(0,1),type="l",lwd=1)
points(x$mig_s2c,ylim=c(0,1),pch=15,col="green")
points(x$mig_c2s,ylim=c(0,1),pch=1,col="blue")

plot(x$migtimeLow,ylim=c(0,1000000),type="l",lwd=1,main="Mig Time",xaxt="n",xlab="")
axis(1,at=1:6,labels=x$priorset,las=1)
points(x$migtimeHi,ylim=c(0,1000000),type="l",lwd=1)
points(x$migtime,ylim=c(0,1000000),pch=15)
dev.off()


par(mfrow=c(2,3))
plot(x$divTimeEst,x$migtime)

plot(x$divTimeEst,x$n_s,pch=15,col="green",ylim=c(0,750000))
points(x$divTimeEst,x$n_c,pch=1,col="blue",ylim=c(0,750000))

plot(x$divTimeEst,x$mig_s2c,pch=15,col="green",ylim=c(0,0.25))
points(x$divTimeEst,x$mig_c2s,pch=1,col="blue",ylim=c(0,0.25))

plot(x$n_s,x$n_c)

plot(x$n_s,x$mig_s2c,pch=15,col="green",ylim=c(0,0.25))
points(x$n_s,x$mig_c2s,pch=1,col="blue",ylim=c(0,0.25))

plot(x$n_c,x$mig_s2c,pch=15,col="green",ylim=c(0,0.25))
points(x$n_c,x$mig_c2s,pch=1,col="blue",ylim=c(0,0.25))


library(corrplot)
corrplot(cor(x[,c(5,8,9,12,13,16)]),method="ellipse",
         type="upper",diag=F)

corrplot(cor(x[2:16]),method="color",
         type="full",diag=T,order="hclust",
         addrect=3)

corrplot(cor(x[2:16])[c(2,3,5,6,9,10,13,14),c(1,4,7,8,11,12,15),drop=F],
         method="ellipse",type="full",diag=T)




#####

y = read.csv("/Users/kprovost/Documents/Dissertation/CHAPTER2_GENOMES/alignment_metrics.output_processed.txt",
             sep="\t",header=T,stringsAsFactors = F)
boxplot(y$PCT_PF_READS_ALIGNED~y$SPECIES,varwidth=T)

hist(y$PCT_PF_READS_ALIGNED[y$SPECIES=="fus"])
hist(y$PCT_PF_READS_ALIGNED,type="n")
