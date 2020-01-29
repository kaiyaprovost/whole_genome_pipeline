file = "/Users/kprovost/Documents/Dissertation/CHAPTER2_GENOMES/status_of_tissues.txt"

table = read.csv(file,sep="\t")
summary(table)

needed = table[table$NEED.SEQ==1,]
summary(needed)

library(reshape)

barplot(needed$STATUS)

mt = melt(needed)
ct = cast(mt, STATUS~SCIENTIFIC)

mct = as.matrix(ct)/2
colnames(mct) = c("A. bil.","A. fla.","C. bru.",
                  "C. sin.","M. fus.","P. nit.",
                  "P. mel.","T. cri.","T. cur.","V. bel.")

par(bg=NA,col="white",col.axis="white",col.lab="white",
    col.main="white",col.sub="white")
barplot(mct, main="Species and Status",las=2,
        col=NA,bor=NA,plot=T,axes=F)
abline(h=0,col="white",lty=1)
abline(h=5,col="white",lty=1)
abline(h=10,col="white",lty=1)
abline(h=15,col="white",lty=1)
abline(h=20,col="white",lty=1)
barplot(mct, main="Species and Status",las=2,
        col=c("blue","red","yellow","grey"),bor="white",add=T)
