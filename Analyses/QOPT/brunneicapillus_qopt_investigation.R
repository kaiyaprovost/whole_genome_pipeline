filename = "/Users/kprovost/Dropbox (AMNH)/Dissertation/bru_longitude_cluster_interrogation.txt"
df = read.table(filename,sep="\t",header=T,strip.white = T,stringsAsFactors = F,fill=T)

df = df[order(df$LONG),]

plot(df$LONG,df$K2_A) ## would need that s-shaped graph
plot(df$LONG,df$P2_A) ## try remove those inds
plot(df$LONG,df$C2_A) ## 3 inds as outliers?

plot(df$P2_A,df$C2_A)

## which column is max
df_sub = df[,c("P3_A","P3_B","P3_C")]
barplot(t(as.matrix(df_sub)),beside=F)

plot(df$LONG,df$P3_A,col="red",pch=1,ylim=c(0,1))
points(df$LONG,df$P3_B,col="cyan",pch=2)
points(df$LONG,df$P3_C,col="black",pch=3)

df$P3_DIF_AB = (df$P3_A - df$P3_B) / (df$P3_A + df$P3_B)
df$P3_DIF_AC = (df$P3_A - df$P3_C) / (df$P3_A + df$P3_C)
df$P3_DIF_BC = (df$P3_B - df$P3_C) / (df$P3_B + df$P3_C)
plot(df$LONG,df$P3_DIF_AB,col="red",pch=1,ylim=c(-1,1)); abline(h=0)
plot(df$LONG,df$P3_DIF_AC,col="cyan",pch=1,ylim=c(-1,1)); abline(h=0)
plot(df$LONG,df$P3_DIF_BC,col="black",pch=1,ylim=c(-1,1)); abline(h=0)
