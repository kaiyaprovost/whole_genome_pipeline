library(RcppCNPy)
file="/Users/kprovost/Dropbox (AMNH)/Phainopepla-nitens-called.Phainopepla-nitens_PCAngsd.1.inbreed.npy"
x=npyLoad(file)
x
file2="/Users/kprovost/Dropbox (AMNH)/Phainopepla-nitens-called.Phainopepla-nitens_PCAngsd.2.inbreed.npy"
y=npyLoad(file2)
y
plot(x,y)
barplot(x)


file3="/Users/kprovost/Dropbox (AMNH)/Phainopepla-nitens-called.Phainopepla-nitens_PCAngsd.1.selection.npy"
z=npyLoad(file3)
head(z)
z2 = z[1:100000,]
z2p=(pchisq(z2, df=1, lower.tail=FALSE))
plot(z2,z2p)
plot(z2,cex=0.3)
abline(h=3.75,col="red")
abline(h=6.6,col="red")
abline(h=10.5,col="red")

lm(log10(z2p)~z2)

plot(z2,log10(z2p))
hist(z2p)

file4="/Users/kprovost/Dropbox (AMNH)/Phainopepla-nitens-called.Phainopepla-nitens_PCAngsd.2.selection.npy"
w=npyLoad(file4)
head(w)
w2 = w[1:1000,]

