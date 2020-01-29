file="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/Vireo-bellii-called.geno.vcf.lmiss.lmiss"

missing = read.table(file,sep="\t",header=T,stringsAsFactors = F)
#colnames(missing) = c("CHR","POS","N_DATA","N_GENOTYPE_FILTERED","N_MISS","F_MISS")
#missing = missing[-1,]

chroms = unique(missing$CHR)
keep = chroms[which(substr(chroms,1,8)=="PseudoNC")]
onlychrom = missing[which(missing$CHR %in% keep),]

## generate windows of missing data
windowsize = 100000
offset = 10000

windows = c()

for (chrom in keep) {
  print(chrom)
  thischrom = onlychrom[onlychrom$CHR==chrom,]
  finalpos = max(thischrom$POS)
  starts = seq(1,finalpos,windowsize)
  stops = (starts+windowsize)-1
  stops = unique(c(stops[stops<finalpos],finalpos))
  starts = starts[1:length(stops)]

  for(i in 1:length(stops)) {
    
    if (i %% 100 == 0 ) { 
      print(paste(stop,"/",finalpos)) 
      plot(windows[,6],col=as.numeric(as.factor(windows[,1])),
           ylim=c(0,1),ylab="Mean Missing Data",xlab="Scaffold")
    }
    start = starts[i]
    stop = stops[i]
    
    middle = round(mean(c(start,stop)))
    meta = cbind(chrom=chrom,middle=middle)
    
    window1 = thischrom[thischrom$POS >= start,]
    window1 = window1[window1$POS <= stop,]
    cols = colnames(window1)[3:6]
    window1 = as.matrix(window1[,3:6])
    windowmean = colMeans(window1)
    names(windowmean) = paste("mean",cols,sep="_")
    windowsd = rbind(matrixStats::colSds(window1))
    colnames(windowsd) = paste("sd",cols,sep="_")
    
    toadd = c(meta,windowmean,windowsd)
    names = c(colnames(meta),names(windowmean),colnames(windowsd))
    names(toadd) = names
    
    if(is.null(windows)) {
      windows = rbind(toadd)
    } else {
      windows = rbind(windows,toadd)
    }
    
  }
  
  #plot(windows[,2],windows[,6],col=as.numeric(as.factor(windows[,1])))
  
  
}

write.csv(windows,paste(file,".windows.temp",sep=""))

#sub = missing[1:100000,]

palette(  c(    "red",    "cyan",    "goldenrod",    
                "green",    "blue",    "purple",    
                "blue",    "black",    "brown",    "magenta"  ))

png(paste(file,".windows.png",sep=""),width = 700, height = 350)
plot(windows[,6],col=as.numeric(as.factor(windows[,1])),
     ylab="Mean Missing Data",xlab="Scaffold",cex=0.1)
dev.off()

plot(windows[,6],col=as.numeric(as.factor(windows[,1])),
     ylab="Mean Missing Data",xlab="Scaffold",cex=0.1)
points(as.numeric(windows[,10])+as.numeric(windows[,6]),col=as.numeric(as.factor(windows[,1])),
       cex=0.05)
points(as.numeric(windows[,6])-as.numeric(windows[,10]),col=as.numeric(as.factor(windows[,1])),
       cex=0.05)


### singletons

singlefile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/Vireo-bellii-called.geno.vcf.singletons.singletons"
single = read.table(singlefile,sep="\t",header=T,stringsAsFactors = F)
chroms = unique(single$CHROM)
keep = chroms[which(substr(chroms,1,8)=="PseudoNC")]
onlychrom2 = single[which(single$CHROM %in% keep),]
df = table(onlychrom2$CHROM,onlychrom2$INDV)

nrow(onlychrom) / nrow(onlychrom2)

onlychrom1.1 = onlychrom[,1:2]
onlychrom2.1 = onlychrom2[,1:2]
colnames(onlychrom2.1)[1] = "CHR"
onlychrom1.1$SINGLETON = 0
onlychrom2.1$SINGLETON = 1
mix=rbind(onlychrom1.1,onlychrom2.1)
#mix=merge(onlychrom1.1,onlychrom2.1,sort=F)
mix = mix[order(mix$CHR, mix$POS),]

agg = aggregate(mix$SINGLETON ~ mix$CHR + mix$POS, FUN=max)
names(agg) = names(mix)

plot(agg$SINGLETON,col=as.numeric(as.factor(agg$CHR)))


write.table(df,paste(singlefile,".byind.temp",sep=""))
