path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/FST/3.slidingWindowJobs/"
setwd(path)

species = c("Amphispiza_bilineata",
            #"Auriparus_flaviceps",
            #"Campylorhynchus_brunneicapillus",
            "Cardinalis_sinuatus",
            "Melozone_fusca"#,
            #"Phainopepla_nitens",
            #"Polioptila_melanura",
            #"Toxostoma_crissale",
            #"Toxostoma_curvirostre"#,
            #"Vireo_bellii"
            )

for (spp in species) {
  print(spp)
  file = paste(path,"SON_CHI_",spp,"_FST_slidingwindow_chrfix.fst",sep="")
  
  lines = readLines(file)
  end = substr(lines[1],start=nchar(lines[1])-2,stop=nchar(lines[1]))
  
  if (end != "Fst") {
    newline1 = paste(lines[1],"\tFst",sep="")
    lines[1] = newline1
    writeLines(lines,file)
  }
  
  windows = read.csv(file,sep="\t",row.names=NULL)
  head(windows)
  
  windows = (windows[order(-(windows$Nsites)),])
  
  #plot(windows$Fst,col=as.numeric(windows$chr))
  
  final = data.frame(Fst=character(0),Chrom=character(0),Freq=character(0))
  
  for (i in 1:length(unique(windows$chr))) {
    print(paste(spp, i,"of",length(unique(windows$chr))))
    chrom = unique(windows$chr)[[i]]
    temp = windows[windows$chr==chrom,]
    frequ = nrow(temp)
    fsts = temp$Fst
    
    eyes = rep(chrom,frequ)
    freqs = rep(frequ,frequ)
    
    toadd = data.frame(Fst=fsts,Chrom=eyes,Freq=freqs)
    #colnames(toadd) = c("Fst","Chrom","Freq")
    
    final = rbind(toadd,final)
    
  }
  
  final=final[order(-final$Freq,final$Chrom),] 
  final$Chrom = factor(final$Chrom,levels(final$Chrom)[unique(final$Chrom)])
  
  head(final)
  
  lvls = levels(final$Chrom)
 
  palette(c("grey","black"))
  
  png(paste(spp,"_windows_chrfix.png",sep=""),
      width=700,height=350)
  par(bg=NA,col.axis="white",fg="white",
      col.lab="white",col.main="white")
  plot(final$Fst,col=as.numeric(as.factor(final$Chrom)),cex=0.2,
       main=spp,xlab="Window (Scaffold)",ylab="FST",
       ylim=c(0,1))
  dev.off()
  
  palette(c("red","orange","goldenrod","green","blue","purple",
            "cyan","grey","brown","magenta"))
  
  png(paste(spp,"_windows_chrfix_zoomed.png",sep=""),
      width=700,height=350)
  par(bg=NA,col.axis="white",fg="white",
      col.lab="white",col.main="white")
  plot(final$Fst,col=as.numeric(as.factor(final$Chrom)),cex=0.2,
       main=spp,xlab="Window (Scaffold)",ylab="FST")
  dev.off()
 
  
  shortchr = substr(final$Chrom ,1,8)
  only_chr = final[which(shortchr %in% "PseudoNC"),]
  
  chrnames = substr(unique(as.factor(only_chr$Chrom)),24,50)
  
  
  png(paste(spp,"_windows_chrfix_chrs.png",sep=""),
      width=700,height=350)
  plot(only_chr$Fst,col=as.numeric(as.factor(only_chr$Chrom)),cex=0.2,
       main=spp,xlab="Window (Scaffold)",ylab="FST",
       ylim=c(0,1))
  text(y=0.8,x=seq(1,100000,round(100000/length(chrnames))),labels=chrnames,cex=0.5)
  
  dev.off()
  
  
  png(paste(spp,"_windows_chrfix_chrs_zoon.png",sep=""),
      width=700,height=350)
  plot(only_chr$Fst,col=as.numeric(as.factor(only_chr$Chrom)),cex=0.2,
       main=spp,xlab="Window (Scaffold)",ylab="FST")
  text(y=mean(only_chr$Fst),x=seq(1,100000,round(100000/length(chrnames))),labels=chrnames,cex=0.5)
  

  dev.off()
   
}



testfile = "/Users/kprovost/Documents/Dissertation/CHAPTER2_GENOMES/ANALYSIS/FST/3.slidingWindowJobs/SON_CHI_Toxostoma_curvirostre_FST_slidingwindow.fst"
test = read.csv(testfile,sep="\t")

summary(test$Fst)
avg = mean(test$Fst)
stdev = sd(test$Fst)
test$zscore = (test$Fst - avg) / stdev
quan1 = quantile(test$Fst,0.05)
quan2 = quantile(test$Fst,0.95)
test$quantile = 1
test$quantile[test$Fst <= quan1] = 2
test$quantile[test$Fst >= quan2] = 3

#hist(test$zscore)
summary(test$zscore)
test$abszscore = abs(test$zscore)
#plot(test$Fst,test$zscore)

#plot(test$Fst,col="black")
plot(test$Fst,col=as.numeric(test$quantile))

