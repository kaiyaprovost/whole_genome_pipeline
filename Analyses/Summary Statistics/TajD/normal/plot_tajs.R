## values for the pestPG are sums for the windows. the theta values you see that are negative are log values. 
## to get an average for a window, need to divide by window size -- does one for each BP

specieslist=c(
#"Amphispiza-bilineata-SON",
#"Amphispiza-bilineata",
#"Auriparus-flaviceps-CHI",
#"Auriparus-flaviceps-SON",
#"Auriparus-flaviceps",
#"Campylorhynchus-brunneicapillus",
"Cardinalis-sinuatus-CHI",
"Cardinalis-sinuatus-SON",
"Cardinalis-sinuatus"#,
#"Melozone-fusca-CHI",
#"Melozone-fusca-SON",
#"Melozone-fusca"#,
#"Phainopepla-nitens",
#"Polioptila-melanura",
#"Toxostoma-crissale",
#"Toxostoma-curvirostre"#,
#"Vireo-bellii-CHI",
#"Vireo-bellii-NOWEIRD",
#"Vireo-bellii-SON",
#"Vireo-bellii"
)

## whole chromosomes at once
for (speciesname in specieslist) {
  
  print(speciesname)
  shortspp=substr(strsplit(speciesname,"-")[[1]][2],1,3)
  
  file = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas/",speciesname,"-taj2.thetas.idx_CHRFIX.pestPG",sep="")
  if(!(file.exists(file))){
    print("NOT EXISTS")
  } else {
    
    tajs = read.table(file,header=T)
    print(head(tajs))
    print(names(tajs)[2:length(tajs)])
    
    tajs = (tajs[order(tajs$Chr,tajs$WinCenter),])
    
    tajs$species=shortspp
    
    
    
    #bigplot = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.txt")
    #bigplotorder = unique(bigplot$chr[bigplot$species==shortspp])
    
    #tajs = tajs[order(match(tajs$Chr, bigplotorder)),]
    
    palette(  c(    "red",    "cyan",    "goldenrod",    "green",    "blue",    "purple",    "blue",    "black",    "brown",    "magenta"  ))
    

    for (colnam in names(tajs)[1:length(tajs)]) {
      print(colnam)
      png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas/",speciesname,"_CHROM_",colnam,".june2020.png",sep=""),width=800,height=300)
      
      toplot=as.numeric(as.character(tajs[,colnam]))
      
      if(sum(!is.na(toplot))!=0) {
        
        barplot(as.numeric(as.character(tajs[,colnam])),
                col=as.numeric(as.factor(tajs$Chr)),
                xlab="Chromosome",
                names=as.character(as.factor(tajs$Chr)),
                las=2)
        
        dev.off()
      }
    }
    
    ## merge tajs with bigplot eventually 
    
    
  }
}



bigtaj = NULL

files = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas",
                   pattern="taj2.thetasWindow.gz_CHRFIX.pestPG",
                   full.names = T,recursive = T)

for (speciesname in specieslist) {
  
  print(speciesname)
  shortspp=substr(strsplit(speciesname,"-")[[1]][2],1,3)
  
  file = files[grepl(speciesname,files)]
  #file = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas/",speciesname,"-taj2.thetasWindow.gz_CHRFIX.pestPG",sep="")
  if(length(file)<=0){
  #if(!(file.exists(file))){
  print("NOT EXISTS")
  } else {
  file=file[1]
  tajs = read.table(file,header=T)
  print(head(tajs))
  print(names(tajs)[2:length(tajs)])
  
  tajs = (tajs[order(tajs$Chr,tajs$WinCenter),])
  
  tajs$species=shortspp
  
  #bigplot = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.txt")
  #bigplotorder = unique(bigplot$chr[bigplot$species==shortspp])
  
  #tajs = tajs[order(match(tajs$Chr, bigplotorder)),]
  
  palette(  c(    "red",    "cyan",    "goldenrod",    "green",    "blue",    "purple",    "blue",    "black",    "brown",    "magenta"  ))
  
  for (colnam in names(tajs)[1:length(tajs)]) {
    print(colnam)
    png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas/",speciesname,colnam,".june2020.png",sep=""),width=800,height=300)
    
    toplot=as.numeric(as.character(tajs[,colnam]))
    
    if(sum(!is.na(toplot))!=0) {
      plot(as.numeric(as.character(tajs[,colnam])),col=as.numeric(as.factor(tajs$Chr)),cex=0.2,
           main=speciesname,xlab="Window (Scaffold)",ylab=colnam)
      
      if (colnam == "Tajima") {
        abline(h=2,col="grey")
        abline(h=-2,col="grey")
      }
      
      dev.off()
    }
  }
  
  ## merge tajs with bigplot eventually 
  
  if(is.null(bigtaj)){
    bigtaj = tajs
  } else {
    bigtaj = rbind(bigtaj,tajs)
  }
}
}

write.csv(bigtaj,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigtaj.new.june2020.temp",row.names=F)

bigtaj = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigtaj.new.june2020.txt")
colnames(bigtaj) = c("windowinfo",
                     "chr",
                     "midPos",
                     "tW",
                     "tP-pi",
                     "tF",
                     "tH",
                     "tL",
                     "Tajima",
                     "fuf",
                     "fud",
                     "fayh",
                     "zeng",
                     "nSites",
                     "species"
)

bigplot = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.june2020.txt")

bigplot2 = merge(bigplot, bigtaj, by = c("chr","midPos","species"),all=T)

write.csv(bigplot2,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot2.june2020.txt",row.names=F)
write.csv(bigplot2,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot2.june2020.txt",row.names=F)


bigplot2 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot2.june2020.txt")

names(bigplot2)
mu=2.21e-9
##pi=4Ne*mu so so pi/4mu = 'Ne'
bigplot2$average_theta_divby_windowsize = rowMeans(bigplot2[,c("tW","tH","tF","tP-pi","tL")]) /50000
bigplot2$est_ne_221e9 = (bigplot2$average_theta_divby_windowsize) / (4*mu)
names(bigplot2)

#write.csv(bigplot2,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot2.txt",row.names=F)

th_calc = (bigplot2$tH / 50000) / (4*mu)
tp_calc = (bigplot2$`tP-pi` / 50000) / (4*mu)
tl_calc = (bigplot2$tL / 50000) / (4*mu)
tf_calc = (bigplot2$tF / 50000) / (4*mu)
tw_calc = (bigplot2$tW / 50000) / (4*mu)


setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas/")

sizes = aggregate(bigplot2$midPos~bigplot2$chr,FUN=max)
names(sizes) = c("chr","midPos")
sizes = unique(sizes[order(sizes$midPos,sizes$chr,decreasing=T),])

bigplotorder = unique(bigplot2$chr)
tajs = bigplot2[order(bigplot2$midPos,decreasing=T),]
nrow(tajs)
tajs = tajs[order(match(tajs$chr,sizes$chr)),]
nrow(tajs)

for (s in unique(bigplot2$species)) {
  print(s)
  
  if(sum(complete.cases(unique(tajs$est_ne_221e9[tajs$species==s])))>0){
  
  
  palette(  c(    "red",    "cyan",    "goldenrod",    "green",    "blue",    "purple",    "blue",    "black",    "brown",    "magenta"  ))
  png(paste(s,"est_ne_221e9_averageThetasc.june2020.png",sep=""),width=800,height=300)
  plot(tajs$est_ne_221e9[tajs$species==s],col=as.numeric(as.factor(tajs$chr[tajs$species==s])),cex=0.2,
       main=paste(s,"average value Ne_est:",
                  formatC(mean(tajs$est_ne_221e9[tajs$species==s],na.rm=T),format="e",digits=2),sep=" "),
       xlab="Window (Scaffold)",ylab="average theta scaled by 50000")
  dev.off()
  
  png(paste(s,"est_ne_221e9_th_calc.june2020.png",sep=""),width=800,height=300)
  plot(th_calc[tajs$species==s],col=as.numeric(as.factor(tajs$chr[tajs$species==s])),cex=0.2,
       main=paste(s,"average value Ne_est:",
                  formatC(mean(th_calc[tajs$species==s],na.rm=T),format="e",digits=2),sep=" "),
       xlab="Window (Scaffold)",ylab="average theta scaled by 50000")
  dev.off()
  
  png(paste(s,"est_ne_221e9_tp_calc.june2020.png",sep=""),width=800,height=300)
  plot(tp_calc[tajs$species==s],col=as.numeric(as.factor(tajs$chr[tajs$species==s])),cex=0.2,
       main=paste(s,"average value Ne_est:",
                  formatC(mean(tp_calc[tajs$species==s],na.rm=T),format="e",digits=2),sep=" "),
       xlab="Window (Scaffold)",ylab="average theta scaled by 50000")
  dev.off()
  
  png(paste(s,"est_ne_221e9_tl_calc.june2020.png",sep=""),width=800,height=300)
  plot(tl_calc[tajs$species==s],col=as.numeric(as.factor(tajs$chr[tajs$species==s])),cex=0.2,
       main=paste(s,"average value Ne_est:",
                  formatC(mean(tl_calc[tajs$species==s],na.rm=T),format="e",digits=2),sep=" "),
       xlab="Window (Scaffold)",ylab="average theta scaled by 50000")
  dev.off()
  
  png(paste(s,"est_ne_221e9_tf_calc.june2020.png",sep=""),width=800,height=300)
  plot(tf_calc[tajs$species==s],col=as.numeric(as.factor(tajs$chr[tajs$species==s])),cex=0.2,
       main=paste(s,"average value Ne_est:",
                  formatC(mean(tf_calc[tajs$species==s],na.rm=T),format="e",digits=2),sep=" "),
       xlab="Window (Scaffold)",ylab="average theta scaled by 50000")
  dev.off()
  
  png(paste(s,"est_ne_221e9_tw_calc.june2020.png",sep=""),width=800,height=300)
  plot(tw_calc[tajs$species==s],col=as.numeric(as.factor(tajs$chr[tajs$species==s])),cex=0.2,
       main=paste(s,"average value Ne_est:",
                  formatC(mean(tw_calc[tajs$species==s],na.rm=T),format="e",digits=2),sep=" "),
       xlab="Window (Scaffold)",ylab="average theta scaled by 50000")
  dev.off()
  }
}

#formatC(numb, format = "e", digits = 2)

write.csv(bigplot2,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot2.june2020.txt",row.names=F)
