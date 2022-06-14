files=list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/tajimas/",pattern="CHRFIX.pestPG",full.names = T)

specieslist = c("bil","fla","bru","sin","fus","nit","mel","cri","cur","bel")

for(species in specieslist) {
  print(species)
  subsetfiles = files[grepl(x=files,pattern=species)]
  noweirds = grepl(x=subsetfiles,pattern="NOWEIRD")
  if(sum(noweirds)>0) {
    subsetfiles = subsetfiles[noweirds]
    species = paste(species,"-NOWEIRD",sep="")
  }
  # if(sum(!(noweirds))>0) {
  #   subsetfiles = subsetfiles[!(noweirds)]
  # }
  
  if(length(subsetfiles)>=3) { 
    
    chi = read.table(subsetfiles[grepl(x=subsetfiles,pattern="CHI")][1],header=T,sep="\t")
    son = read.table(subsetfiles[grepl(x=subsetfiles,pattern="SON")][1],header=T,sep="\t")
    both = read.table(subsetfiles[!(grepl(x=subsetfiles,pattern="SON")) & !(grepl(x=subsetfiles,pattern="CHI"))][1],header=T,sep="\t")
    
    chi = chi[,c(2:3,9)]
    son = son[,c(2:3,9)]
    both = both[,c(2:3,9)]
    
    colnames(chi)[3] = "chi_Tajima"
    colnames(son)[3] = "son_Tajima"
    colnames(both)[3] = "both_Tajima"
    
    merged = merge(chi,son,all=T)
    merged = merge(merged,both,all=T)
    # correlation = cor(merged[,3:5],use="pairwise.complete.obs")
    # write.table(correlation,file=paste(species,"_tajimas_correlation_son_chi_both.txt",sep=""))
    
    # png(paste(species,"_tajimas_correlation_son_chi_both.png",sep=""),height=1200)
    # par(mfrow=c(3,1))
    # plot(merged[,c(3:4)])
    # abline(a=0,b=1,col="red")
    # plot(merged[,c(3,5)])
    # abline(a=0,b=1,col="red")
    # plot(merged[,c(4:5)])
    # abline(a=0,b=1,col="red")
    # dev.off()
    
    # png(paste(species,"_tajimas_hist.png",sep=""),height=1200)
    # par(mfrow=c(3,1))
    # hist(merged[,3]-merged[,4])
    # abline(v=0,col="red")
    # hist(merged[,3]-merged[,5])
    # abline(v=0,col="red")
    # hist(merged[,4]-merged[,5])  
    # abline(v=0,col="red")
    # dev.off()
    
    png(paste(species,"_tajimas_genome.png",sep=""),height=9,width=12,units="in",res=600)
    par(mfrow=c(3,1))
    plot(merged[,3],ylim=c(min(merged[,3:5],na.rm=T),max(merged[,3:5],na.rm=T)))
    plot(merged[,4],ylim=c(min(merged[,3:5],na.rm=T),max(merged[,3:5],na.rm=T)))
    plot(merged[,5],ylim=c(min(merged[,3:5],na.rm=T),max(merged[,3:5],na.rm=T)))
    dev.off()
    
  }
}
