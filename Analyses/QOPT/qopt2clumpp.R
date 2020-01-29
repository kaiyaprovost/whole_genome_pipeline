## qopt to clumpp format 

##1 1 (x) 1 : 0.315 0.002 0.683
##2 2 (x) 1 : 0.315 0.002 0.683

path ="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/pca qopt files/chroms/" 
setwd(path)
#spp = "Auriparus-flaviceps"

filelist = list.files(pattern="_PCAngsd",recursive = T)

#qopt="PseudoNC_011465.1_Tgut_2.0.F.Auriparus-flaviceps_PCAngsd.4.K5.a0.qopt"

for (i in 1:length(filelist)) {
  
  file = filelist[i]
  
  print(file)
  
  split = strsplit(file,"/")[[1]]
  
  spp = split[1]
  K = split[2]
  qopt=split[3]
  
  outfile = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/CLUMPP_MacOSX.1.1.2/",spp,".",K,".qopt.all.clumpp.txt",sep="")
  listfile = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/CLUMPP_MacOSX.1.1.2/",spp,".",K,".chromosomes.order.txt",sep="")
  
  
  x = read.delim(file,header=F,sep=" ")
  
  k = ncol(x)
  inds = nrow(x)
  indcol = 1:inds
  
  add=cbind(indcol,indcol,rep("(0)",inds),rep("1 :",inds),x)
  
  if(i==1) {
    write.table(add,outfile,append=F,quote=F,sep=" ",col.names=F,row.names=F)
    cat(qopt,"\n",file=listfile,append=T)
    
  } else {
    write.table(add,outfile,append=T,quote=F,sep=" ",col.names=F,row.names=F)
    cat(qopt,"\n",file=listfile,append=T)
    
  }
  
  cat("\n",file=outfile,append=T)
  
}

