library(hzar)
library(RcppCNPy)

setwd("~")

## take in arguments

species="bellii"
plotWholeSpecies=T

## set up functions
preprocess_hzar= function(qopt_ada,chain=1e5,burn=5e2,low=0,high=20,scaling="fixed",tails="none",verbose=F){
  qopt_model <- hzar.makeCline1DFreq(data=qopt_ada, scaling=scaling,tails=tails);
  qopt_model <- hzar.model.addBoxReq(meta.model=qopt_model,low=low,high=high)
  qopt_model_FitR <- hzar.first.fitRequest.old.ML(model=qopt_model,obsData=qopt_ada,verbose=verbose);
  qopt_model_FitR$mcmcParam$chainLength <- chain;
  qopt_model_FitR$mcmcParam$burnin <- burn;
  
  return(qopt_model_FitR)
  
}

postprocess_hzar = function(qopt_ada,qopt_model_data,add=F,low=0,high=20,index=1) {
  
  qopt_model_null <- hzar.dataGroup.null(qopt_ada);
  qopt_null_ada <- list(clineModel = qopt_model_data,nullModel  = qopt_model_null);
  qopt_null_ada_2 <- hzar.make.obsDataGroup(qopt_null_ada);
  qopt_null_ada_2 <- hzar.copyModelLabels(qopt_null_ada,qopt_null_ada_2);
  
  aicc=(hzar.AICc.hzar.obsDataGroup(qopt_null_ada_2))
  delta = aicc - min(aicc,na.rm=T)
  cline_delta = delta[1,]
  null_delta = delta[2,]
  if (cline_delta > null_delta){
    ## print the null line
    hzar.plot.cline(qopt_model_null,xlim=c(low,high),ylim=c(0,1),col=index,add=add,lty=3)
  } else if (null_delta - cline_delta > 10){
    ## print the line bold
    hzar.plot.cline(qopt_model_data,xlim=c(low,high),ylim=c(0,1),col=index,add=add,lwd=2); ## i think this is the same as previous curve
  } else if (null_delta - cline_delta > 2){
    ## print the line regular
    hzar.plot.cline(qopt_model_data,xlim=c(low,high),ylim=c(0,1),col=index,add=add); ## i think this is the same as previous curve
  } else {
    ## print the line dashed 
    hzar.plot.cline(qopt_model_data,xlim=c(low,high),ylim=c(0,1),col=index,add=add,lty=3); ## i think this is the same as previous curve
    
  }
}

run_hzar=function(qopt_ada,plot=T,index=1,add=F){
  if(index==1 | add==F) { add=F } else { add=T}
  
  qopt_model_FitR=preprocess_hzar(qopt_ada)
  if(plot==T){hzar.plot.obsData(qopt_ada)}
  
  qopt_model_Fit <- hzar.doFit(qopt_model_FitR)
  if(plot==T){plot(hzar.mcmc.bindLL(qopt_model_Fit))}
  
  qopt_model_data <- hzar.dataGroup.add(qopt_model_Fit);
  qopt_model_data <- hzar.dataGroup.add(qopt_model_data,hzar.chain.doSeq(hzar.next.fitRequest(qopt_model_Fit)));
  hzar.plot.cline(qopt_model_data,xlim=c(0,20),ylim=c(0,1),col=index,add=add);
  if(plot==T){
    hzar.plot.fzCline(qopt_model_data); ## won't work with just raw longitude
    print(hzar.getLLCutParam(qopt_model_data,c("center","width")))
  };
  
  postprocess_hzar(qopt_ada,qopt_model_data,add,index=index)
  
  res<-list()
  res$mean_width<-mean(qopt_model_Fit$mcmcRaw[,2])
  # res$mean_width<-mean(qopt_model_Fit$mcmcRaw[,2])
  res$out<-qopt_model_Fit 
  res$dat<-qopt_model_data
  return(res)
}

no_mcmc_hzar = function(qopt_ada,index=1,add=F,plot=F){
  qopt_model_FitR=preprocess_hzar(qopt_ada)
  if(index==1 | add==F) { add=F } else { add=T}
  
  qopt_model_data <- hzar.dataGroup.add(qopt_model_FitR);
  qopt_model_data <- hzar.dataGroup.add(qopt_model_data,hzar.chain.doSeq(hzar.next.fitRequest(qopt_model_FitR)));
  
  postprocess_hzar(qopt_ada,qopt_model_data,add)
  
}

make_qopt_ada = function(filename,spp,desert="",weird="",removeblank=F) {
  
  if(tools::file_ext(filename)=="npy") {
    individual_qopt=npyLoad(x)
  } else {
    individual_qopt=read.table(filename,header=F)
  }
  
  
  colnames(individual_qopt)=paste("CLUSTER",seq(1:ncol(individual_qopt)),sep="")
  all_qopt = all_qopt[order(all_qopt$RG),]
  
  if(removeblank==T) {
    add_qopt = all_qopt[!(is.na(all_qopt$K2_A)),]
  } else {
    add_qopt = all_qopt
  }
  
  if(desert != ""){
    add_qopt = add_qopt[add_qopt$DES==desert,]
  } 
 
  
  #add_qopt=add_qopt[add_qopt$SP==paste(spp,weird,sep=""),] 
  spp = paste(spp,weird,sep="")
  add_qopt=add_qopt[add_qopt$SP==spp,] 
  
  add_qopt = cbind(add_qopt[add_qopt$SP==spp,c("RG","SP","LAT","LONG","distance","nSamples")],
                   individual_qopt)
  
  qopt_ada = hzar.doMolecularData1DPops(add_qopt$distance,
                                        add_qopt$CLUSTER1,
                                        add_qopt$nSamples)
  
  is_positive=(cor(qopt_ada$frame[,1:2],use="pairwise.complete.obs")[2])>0
  if(is_positive==F){
    qopt_ada = hzar.doMolecularData1DPops(add_qopt$distance,
                                          add_qopt$CLUSTER2,
                                          add_qopt$nSamples)
  }
  return(qopt_ada)
}

npyOutputFunction = function(i,smallsppfiles,smallspp=NULL,all_qopt){
  
  x=smallsppfiles[i]
  print(x)
  
  if(is.null(smallspp)){
    smallspp=strsplit(basename(x),"\\.")[[1]][1]
  }
  
  #spp=paste(strsplit(strsplit(strsplit(basename(x),"\\.")[[1]][5],"_")[[1]][1],"-")[[1]][-1],sep="-")
  
  splits=strsplit(basename(x),"\\.")[[1]]
  k = as.numeric(splits[length(splits)-3])+1
  #k=as.numeric(splits[3])+1
  spp=paste(strsplit(smallspp,"-")[[1]][1:2],sep="-",collapse="-")
  
  if(grepl("CHI",smallspp,fixed=T)){
    desert="CHI"
  } else if(grepl("SON",smallspp,fixed=T)){
    desert="SON"
  } else{
    desert=""
  }
  
  if (grepl("WEIRD",smallspp,fixed=T)){weird="-NOWEIRD"} else {weird=""}
  
  qopt_ada=make_qopt_ada(filename=x,spp=species,desert=desert,weird=weird)
  
  res=run_hzar(qopt_ada,plot=F,index=i)
  
  output=list()
  mean_width=res$mean_width
  center = mean(res$out$mcmcRaw[,1],na.rm=T)
  output$mean_width = mean_width
  output$center = center
  #output=c(mean_width,center)
  output2=as.data.frame(do.call(cbind, output))
  rownames(output2)=basename(x)
  textoutput=paste(smallspp,"_hzar_chroms.txt",sep="")
  write.table(output2,file=textoutput,append=T,quote=F,sep="\t")
  return(output)
  
  
  
  
  
  
}

## load lat long etc data
#all_qopt=read.table("AllSpeciesMetadata_allK.csv",
#                    sep=",",header=T)
## /Users/kprovost/Downloads/AllSpeciesMetadata_allK_9june2020.csv
all_qopt=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/QOPT/AllSpeciesMetadata_allK_5january2021.csv",
                    sep="\t",header=T)

all_qopt = all_qopt[order(all_qopt$RG),]
all_qopt$distance = all_qopt$LONG-min(all_qopt$LONG,na.rm=T) ## this may be degrees?
all_qopt$nSamples = 1

## do for full species values 
out = data.frame()
for(spp in unique(all_qopt$SP)) {

for(i in 1:3) {
  
  if(i==1){colnames=c("P2_A","P2_B")}else if(i==2){colnames=c("K2_A","K2_B")} else{colnames=c("C2_A","C2_B")}
  
  j=c("P","K","C")[i]
  
  add_qopt = all_qopt[all_qopt$SP==spp,c("RG","SP","LAT","LONG","distance","nSamples",colnames)]
  
  colnames(add_qopt) = c("RG","SP","LAT","LONG","distance","nSamples","CLUSTER1","CLUSTER2")
  
  add_qopt = add_qopt[complete.cases(add_qopt),]
  
  if(nrow(add_qopt)>0) {
    qopt_ada = hzar.doMolecularData1DPops(add_qopt$distance,
                                          add_qopt$CLUSTER2,
                                          add_qopt$nSamples)
    
    print(paste(nrow(add_qopt),spp,j))
    #res=run_hzar(qopt_ada,plot=F,index=1,add=T)
    #mean_width=res$mean_width
    #center = mean(res$out$mcmcRaw[,1],na.rm=T)
    #row = cbind(spp,j,mean_width,center)
    #out = rbind(out,row)
  }
  
}

}

out

## iterate through species files
#files = list.files("ADMIX_NUMPY_AUG_2020/",
#                   pattern="1.admix.Q.npy",full.names = T) ## or K2.a0.qopt
#files = list.files("/Users/kprovost/Downloads/K2/",
#                   pattern="K2.a0.qopt",full.names = T) ## or K2.a0.qopt

#files=list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES",
#                 pattern="1.admix.Q.npy",full.names = T)

files=list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES",
                 pattern="1.admix.Q.npy.csv",full.names = T,recursive=T)

sppfiles=(files[grepl(species,files)])

for(species in (c("flaviceps-NOWEIRD"))){
  sppfiles=(files[grepl(species,files)])
  #sppfiles=sppfiles[1:5]
  pdf(paste(species,"_test_nomcmc.hzar.pdf",sep=""))
  for(i in 1:length(sppfiles)){
    print(i)
    chrom=strsplit(strsplit(basename(sppfiles[i]),"\\.")[[1]][3], "_")[[1]][3]
    qopt_ada=make_qopt_ada(filename=sppfiles[i],spp=species)
    no_mcmc_hzar(qopt_ada,i)
    title(chrom)
  }
  dev.off()
}; dev.off()

for(species in (c("bellii-NOWEIRD"))){
  spp=species
  sppfiles=(files[grepl(species,files)])
  #sppfiles=sppfiles[1:5]
  textoutput=paste(species,"_hzar_chroms_npy.txt",sep="")
  
  pdf(paste(species,"_hzar_chroms_add.pdf",sep=""))
  output=lapply(1:length(sppfiles),FUN = function(i){
    x=sppfiles[i]
    print(x)
    spp=paste(strsplit(strsplit(strsplit(basename(x),"\\.")[[1]][5],"_")[[1]][1],"-")[[1]][-1],sep="-")
    #individual_qopt = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/pca qopt files/chroms/BEL/K2/PseudoNC_007897.1_Tgut_mtDNA.0.F.Vireo-bellii_PCAngsd.1.K2.a0.qopt",
    #                             header=F)
    qopt_ada=make_qopt_ada(filename=x,spp=species,removeblank = T)
    res=run_hzar(qopt_ada,plot=F,index=i,add=T)
    output=list()
    mean_width=res$mean_width
    center = mean(res$out$mcmcRaw[,1],na.rm=T)
    output$mean_width = mean_width
    output$center = center
    #output=c(mean_width,center)
    output2=as.data.frame(do.call(cbind, output))
    rownames(output2)=basename(x)
    
    write.table(output2,file=textoutput,append=T,quote=F,sep="\t")
    return(output)
  }); 
  print(output)
  dev.off()
  
  pdf(paste(species,"_hzar_width_center.pdf",sep=""))
  chrom=sapply(sppfiles,FUN=function(x){strsplit(strsplit(basename(x),"\\.")[[1]][2], "_")[[1]][3]})
  plot(0,xlim=c(0,20),ylim=c(1,length(output)+1),type="n",yaxt="n",ylab="",xlab="Longitude")
  title(spp)
  axis(2,1:(length(output)+1),c(chrom,"MEAN"),las=2,cex.axis=0.5)
  for(j in 1:length(output)){
    center=output[[j]]$center
    width=output[[j]]$mean_width
    start=center+(width/2)
    end=center-(width/2)
    lines(c(start,center,end),c(j,j,j))
    points(center,j)
  }
  output2=as.data.frame(do.call(rbind, output))
  width_m=mean(as.numeric(output2$mean_width))
  center_m=mean(as.numeric(output2$center))
  width_s =sd(as.numeric(output2$mean_width))
  center_s=sd(as.numeric(output2$center))
  
  start_m=center_m+(width_m/2)
  end_m=center_m-(width_m/2)
  lines(c(start_m,center_m,end_m),rep(length(output)+1,3),col="red")
  points(center_m,length(output)+1,col="red",pch=0)
  #abline(v=start_m,col="red",lty=3)
  #abline(v=end_m,col="red",lty=3)
  #abline(v=center_m,col="red",lty=2)
  
  abline(v=center_m+center_s,col="red",lty=3)
  abline(v=center_m-center_s,col="red",lty=3)
  
  abline(v=start_m+(width_s/2),col="red",lty=3)
  abline(v=start_m-(width_s/2),col="red",lty=3)
  
  abline(v=end_m-(width_s/2),col="red",lty=3)
  abline(v=end_m+(width_s/2),col="red",lty=3)
  
  
  par(mfrow=c(2,1))
  hist(as.numeric(output2$center),breaks=20,main=spp,xlab="Center Location")
  abline(v=center_m,col="red")
  abline(v=center_m+center_s,col="red",lty=3)
  abline(v=center_m-center_s,col="red",lty=3)
  
  hist(as.numeric(output2$mean_width),breaks=20,main=spp,xlab="Width")
  abline(v=width_m,col="red")
  abline(v=width_m+width_s,col="red",lty=3)
  abline(v=width_m-width_s,col="red",lty=3)
  
  par(mfrow=c(1,1))
  plot(as.numeric(output2$mean_width),as.numeric(output2$center),main=spp,xlab="Width",ylab="Center")
  
  dev.off()
  
}; dev.off()



final_data=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/QOPT/HZAR/Hzar_cline_overlap_data/data_tables/ALLSPECIES_hzar_chroms.txt",
                      sep="\t",header=T)

final_data$chromosome=factor(final_data$chromosome,levels=c(1,"1A","1B",2:4,"4A",5:28,"LG2","LG5","LGE22","mtDNA","Z"))
final_data$species=factor(final_data$species,levels=c("Vireo-bellii-NOWEIRD","Amphispiza-bilineata","Campylorhynchus-brunneicapillus",
                                                      "Toxostoma-crissale","Toxostoma-curvirostre","Auriparus-flaviceps-NOWEIRD",
                                                      "Melozone-fusca","Polioptila-melanura","Phainopepla-nitens","Cardinalis-sinuatus",
                                                      "Vireo-bellii","Auriparus-flaviceps"))

pdf("width_of_chromosome_clines_all_species.pdf",width=10)
boxplot(final_data$mean_width~final_data$chromosome,las=2)
boxplot(final_data$center~final_data$chromosome,las=2)

boxplot(final_data$mean_width~final_data$species,las=2,names=c("bel","bil","bru","cri","cur","fla","fus","mel","nit","sin","bel-W","fla-W"))
boxplot(final_data$center~final_data$species,las=2,names=c("bel","bil","bru","cri","cur","fla","fus","mel","nit","sin","bel-W","fla-W"))
dev.off()


if(length(sppfiles)>=1){
  
  allsplits=unlist(lapply(sppfiles,FUN=function(x){strsplit(basename(x),"\\.")[[1]][1]}))
  
  if(plotWholeSpecies==T){pdf(paste(species,"_hzar_npy_chroms.pdf",sep=""))
    
    output=lapply(1:length(sppfiles),FUN = function(i){npyOutputFunction(i,smallsppfiles=sppfiles,smallspp=NULL,all_qopt)})
    dev.off()
    
    pdf(paste(smallspp,"_hzar_width_center.pdf",sep=""))
    chrom=sapply(sppfiles,FUN=function(x){strsplit(strsplit(basename(x),"\\.")[[1]][2], "_")[[1]][3]})
    plot(0,xlim=c(0,20),ylim=c(1,length(output)+1),type="n",yaxt="n",ylab="",xlab="Longitude")
    title(spp)
    axis(2,1:(length(output)+1),c(chrom,"MEAN"),las=2,cex.axis=0.5)
    for(j in 1:length(output)){
      center=output[[j]]$center
      width=output[[j]]$mean_width
      start=center+(width/2)
      end=center-(width/2)
      lines(c(start,center,end),c(j,j,j))
      points(center,j)
    }
    output2=as.data.frame(do.call(rbind, output))
    width_m=mean(as.numeric(output2$mean_width))
    center_m=mean(as.numeric(output2$center))
    width_s =sd(as.numeric(output2$mean_width))
    center_s=sd(as.numeric(output2$center))
    
    start_m=center_m+(width_m/2)
    end_m=center_m-(width_m/2)
    lines(c(start_m,center_m,end_m),rep(length(output)+1,3),col="red")
    points(center_m,length(output)+1,col="red",pch=0)
    #abline(v=start_m,col="red",lty=3)
    #abline(v=end_m,col="red",lty=3)
    #abline(v=center_m,col="red",lty=2)
    
    abline(v=center_m+center_s,col="red",lty=3)
    abline(v=center_m-center_s,col="red",lty=3)
    
    abline(v=start_m+(width_s/2),col="red",lty=3)
    abline(v=start_m-(width_s/2),col="red",lty=3)
    
    abline(v=end_m-(width_s/2),col="red",lty=3)
    abline(v=end_m+(width_s/2),col="red",lty=3)
    
    
    par(mfrow=c(2,1))
    hist(as.numeric(output2$center),breaks=20,main=spp,xlab="Center Location")
    abline(v=center_m,col="red")
    abline(v=center_m+center_s,col="red",lty=3)
    abline(v=center_m-center_s,col="red",lty=3)
    
    hist(as.numeric(output2$mean_width),breaks=20,main=spp,xlab="Width")
    abline(v=width_m,col="red")
    abline(v=width_m+width_s,col="red",lty=3)
    abline(v=width_m-width_s,col="red",lty=3)
    
    plot(as.numeric(output2$mean_width),as.numeric(output2$center),main=spp,xlab="Width",ylab="Center")
    
    dev.off()
    
    
    
  } else if(plotWholeSpecies==F){
    
    for(smallspp in sample(unique(allsplits))){
      print(smallspp)
      smallsppfiles=sppfiles[allsplits==smallspp]
      
      pdf(paste(species,"_",smallspp,"_hzar_npy_chroms.pdf",sep=""))
      
      
      if(length(smallsppfiles)>1){
        output=lapply(1:length(smallsppfiles),FUN = function(i){npyOutputFunction(i,smallsppfiles,smallspp,all_qopt)})
      } else {
        output=npyOutputFunction(i,smallsppfiles,smallspp,all_qopt)
      }
      dev.off()
      
      
      pdf(paste(smallspp,"_hzar_width_center.pdf",sep=""))
      chrom=sapply(sppfiles,FUN=function(x){strsplit(strsplit(basename(x),"\\.")[[1]][2], "_")[[1]][3]})
      plot(0,xlim=c(0,20),ylim=c(1,length(output)+1),type="n",yaxt="n",ylab="",xlab="Longitude")
      title(spp)
      axis(2,1:(length(output)+1),c(chrom,"MEAN"),las=2,cex.axis=0.5)
      for(j in 1:length(output)){
        center=output[[j]]$center
        width=output[[j]]$mean_width
        start=center+(width/2)
        end=center-(width/2)
        lines(c(start,center,end),c(j,j,j))
        points(center,j)
      }
      output2=as.data.frame(do.call(rbind, output))
      width_m=mean(as.numeric(output2$mean_width))
      center_m=mean(as.numeric(output2$center))
      width_s =sd(as.numeric(output2$mean_width))
      center_s=sd(as.numeric(output2$center))
      
      start_m=center_m+(width_m/2)
      end_m=center_m-(width_m/2)
      lines(c(start_m,center_m,end_m),rep(length(output)+1,3),col="red")
      points(center_m,length(output)+1,col="red",pch=0)
      #abline(v=start_m,col="red",lty=3)
      #abline(v=end_m,col="red",lty=3)
      #abline(v=center_m,col="red",lty=2)
      
      abline(v=center_m+center_s,col="red",lty=3)
      abline(v=center_m-center_s,col="red",lty=3)
      
      abline(v=start_m+(width_s/2),col="red",lty=3)
      abline(v=start_m-(width_s/2),col="red",lty=3)
      
      abline(v=end_m-(width_s/2),col="red",lty=3)
      abline(v=end_m+(width_s/2),col="red",lty=3)
      
      
      par(mfrow=c(2,1))
      hist(as.numeric(output2$center),breaks=20,main=spp,xlab="Center Location")
      abline(v=center_m,col="red")
      abline(v=center_m+center_s,col="red",lty=3)
      abline(v=center_m-center_s,col="red",lty=3)
      
      hist(as.numeric(output2$mean_width),breaks=20,main=spp,xlab="Width")
      abline(v=width_m,col="red")
      abline(v=width_m+width_s,col="red",lty=3)
      abline(v=width_m-width_s,col="red",lty=3)
      
      plot(as.numeric(output2$mean_width),as.numeric(output2$center),main=spp,xlab="Width",ylab="Center")
      
      dev.off()
      
      
      
    }
  }
  
} 
dev.off()





##### 
## OLD
# run_hzar=function(qopt_ada,plot=T,index=1){
#   
#   
#   if(plot==T){hzar.plot.obsData(qopt_ada)}
#   
#   qopt_model <-
#     hzar.makeCline1DFreq(data=qopt_ada, scaling="fixed",tails="none");
#   
#   qopt_model <-
#     hzar.model.addBoxReq(meta.model=qopt_model,
#                          low=0,high=20);
#   
#   qopt_model_FitR <-
#     hzar.first.fitRequest.old.ML(model=qopt_model ,
#                                  obsData=qopt_ada,
#                                  verbose=FALSE);
#   
#   qopt_model_FitR$mcmcParam$chainLength <- 1e5;
#   qopt_model_FitR$mcmcParam$burnin <- 5e2;
#   qopt_model_Fit <- hzar.doFit(qopt_model_FitR)
#   if(plot==T){plot(hzar.mcmc.bindLL(qopt_model_Fit))}
#   qopt_model_data <-
#     hzar.dataGroup.add(qopt_model_Fit);
#   
#   qopt_model_data <-
#     hzar.dataGroup.add(
#       qopt_model_data,
#       hzar.chain.doSeq(hzar.next.fitRequest(qopt_model_Fit)));
#   hzar.plot.cline(qopt_model_data,xlim=c(0,20),ylim=c(0,1),col=index,
#                   add=index!=1);
#   if(plot==T){hzar.plot.fzCline(qopt_model_data); ## won't work with just raw longitude
#   }
#   
#   if(plot==T){print(hzar.getLLCutParam(qopt_model_data,c("center","width")))};
#   qopt_model_null <- hzar.dataGroup.null(qopt_ada);
#   qopt_null_ada <- list(clineModel = qopt_model_data,
#                         nullModel  = qopt_model_null);
#   qopt_null_ada_2 <- hzar.make.obsDataGroup(qopt_null_ada);
#   qopt_null_ada_2 <- hzar.copyModelLabels(qopt_null_ada,qopt_null_ada_2);
#   if(plot==T){hzar.plot.cline(qopt_null_ada_2,xlim=c(0,20),ylim=c(0,1),col=index,
#                               add=index!=1)
#     print(hzar.AICc.hzar.obsDataGroup(qopt_null_ada_2))
#   }
#   
#   res<-list()
#   res$mean_width<-mean(qopt_model_Fit$mcmcRaw[,2])
#   res$out<-qopt_model_Fit 
#   res$dat<-qopt_model_data
#   return(res)
#   
# }

