library(hzar)
setwd("~")

## EDIT RUN_HZAR TO MAKE IT PRINT ACCORDING TO AICC
## JUST NEED BARE ESTIMATE FROM DATA IN FAST_HZAR

#all_qopt=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_allK.csv",
#                    sep=",",header=T)
all_qopt=read.table("/Users/kprovost/Downloads/AllSpeciesMetadata_allK_9june2020.csv",
                    sep="\t",header=T)
all_qopt = all_qopt[order(all_qopt$RG),]

all_qopt$distance = all_qopt$LONG-min(all_qopt$LONG,na.rm=T)
all_qopt$nSamples = 1

preprocess_hzar= function(qopt_ada,chain=1e5,burn=5e2,low=0,high=20,scaling="fixed",tails="none",verbose=F){
  qopt_model <- hzar.makeCline1DFreq(data=qopt_ada, scaling=scaling,tails=tails);
  qopt_model <- hzar.model.addBoxReq(meta.model=qopt_model,low=low,high=high)
  qopt_model_FitR <- hzar.first.fitRequest.old.ML(model=qopt_model,obsData=qopt_ada,verbose=verbose);
  qopt_model_FitR$mcmcParam$chainLength <- chain;
  qopt_model_FitR$mcmcParam$burnin <- burn;
  
  return(qopt_model_FitR)
  
}

postprocess_hzar = function(qopt_ada,qopt_model_data,add=F,low=0,high=20) {
  
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
  
  postprocess_hzar(qopt_ada,qopt_model_data,add)
  
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

make_qopt_ada = function(filename,spp) {
  individual_qopt=read.table(filename,header=F)
  colnames(individual_qopt)=paste("CLUSTER",seq(1:ncol(individual_qopt)),sep="")
  all_qopt = all_qopt[order(all_qopt$RG),]
  add_qopt = cbind(all_qopt[all_qopt$SP==spp,c("RG","SP","LAT","LONG","distance","nSamples")],
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

#files=list.files("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/pca qopt files/chroms/",
#                 full.names = T,recursive=T,pattern="qopt$")

files=list.files("/Users/kprovost/Downloads/K2/",
                 full.names = T,recursive=T,pattern="qopt$")
## whichever one is 32 is the one that is weird 

files=files[grepl("0.F",files)]
files=files[grepl("K2",files)]

for(species in (c(
  "bellii",
  "bilineata","brunneicapillus",
  "crissale",
  "curvirostre","flaviceps","fusca",
  "melanura","nitens",
  "sinuatus"
))){
  
  sppfiles=(files[grepl(species,files)])
  #sppfiles=sppfiles[1:5]
  pdf(paste(species,"_test.hzar.pdf",sep=""))
  for(i in 1:length(sppfiles)){
    print(i)
    chrom=strsplit(strsplit(basename(sppfiles[i]),"\\.")[[1]][2], "_")[[1]][3]
    qopt_ada=make_qopt_ada(filename=sppfiles[i],spp=species)
    no_mcmc_hzar(qopt_ada,i)
    title(chrom)
  }
  dev.off()
}; dev.off()

for(species in (c(
  #"bellii",
  #"bellii-NOWEIRD","flaviceps-NOWEIRD",
  "bilineata","brunneicapillus",
  "crissale",
  "curvirostre","flaviceps","fusca",
  "melanura","nitens",
  "sinuatus"
))){
  spp=species
  sppfiles=(files[grepl(species,files)])
  #sppfiles=sppfiles[1:5]
  textoutput=paste(species,"_hzar_chroms.txt",sep="")
  
  pdf(paste(species,"_hzar_chroms.pdf",sep=""))
  output=lapply(1:length(sppfiles),FUN = function(i){
    x=sppfiles[i]
    print(x)
    spp=paste(strsplit(strsplit(strsplit(basename(x),"\\.")[[1]][5],"_")[[1]][1],"-")[[1]][-1],sep="-")
    #individual_qopt = read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/pca qopt files/chroms/BEL/K2/PseudoNC_007897.1_Tgut_mtDNA.0.F.Vireo-bellii_PCAngsd.1.K2.a0.qopt",
    #                             header=F)
    qopt_ada=make_qopt_ada(filename=x,spp=species)
    res=run_hzar(qopt_ada,plot=F,index=i)
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
  
  plot(as.numeric(output2$mean_width),as.numeric(output2$center),main=spp,xlab="Width",ylab="Center")
  
  dev.off()
  
}; dev.off()


npyOutputFunction = function(i,smallsppfiles,smallspp=NULL,all_qopt,fast=T){
  
  x=smallsppfiles[i]
  #print(x)
  print(paste(i,length(smallsppfiles),sep="/"))
  
  if(is.null(smallspp)){
    smallspp=strsplit(basename(x),"\\.")[[1]][1]
  }
  
  #spp=paste(strsplit(strsplit(strsplit(basename(x),"\\.")[[1]][5],"_")[[1]][1],"-")[[1]][-1],sep="-")
  
  splits=strsplit(basename(x),"\\.")[[1]]
  k=as.numeric(splits[length(splits)-3])+1
  #k=as.numeric(splits[3])+1
  spp=paste(strsplit(smallspp,"-")[[1]][1:2],sep="-",collapse="-")
  
  if(grepl("CHI",smallspp,fixed=T)){desert="CHI"} else if(grepl("SON",smallspp,fixed=T)){desert="SON"} else {desert=""}
  
  if (grepl("WEIRD",smallspp,fixed=T)){weird="-NOWEIRD"} else {weird=""}
  
  individual_qopt=npyLoad(x)
  
  colnames(individual_qopt)=paste("CLUSTER",seq(1:ncol(individual_qopt)),sep="")
  all_qopt = all_qopt[order(all_qopt$RG),]
  add_qopt = all_qopt[!(is.na(all_qopt$K2_A)),]
  #add_qopt = all_qopt[,c("RG","SP","LAT","LONG","distance","nSamples","DES")]
  
  if(desert != ""){add_qopt = add_qopt[add_qopt$DES==desert,]} 
  
  add_qopt=add_qopt[add_qopt$SP==paste(species,weird,sep=""),]
  
  add_qopt = cbind(add_qopt,individual_qopt)
  
  qopt_ada = hzar.doMolecularData1DPops(add_qopt$distance,add_qopt$CLUSTER1,add_qopt$nSamples)
  
  
  is_positive=(cor(qopt_ada$frame[,1:2],use="pairwise.complete.obs")[2])>0
  if(is_positive==F){qopt_ada = hzar.doMolecularData1DPops(add_qopt$distance,add_qopt$CLUSTER2,add_qopt$nSamples)}
  
  if(fast==F){
    res=run_hzar(qopt_ada,plot=F,index=i)
    output=list()
    mean_width=res$mean_width
    center = mean(res$out$mcmcRaw[,1],na.rm=T)
    output$mean_width = mean_width
    output$center = center
    #output=c(mean_width,center)
    return(output)
  } else {
    no_mcmc_hzar(qopt_ada,index=i)
    return(NULL)
  }
  
}
library(RcppCNPy)
files = list.files("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/ADMIX_NUMPY_AUG_2020/",
                   pattern="1.admix.Q.npy",full.names = T,recursive = T) ## PCA1 = k2
all_qopt=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_allK.csv",
                    sep=",",header=T)
all_qopt = all_qopt[order(all_qopt$RG),]
all_qopt$distance = all_qopt$LONG-min(all_qopt$LONG,na.rm=T)
all_qopt$nSamples = 1
setwd("~")
for(fast in c(F)){
  for(plotWholeSpecies in c(F)){
    
    for(species in sort(c(
      "bellii"#,
      #"bilineata","brunneicapillus",
      #"crissale","curvirostre","flaviceps","fusca",
      #"melanura","nitens","sinuatus"
    ))){
      
      sppfiles=(files[grepl(species,files)])
      
      if(length(sppfiles)>=1){
        
        allsplits=unlist(lapply(sppfiles,FUN=function(x){strsplit(basename(x),"\\.")[[1]][1]}))
        ## need to fix this so that if its got pseudonc in it its pulling from tgut instead
        
        if(plotWholeSpecies==T){pdf(paste(species,"_fast",as.character(fast),"_hzar_npy_chroms.pdf",sep=""))
          
          output=lapply(1:length(sppfiles),FUN = function(i){npyOutputFunction(i,smallsppfiles=sppfiles,smallspp=NULL,all_qopt,fast=fast)})
          dev.off()
          
        } else if(plotWholeSpecies==F){
          
          for(smallspp in sample(unique(allsplits))){
            print(smallspp)
            smallsppfiles=sppfiles[allsplits==smallspp]
            
            pdf(paste(species,"_",smallspp,"_fast",as.character(fast),"_hzar_npy_chroms.pdf",sep=""))
            
            
            if(length(smallsppfiles)>1){
              output=lapply(1:length(smallsppfiles),FUN = function(i){npyOutputFunction(i,smallsppfiles,smallspp,all_qopt,fast=fast)})
            } else {
              output=npyOutputFunction(i,smallsppfiles,smallspp,all_qopt,fast=fast)
            }
            dev.off()
          }
        }
      }
      
    } 
  }
}
dev.off()


#output = (do.call(rbind, output))

#beeswarm::beeswarm(unlist(output[,1]),cex=0.6,col="blue")
#beeswarm::beeswarm(unlist(output[,2]),cex=0.6,col="blue")

#plot(unlist(output[,1]),unlist(output[,2]),
#     xlab="Location",ylab="Width")






#####
for(spp in (c("bellii"#,"bellii-NOWEIRD","bilineata",
              #"brunneicapillus","crissale","curvirostre",
              #"flaviceps","flaviceps-NOWEIRD"#,
              #  "fusca","melanura","nitens","sinuatus"
))) {
  print(spp)
  
  for(method in rev(c("K","P","C"))) {
    print(method)
    
    
    
    pdf(paste(spp,"_",method,"_hzar.pdf",sep=""))
    
    qopt=all_qopt[all_qopt$SP==spp,]
    
    if(method=="K"){
      qopt=qopt[,c("SP","distance","K2_A","nSamples")]
      qopt=qopt[complete.cases(qopt),]
      qopt_ada = hzar.doMolecularData1DPops(qopt$distance[qopt$SP==spp],
                                            qopt$K2_A[qopt$SP==spp],
                                            qopt$nSamples[qopt$SP==spp])
    } else if(method=="P"){
      qopt=qopt[,c("SP","distance","P2_A","nSamples")]
      qopt=qopt[complete.cases(qopt),]
      qopt_ada = hzar.doMolecularData1DPops(qopt$distance[qopt$SP==spp],
                                            qopt$P2_A[qopt$SP==spp],
                                            qopt$nSamples[qopt$SP==spp])
    } else if(method=="C"){
      qopt=qopt[,c("SP","distance","C2_A","nSamples")]
      qopt=qopt[complete.cases(qopt),]
      qopt_ada = hzar.doMolecularData1DPops(qopt$distance[qopt$SP==spp],
                                            qopt$C2_A[qopt$SP==spp],
                                            qopt$nSamples[qopt$SP==spp])
    } else {
      next
    }
    
    res=run_hzar(qopt_ada)
    
    
    dev.off()
    
  }
  
}


#####

## old


# run_hzar=function(qopt_ada,plot=T,index=1,add=F){
#   if(plot==T){hzar.plot.obsData(qopt_ada)}
#   qopt_model <-
#     hzar.makeCline1DFreq(data=qopt_ada, scaling="fixed",tails="none");
#   qopt_model <-
#     hzar.model.addBoxReq(meta.model=qopt_model,
#                          low=0,high=20)
#   qopt_model_FitR <-
#     hzar.first.fitRequest.old.ML(model=qopt_model ,
#                                  obsData=qopt_ada,
#                                  verbose=FALSE);
#   qopt_model_FitR$mcmcParam$chainLength <- 1e5;
#   qopt_model_FitR$mcmcParam$burnin <- 5e2;
#   qopt_model_Fit <- hzar.doFit(qopt_model_FitR)
#   if(plot==T){plot(hzar.mcmc.bindLL(qopt_model_Fit))}
#   qopt_model_data <- hzar.dataGroup.add(qopt_model_Fit);
#   qopt_model_data <-
#     hzar.dataGroup.add(qopt_model_data,
#                        hzar.chain.doSeq(hzar.next.fitRequest(qopt_model_Fit)));
#   hzar.plot.cline(qopt_model_data,xlim=c(0,20),ylim=c(0,1),col=index,
#                   add=index!=1);
#   if(plot==T){hzar.plot.fzCline(qopt_model_data); ## won't work with just raw longitude
#   }
#   if(plot==T){print(hzar.getLLCutParam(qopt_model_data,c("center","width")))};
#   qopt_model_null <- hzar.dataGroup.null(qopt_ada);
#   qopt_null_ada <- list(clineModel = qopt_model_data,
#                         nullModel  = qopt_model_null);
#   qopt_null_ada_2 <- hzar.make.obsDataGroup(qopt_null_ada);
#   qopt_null_ada_2 <- hzar.copyModelLabels(qopt_null_ada,qopt_null_ada_2);
#   aicc=(hzar.AICc.hzar.obsDataGroup(qopt_null_ada_2))
#   delta = aicc - min(aicc,na.rm=T)
#   cline_delta = delta[1,]
#   null_delta = delta[2,]
#   if (cline_delta > null_delta){
#     ## print the null line
#     hzar.plot.cline(qopt_model_null,xlim=c(0,20),ylim=c(0,1),col=index,
#                     add=index!=1,lty=3)
#   } else if (null_delta - cline_delta > 10){
#     ## print the line bold
#     hzar.plot.cline(qopt_model_data,xlim=c(0,20),ylim=c(0,1),col=index,
#                     add=index!=1,lwd=2); ## i think this is the same as previous curve
#   } else if (null_delta - cline_delta > 2){
#     ## print the line regular
#     hzar.plot.cline(qopt_model_data,xlim=c(0,20),ylim=c(0,1),col=index,
#                     add=index!=1); ## i think this is the same as previous curve
#   } else {
#     ## print the line dashed 
#     hzar.plot.cline(qopt_model_data,xlim=c(0,20),ylim=c(0,1),col=index,
#                     add=index!=1,lty=3); ## i think this is the same as previous curve
#     
#   }
#   res<-list()
#   res$mean_width<-mean(qopt_model_Fit$mcmcRaw[,2])
#   # res$mean_width<-mean(qopt_model_Fit$mcmcRaw[,2])
#   res$out<-qopt_model_Fit 
#   res$dat<-qopt_model_data
#   return(res)
# }
# 
# no_mcmc_hzar = function(qopt_ada,index=1,add=F){
#   if(index==1 | add==F) { add=F } else { add=T}
#   qopt_model <-    hzar.makeCline1DFreq(data=qopt_ada, scaling="fixed",tails="none");
#   qopt_model <-    hzar.model.addBoxReq(meta.model=qopt_model,
#                                         low=0,high=20);
#   qopt_model_FitR <-
#     hzar.first.fitRequest.old.ML(model=qopt_model ,
#                                  obsData=qopt_ada,
#                                  verbose=FALSE);
#   qopt_model_FitR$mcmcParam$chainLength <- 1e5;
#   qopt_model_FitR$mcmcParam$burnin <- 5e2;
#   qopt_model_data2 <- hzar.dataGroup.add(qopt_model_FitR);
#   qopt_model_data2 <-
#     hzar.dataGroup.add(
#       qopt_model_data2,
#       hzar.chain.doSeq(hzar.next.fitRequest(qopt_model_FitR)));
#   qopt_model_null <- hzar.dataGroup.null(qopt_ada);
#   qopt_null_ada_A <- list(clineModel = qopt_model_data2,
#                           nullModel  = qopt_model_null);
#   qopt_null_ada_A_2 <- hzar.make.obsDataGroup(qopt_null_ada_A);
#   qopt_null_ada_A_2 <- hzar.copyModelLabels(qopt_null_ada_A,qopt_null_ada_A_2);
#   aicc=(hzar.AICc.hzar.obsDataGroup(qopt_null_ada_A_2))
#   delta = aicc - min(aicc,na.rm=T)
#   cline_delta = delta[1,]
#   null_delta = delta[2,]
#   if (cline_delta > null_delta){
#     ## print the null line
#     hzar.plot.cline(qopt_model_null,xlim=c(0,20),ylim=c(0,1),col=index,
#                     add=add)
#   } else if (null_delta - cline_delta > 10){
#     ## print the line bold
#     hzar.plot.cline(qopt_model_data2,xlim=c(0,20),ylim=c(0,1),col=index,
#                     add=add,lwd=2); ## i think this is the same as previous curve
#   } else if (null_delta - cline_delta > 2){
#     ## print the line regular
#     hzar.plot.cline(qopt_model_data2,xlim=c(0,20),ylim=c(0,1),col=index,
#                     add=add); ## i think this is the same as previous curve
#   } else {
#     ## print the line dashed 
#     hzar.plot.cline(qopt_model_data2,xlim=c(0,20),ylim=c(0,1),col=index,
#                     add=add,lty=3); ## i think this is the same as previous curve
#     
#   }
#   
# }
## this is bad because only fits the default cline
# fast_hzar = function(qopt_ada,index=1,add=F){
#   qopt_model <-    hzar.makeCline1DFreq(data=qopt_ada, scaling="fixed",tails="none");
#   qopt_model <-    hzar.model.addBoxReq(meta.model=qopt_model,
#                                         low=0,high=20);
#   qopt_model_FitR <-
#     hzar.first.fitRequest.old.ML(model=qopt_model ,
#                                  obsData=qopt_ada,
#                                  verbose=FALSE);
#   ## this does not work because this just fits a standard cline to the data. it never actually fits it?
#   
#   if(index==1 | add==F){
#     hzar.plot.cline(qopt_model_FitR,xlim=c(0,20),ylim=c(0,1),col=index,
#                     add=F);
#   } else {
#     hzar.plot.cline(qopt_model_FitR,xlim=c(0,20),ylim=c(0,1),col=index,
#                     add=T);
#   }
#   
#   
# }
