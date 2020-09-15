library(hzar)
library(RcppCNPy)

setwd("~")

## take in arguments

species="bellii"
plotWholeSpecies=T

## set up functions

run_hzar=function(qopt_ada,plot=T,index=1){
  
  
  if(plot==T){hzar.plot.obsData(qopt_ada)}
  
  qopt_model <-
    hzar.makeCline1DFreq(data=qopt_ada, scaling="fixed",tails="none");
  
  qopt_model <-
    hzar.model.addBoxReq(meta.model=qopt_model,
                         low=0,high=20);
  
  qopt_model_FitR <-
    hzar.first.fitRequest.old.ML(model=qopt_model ,
                                 obsData=qopt_ada,
                                 verbose=FALSE);
  
  qopt_model_FitR$mcmcParam$chainLength <- 1e5;
  qopt_model_FitR$mcmcParam$burnin <- 5e2;
  qopt_model_Fit <- hzar.doFit(qopt_model_FitR)
  if(plot==T){plot(hzar.mcmc.bindLL(qopt_model_Fit))}
  qopt_model_data <-
    hzar.dataGroup.add(qopt_model_Fit);
  
  qopt_model_data <-
    hzar.dataGroup.add(
      qopt_model_data,
      hzar.chain.doSeq(hzar.next.fitRequest(qopt_model_Fit)));
  hzar.plot.cline(qopt_model_data,xlim=c(0,20),ylim=c(0,1),col=index,
                  add=index!=1);
  if(plot==T){hzar.plot.fzCline(qopt_model_data); ## won't work with just raw longitude
  }
  
  if(plot==T){print(hzar.getLLCutParam(qopt_model_data,c("center","width")))};
  qopt_model_null <- hzar.dataGroup.null(qopt_ada);
  qopt_null_ada <- list(clineModel = qopt_model_data,
                        nullModel  = qopt_model_null);
  qopt_null_ada_2 <- hzar.make.obsDataGroup(qopt_null_ada);
  qopt_null_ada_2 <- hzar.copyModelLabels(qopt_null_ada,qopt_null_ada_2);
  if(plot==T){hzar.plot.cline(qopt_null_ada_2,xlim=c(0,20),ylim=c(0,1),col=index,
                              add=index!=1)
    print(hzar.AICc.hzar.obsDataGroup(qopt_null_ada_2))
  }
  
  res<-list()
  res$mean_width<-mean(qopt_model_Fit$mcmcRaw[,2])
  res$out<-qopt_model_Fit 
  res$dat<-qopt_model_data
  return(res)
  
}

npyOutputFunction = function(i,smallsppfiles,smallspp=NULL,all_qopt){
  
  x=smallsppfiles[i]
  print(x)
  
  if(is.null(smallspp)){
    smallspp=strsplit(basename(x),"\\.")[[1]][1]
  }
  
  #spp=paste(strsplit(strsplit(strsplit(basename(x),"\\.")[[1]][5],"_")[[1]][1],"-")[[1]][-1],sep="-")
  
  splits=strsplit(basename(x),"\\.")[[1]]
  k=as.numeric(splits[3])+1
  spp=paste(strsplit(smallspp,"-")[[1]][1:2],sep="-",collapse="-")
  
  if(grepl("CHI",smallspp,fixed=T)){
    desert="CHI"
  } else if(grepl("SON",smallspp,fixed=T)){
    desert="SON"
  } else{
    desert=""
  }
  
  if (grepl("WEIRD",smallspp,fixed=T)){weird="-NOWEIRD"} else {weird=""}
  
  individual_qopt=npyLoad(x)
  
  colnames(individual_qopt)=paste("CLUSTER",seq(1:ncol(individual_qopt)),sep="")
  all_qopt = all_qopt[order(all_qopt$RG),]
  add_qopt = all_qopt[!(is.na(all_qopt$K2_A)),]
  #add_qopt = all_qopt[,c("RG","SP","LAT","LONG","distance","nSamples","DES")]
  
  if(desert != ""){
    add_qopt = add_qopt[add_qopt$DES==desert,]
  } 
  
  add_qopt=add_qopt[add_qopt$SP==paste(species,weird,sep=""),]
  
  add_qopt = cbind(add_qopt,individual_qopt)
  
  qopt_ada = hzar.doMolecularData1DPops(add_qopt$distance,
                                        add_qopt$CLUSTER1,
                                        add_qopt$nSamples)
  
  
  is_positive=(cor(qopt_ada$frame[,1:2],use="pairwise.complete.obs")[2])>0
  if(is_positive==F){
    qopt_ada = hzar.doMolecularData1DPops(add_qopt$distance,
                                          add_qopt$CLUSTER2,
                                          add_qopt$nSamples)
  }
  
  res=run_hzar(qopt_ada,plot=F,index=i)
  
  output=list()
  mean_width=res$mean_width
  center = mean(res$out$mcmcRaw[,1],na.rm=T)
  output$mean_width = mean_width
  output$center = center
  #output=c(mean_width,center)
  return(output)
  
}

## load lat long etc data

all_qopt=read.table("AllSpeciesMetadata_allK.csv",
                    sep=",",header=T)
all_qopt = all_qopt[order(all_qopt$RG),]
all_qopt$distance = all_qopt$LONG-min(all_qopt$LONG,na.rm=T)
all_qopt$nSamples = 1

## iterate through species files

files = list.files("ADMIX_NUMPY_AUG_2020/",
                   pattern="1.admix.Q.npy",full.names = T)
sppfiles=(files[grepl(species,files)])

if(length(sppfiles)>=1){
  
  allsplits=unlist(lapply(sppfiles,FUN=function(x){strsplit(basename(x),"\\.")[[1]][1]}))
  
  if(plotWholeSpecies==T){pdf(paste(species,"_hzar_npy_chroms.pdf",sep=""))
    
    output=lapply(1:length(sppfiles),FUN = function(i){npyOutputFunction(i,smallsppfiles=sppfiles,smallspp=NULL,all_qopt)})
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
    }
  }
  
} 
dev.off()






