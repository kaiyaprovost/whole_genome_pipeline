rm(list=ls())

library(ape)
library(vcfR)
library(adegenet)
library(poppr)
library(phangorn)
library(plotrix)
library(rgdal)
library(RColorBrewer)
library(tools)
library(R.utils)

mypath="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/"
setwd(mypath)
overwrite=T
doPlots=F

specieslist =rev(sort(c("cur","cri","bru","bil","fla",
                        "fus","mel","nit","sin","bel")))
#specieslist=c("sin")


## possibly amphispiza is too large to read 

colors=c("red","blue","orange","green","goldenrod","magenta",
         "black","grey","brown","pink","purple","darkgreen","lightblue")

calcDist = T
reloadDist = F
combineVcfs = F
compareNewick = F

## import color scheme
col = brewer.pal(12,"Paired")
red=col[6]
blue=col[2]
green=col[3]
yellow=col[7]
grey=col[9]

chromlist=c("1","1A","1B","2","3","4","4A","5","6","7","8","9","10",
            "11","12","13","14","15","16","17","18","19","20","21",
            "22","23","24","25","26","27","28","LG2","LGE22","mtDNA","Z")

if (calcDist == T) {
  
  for(rep in c(3,5,6,2)){
    print(rep)
    if(rep==1) {
      allfiles=list.files(path=paste(mypath,"VCFS/",sep=""),pattern="vcf.+converted$",recursive = T,full.names = T) 
      
    } else if (rep==2) {
      allfiles=list.files(path=paste(mypath,"VCFS/",sep=""),pattern="vcf.+converted.+vcf$",recursive = T,full.names = T)
      
    } else if (rep==3) {
      allfiles=list.files(path=paste(mypath,"VCFS/",sep=""),pattern="vcf.+converted.+vcf.gz",recursive = T,full.names = T)
      
    } else if(rep==4) {
      allfiles=list.files(path=paste(mypath,"GZIPPED_VCFS/",sep=""),pattern="vcf.+converted$",recursive = T,full.names = T) 
      
    } else if (rep==5) {
      allfiles=list.files(path=paste(mypath,"GZIPPED_VCFS/",sep=""),pattern="vcf.+converted.+vcf$",recursive = T,full.names = T)
      
    } else {
      allfiles=list.files(path=paste(mypath,"GZIPPED_VCFS/",sep=""),pattern="vcf.+converted.+vcf.gz",recursive = T,full.names = T)
      
    }
    
    #files=c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/GZIPPED_VCFS/SINUATUS/Cardinalis-sinuatus-called.geno.PseudoNC_007897.1_Tgut_mtDNA.vcf.fixedchroms.converted_w100000_o100000_0.window.vcf.gz")
    x <- file.info(allfiles)
    allfiles = allfiles[match(1:length(allfiles),rank(x$size))]
    allfiles = allfiles[complete.cases(allfiles)]
    #allfiles=sample(allfiles)
    
    head(allfiles)
    print(paste("FILES TO RUN:",length(allfiles)))
    
    for (spp in specieslist) {
      print(spp)
      
      files=allfiles[grep(spp,allfiles)]
      
      if (length(files)>0) {
        
        for (vcffile_index in 1:length(files)) {
          vcffile=files[vcffile_index]
          #print(vcffile)
          
          now=Sys.time()

          extension=tools::file_ext(vcffile)
          if(extension == "gz") {
            R.utils::gunzip(vcffile,skip=T)
            vcffile=tools::file_path_sans_ext(vcffile)
          }
          zippedfile = paste(vcffile,".gz",sep="")
          
          print(vcffile)
          print(paste(vcffile_index,"of",length(files)))
          
          outfile = paste(vcffile,"_distances.dist",sep="")
          
          nopathname=basename(vcffile)
          split = strsplit(vcffile,"/")[[1]]
          species=split[11]
          
          print("looking for matches")
          matching=list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/",
                     pattern=paste(species,"/",nopathname,"_distances.dist",sep=""),
                     recursive = T,full.names = T)
          
          distexists=as.logical(length(matching))
                    
          # distexists = as.logical(sum(c(
          #   file.exists(outfile),
          #   file.exists(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/DISTS/",
          #                     species,"/",nopathname,"_distances.dist",sep="")),
          #   file.exists(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/DISTS/",
          #                     nopathname,"_distances.dist",sep="")),
          #   file.exists(paste(outfile,".gz",sep="")),
          #   file.exists(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/DISTS/",
          #                     species,"/",nopathname,"_distances.dist.gz",sep="")),
          #   file.exists(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/DISTS/",
          #                     nopathname,"_distances.dist.gz",sep=""))
          # )))
          
          if(distexists && overwrite==F){
            print("-----distances already calculated")
            file.rename(vcffile,paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/FINISH_DISTS/",nopathname,sep=""))
            
          } else {
            
            outfile=paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/DISTS/",
                          species,"/",nopathname,"_distances.dist",sep="")
            
            print("Reading")
            vcf <- vcfR::read.vcfR(vcffile, verbose = TRUE,limit=1e08)#nrows=3000000)
            
            
            print("Genlight")
            x <- vcfR::vcfR2genlight(vcf)
            
            if (!(max(ploidy(x))==min(ploidy(x)) && !(is.na(max(ploidy(x)))))) {
              print("XXXXXXX SOMETHING WRONG WITH FILE -- SET PLOIDY MANUALLY XXXXXXXXX")
              cat("\n[ Manually Setting Ploidy ]\n", file=outfile,sep="\n",append=TRUE)
              ploidy(x)=rep(2,length(ploidy(x)))
            }
            
            print("Dist 1")
            x.dist <- dist(x)
            cat(vcffile, file=outfile,sep="\n",append=FALSE)
            cat("\nRaw Distance Matrix\n", file=outfile,sep="\n",append=TRUE)
            write.table(as.matrix(x.dist), outfile, row.names=FALSE, col.names=FALSE,append = T)
            print("Dist 2")
            x.dist2 <- poppr::bitwise.dist(x,missing_match = T,percent=F,scale_missing = F)
            print("Dist 3")
            x.dist3 <- poppr::bitwise.dist(x,missing_match = T,percent=F,scale_missing = T)
            print("Dist 4")
            x.dist4 <- poppr::bitwise.dist(x,missing_match = F,percent=F,scale_missing = F)
            print("Dist 5")
            x.dist5 <- poppr::bitwise.dist(x,missing_match = F,percent=F,scale_missing = T)
            
            cat("\nBitwise Matrix Missing Match True Scale Missing False\n", file=outfile,sep="\n",append=TRUE)
            write.table(as.matrix(x.dist2), outfile, row.names=FALSE, col.names=FALSE,append = T)
            cat("\nBitwise Matrix Missing Match True Scale Missing True\n", file=outfile,sep="\n",append=TRUE)
            write.table(as.matrix(x.dist3), outfile, row.names=FALSE, col.names=FALSE,append = T)
            cat("\nBitwise Matrix Missing Match False Scale Missing False\n", file=outfile,sep="\n",append=TRUE)
            write.table(as.matrix(x.dist4), outfile, row.names=FALSE, col.names=FALSE,append = T)
            cat("\nBitwise Matrix Missing Match False Scale Missing True\n", file=outfile,sep="\n",append=TRUE)
            write.table(as.matrix(x.dist5), outfile, row.names=FALSE, col.names=FALSE,append = T)
            
            
            cat("\n#####\n", file=outfile,sep="\n",append=TRUE)
            
            tree1 =tryCatch(njs(x.dist),error=function(cond){return(NA)})
            tree2 =tryCatch(njs(x.dist2),error=function(cond){return(NA)})
            tree3 =tryCatch(njs(x.dist3),error=function(cond){return(NA)})
            tree4 =tryCatch(njs(x.dist4),error=function(cond){return(NA)})
            tree5 =tryCatch(njs(x.dist5),error=function(cond){return(NA)})
            
            if(doPlots==T){
              png(paste(outfile,"_distances.png",sep=""),
                  width=1200,height=800)
              par(mfrow=c(1,5),mar=c(0,0,0,0))
              
              if(!(is.na(tree1))) {
                plot(tree1,tip.color=colors)
              } else {
                barplot(0)
              }
              if(!(is.na(tree2))) {
                plot(tree2,tip.color=colors)
              } else {
                barplot(0)
              }
              if(!(is.na(tree3))) {
                plot(tree3,tip.color=colors)
              } else {
                barplot(0)
              }
              if(!(is.na(tree4))) {
                plot(tree4,tip.color=colors)
              } else {
                barplot(0)
              }
              if(!(is.na(tree5))) {
                plot(tree5,tip.color=colors)
              } else {
                barplot(0)
              }
              dev.off()
              gzip(paste(outfile,"_distances.png",sep=""),paste(outfile,"_distances.png.gz",sep=""),skip=T)
            }
            
            gzip(outfile,paste(outfile,".gz",sep=""),skip=T)
            
          }
          gzip(vcffile,zippedfile,skip=T)    
          
          file.rename(zippedfile,paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/FINISH_DISTS/",nopathname,".gz",sep=""))
          
          rm(vcf)
          rm(x)
          rm(x.dist)
          rm(x.dist2)
          rm(x.dist3)
          rm(x.dist4)
          rm(x.dist5)
          rm(tree1)
          rm(tree2)
          rm(tree3)
          rm(tree4)
          rm(tree5)
        
          print(now-Sys.time())
          }
        
      }
      
      rm(files)
    }
    rm(allfiles)
    
  }
}

if(reloadDist == T) {
  
  path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/DISTS/"
  distfilelist = list.files(path,pattern=".+all.vcf.+distances.dist",recursive = T,full.names = T)
  #distfilelist=c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/Amphispiza-bilineata-called.geno/DISTS/Amphispiza-bilineata-called.geno.PseudoNC.all.vcf.gz_distances.dist")
  
  for (distfile in distfilelist[1:length(distfilelist)]) {
    
    print(distfile)
    
    outpng = paste(substr(distfile,1,nchar(distfile)-5),"_consensus.png",sep="")
    threepng = paste(substr(distfile,1,nchar(distfile)-5),"_tree3_dots.png",sep="")
    threenewick = paste(substr(distfile,1,nchar(distfile)-5),"_tree3.newick",sep="")
    
    if (file.exists(outpng) && overwrite==F) {
      print("---skipping")
    } else {
      
      #distfile = "Vireo-bellii-called.geno/Vireo-bellii-called.geno.PseudoNC_007897.1_Tgut_mtDNA.vcf_distances.dist"
      distlines = readLines(distfile)
      
      ## get the labels
      spp = (strsplit(strsplit(distfile,"/")[[1]][12],"\\.")[[1]][1])
      spp = substr(spp,1,nchar(spp)-7)
      #spp="Vireo-bellii"
      #spp="Amphispiza-bilineata"
      region = (strsplit(strsplit(distfile,"/")[[1]][12],"\\.")[[1]][4])
      
      labels=readLines(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/",spp,".bamlist",sep=""))
      labels=sapply(1:length(labels),FUN=function(x) {
        strsplit(labels[x],"/")[[1]][8]
      }
      )
      labels=sapply(1:length(labels),FUN=function(x) {
        strsplit(labels[x],"\\.")[[1]][1]
      }
      )
      
      spp=sub("-","_",spp)
      
      sonlabels=readLines(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/son/SON_",spp,".indlist",sep=""))
      sonlabels=sapply(1:length(sonlabels),FUN=function(x) {
        strsplit(sonlabels[x],"/")[[1]][8]
      }
      )
      sonlabels=sapply(1:length(sonlabels),FUN=function(x) {
        strsplit(sonlabels[x],"\\.")[[1]][1]
      }
      )
      
      
      ## alter the names
      
      ## match the colors to the red-yellow-blue spectrum QGIS uses
      #d7191c
      #fdae61
      #ffffbf
      #abdda4
      #2b83ba
      
      color=labels
      color[which(!(labels %in% sonlabels))] = "#d7191c" ## before was variable blue and green
      color[which((labels %in% sonlabels))] = "#2b83ba"
      
      labels[which(!(labels %in% sonlabels))] = paste("C",labels[which(!(labels %in% sonlabels))],"C",sep="-")
      labels[which(labels %in% sonlabels)] = paste("S",labels[which(labels %in% sonlabels)],"S",sep="-")
      
      labels = sub("AMN_245","",labels)
      labels = sub("_P002_","-",labels)
      labels = sub("_P001_","-",labels)
      labels = sub("_P01_","-",labels)
      labels = sub(paste("_",strsplit(spp,"_")[[1]][2],sep=""),"",labels)
      labels = sub("W","",labels)
      
      breaks = which(distlines == "")
      ends = breaks-1
      starts = c(1,breaks+1)[1:length(breaks)]
      singlelines = intersect(starts,ends)
      matstarts = starts[!(starts %in% singlelines)]
      matends = ends[!(ends %in% singlelines)]
      
      inds=unique(matends-matstarts)[1]+1
      #toname = paste("indiv",0:(inds-1),sep="")
      
      mat1 = read.table(textConnection(distlines[matstarts[1]:matends[1]])); colnames(mat1)=labels; rownames(mat1) = labels
      mat2 = read.table(textConnection(distlines[matstarts[2]:matends[2]])); colnames(mat2)=labels; rownames(mat2) = labels
      mat3 = read.table(textConnection(distlines[matstarts[3]:matends[3]])); colnames(mat3)=labels; rownames(mat3) = labels
      mat4 = read.table(textConnection(distlines[matstarts[4]:matends[4]])); colnames(mat4)=labels; rownames(mat4) = labels
      mat5 = read.table(textConnection(distlines[matstarts[5]:matends[5]])); colnames(mat5)=labels; rownames(mat5) = labels
      
      tree1 =tryCatch(njs(as.matrix(mat1)),error=function(cond){return(NA)})
      tree2 =tryCatch(njs(as.matrix(mat2)),error=function(cond){return(NA)})
      tree3 =tryCatch(njs(as.matrix(mat3)),error=function(cond){return(NA)})
      tree4 =tryCatch(njs(as.matrix(mat4)),error=function(cond){return(NA)})
      tree5 =tryCatch(njs(as.matrix(mat5)),error=function(cond){return(NA)})
      
      treelist = (list(tree1,tree2,tree3,tree4,tree5))
      
      ## write out tree3 to the distance matrix folder for GDMS
      labels2 = paste("Ind_",seq(0,length(labels)-1),sep="")
      spp= toupper(strsplit(spp,as.character("_"))[[1]][2])
      colnames(mat3)=labels2; rownames(mat3) = labels2
      write.csv(mat3,paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",spp,"/",spp,"_distancematrix_",region,".csv",sep = ""))
      
      
      if (unique(is.na(treelist))==T) {
        print("FAILURE")
      } else {
        class(treelist) = "multiPhylo"
        
        rfdist=RF.dist(treelist)
        rfdist
        
        
        cons = consensus(treelist,p=0.5)
        
        if(doPlots==T){
          png(outpng)
          par(mar=c(0,0,0,0))
          plot(cons,tip.color=color,type="radial")
          dev.off()
          
          png(threepng)
          par(mar=c(0,0,0,0))
          plot(tree3,tip.color=color,type="radial")
          dev.off()
          
          png(threepng)
          par(mar=c(0,0,0,0),
              bg = 'black', fg = 'white',
              col="white",
              font=2)
          plot(tree3,tip.color="black",type="radial",
               edge.color="#ffffbf",show.tip.label=T,edge.width=3)
          tiplabels(pch=16,col=color,cex=3)
          text(0,0,spp,col="darkgrey")
          dev.off()
        }
        
        write.tree(tree3,threenewick)
        
        
      } 
    }
  }
}

if(combineVcfs==T) {
  
  #for (chrom in chromlist) {
  for (chrom in c("mtDNA")) {
    files=list.files(pattern=paste("Tgut_",chrom,".+vcf$",sep=""),recursive = T) 
    #files=c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/Amphispiza-bilineata-called.geno/VCFS/Amphispiza-bilineata-called.geno.PseudoNC.all.vcf.gz")
    x <- file.info(files)
    files = files[match(1:length(files),rank(x$size))]
    #head(files)
    print(paste("NUM FILES:",length(files),sep=" "))
    
    mastervcf = c()
    
    for (vcffile in files) {
      print(vcffile)
      sppname = strsplit(vcffile,"/")[[1]][2]
      sppname = strsplit(sppname,"\\.")[[1]][1]
      sppname = strsplit(sppname,"-")[[1]][1:2]
      sppname = paste(sppname[1],"-",sppname[2],sep="")
      vcf <- vcfR::read.vcfR(vcffile, verbose = TRUE,limit=1e08)#nrows=3000000)
      
      if(is.null(mastervcf)) {
        
      }
      
    }
    
  }
  
}

if(compareNewick==T) {
  print("yes")
  
  path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES"
  setwd(path)
  
  for (spp in specieslist) {
    print(spp)
    outplot = paste(spp,"_corrplot.png",sep="")
    newickfiles1 = list.files(path,pattern="newick",recursive = T,full.names = T)
    newickfiles2 = list.files(path,pattern=spp,recursive = T,full.names = T)
    newickfiles = intersect(newickfiles1,newickfiles2)
    ## need a better system for this
    
    names=sapply(newickfiles,FUN=function(x){strsplit(x,"/")[[1]][11]}) ## before was 12
    names=sapply(names,FUN=function(x){strsplit(x,"\\.")[[1]][4]})
    names=as.character(names)
    
    trees = lapply(newickfiles,FUN=read.tree)
    class(trees) = "multiPhylo"
    phydist=phytools::multiRF(trees)
    colnames(phydist) = names
    rownames(phydist) = names
    phydist=as.matrix(phydist)
    
    phydist=unique(phydist)
    phydist=unique(t(phydist))
    phydist=t(phydist)
    
    spectralpal=colorRampPalette(c("#d7191c","#fdae61","#ffffbf","#abdda4","#2b83ba"))
    
    if(doPlot==T){
      png(outplot)
      corrplot::corrplot(phydist,is.corr=F,diag=T,order="alphabet",method="color",cl.lim=c(0,44),
                         #col=colorRampPalette(c(viridis::viridis(5)))(45),
                         col=c(rep(rgb(0,0,0,0),44),spectralpal(45)),
                         main=paste("\n\n",spp,sep=""))
      dev.off()
    }
    #"Toxostoma-curvirostre-called.geno/DISTS/Toxostoma-curvirostre-called.geno.PseudoNC.all.vcf_distances_tree3.newick"
  }
  
}
