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

specieslist =(sort(c("cur","cri","bru","bil","fla",
                     "fus","mel","nit","sin","bel")))
#specieslist=c("bel")


## possibly amphispiza is too large to read 

colors=c("red","blue","orange","green","goldenrod","magenta",
         "black","grey","brown","pink","purple","darkgreen","lightblue")

calcDist = F
reloadDist = F
combineVcfs = F
compareNewick = T

## import color scheme
col = brewer.pal(12,"Paired")
red=col[6]
blue=col[2]
green=col[3]
yellow=col[7]
grey=col[9]

chromlist=c("1","1A","1B","LG5","mtDNA","16")

if (calcDist == T) {
  
  #  for(rep in c(3,2)){
  
  for (spp in specieslist) {
    print(spp)
    
    for(chrom in chromlist){
      print(chrom)
      #files=allfiles[grep(paste("Tgut_",chrom,".vcf",sep=""),allfiles)]
      
      upperspp = toupper(spp)
      
      paths=c(list.dirs(path=paste(mypath,"VCFS",sep=""),recursive=F),
              list.dirs(path=paste(mypath,"FINISH_DISTS",sep=""),recursive=F) )
      spppath = paths[grepl(upperspp,paths)]
      
      paths2=list.dirs(path=spppath,recursive=F,full.names = F)
      paths3=list.dirs(path=spppath,recursive=F,full.names = T)
      sppchrompath = paths3[paths2==chrom]
      #files=files[grep(spp,files)]
      
      #print(rep)
      #if(rep==1) {
      #allfiles=list.files(path=sppchrompath,pattern=".vcf.+converted",recursive = T,full.names = T) 
      
      #} else if (rep==2) {
      #allfiles=list.files(path=sppchrompath,pattern=paste(spp,".+Tgut_",chrom,".vcf.+converted.+vcf$",sep=""),recursive = T,full.names = T)
      
      #} else if (rep==3) {
      #allfiles=list.files(path=sppchrompath,pattern=paste(spp,".+Tgut_",chrom,".vcf.+converted.+vcf.gz",sep=""),recursive = T,full.names = T)
      
      #} 
      
      files=c(list.files(path=sppchrompath,pattern="vcf.fixedchroms",recursive = T,full.names = T))
      
      ## remove windows
      files = files[!(grepl("window",files))]
      
      #files=c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/GZIPPED_VCFS/SINUATUS/Cardinalis-sinuatus-called.geno.PseudoNC_007897.1_Tgut_mtDNA.vcf.fixedchroms.converted_w100000_o100000_0.window.vcf.gz")
      x <- file.info(files)
      files = files[match(1:length(files),rank(x$size))]
      files = files[complete.cases(files)]
      #files=sample(files)
      
      head(files)
      print(paste("FILES TO RUN:",length(files)))
      
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
          
          if(overwrite == F){
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
            
          } else {distexists=F}
          
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
    
  }
  
  rm(allfiles)
  
}

if(reloadDist == T) {
  
  path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/DISTS"
  distfilelist = list.files(path,pattern="distances.dist",recursive = T,full.names = T)
  #distfilelist=c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/Amphispiza-bilineata-called.geno/DISTS/Amphispiza-bilineata-called.geno.PseudoNC.all.vcf.gz_distances.dist")
  
  distfilelist = distfilelist[!(grepl("window",distfilelist))]

  
  for (distfile in distfilelist[1:length(distfilelist)]) {
    
    print(distfile)
    
    extension=tools::file_ext(distfile)
    if(extension == "gz") {
      R.utils::gunzip(distfile,skip=T)
      distfile=tools::file_path_sans_ext(distfile)
    }
    zippedfile = paste(distfile,".gz",sep="")
    
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
      
      spp =strsplit(basename(distfile),"\\.")[[1]][1]
      #spp = (strsplit(strsplit(distfile,"/")[[1]][12],"\\.")[[1]][1])
      spp = substr(spp,1,nchar(spp)-7)
      #spp="Vireo-bellii"
      #spp="Amphispiza-bilineata"
      
      region=(strsplit(basename(distfile),"\\.")[[1]][4])
      #region = (strsplit(strsplit(distfile,"/")[[1]][12],"\\.")[[1]][4])
      
      
      
      
      
      if(spp!="Auriparus-flaviceps"){
        labels=readLines(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/",spp,".bamlist",sep=""))
        spp=sub("-","_",spp)
        sonlabels=readLines(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/son/SON_",spp,".indlist",sep=""))
        
      } else {
        labels=readLines(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/old/",spp,"-WITHFAILS.bamlist",sep=""))
        spp=sub("-","_",spp)
        sonlabels=readLines(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/son/testing_flaviceps/SON_",spp,"-WITHFAILS.indlist",sep=""))
        
      }
      
      
      labels=sapply(1:length(labels),FUN=function(x) {
        strsplit(labels[x],"/")[[1]][8]
      }
      )
      labels=sapply(1:length(labels),FUN=function(x) {
        strsplit(labels[x],"\\.")[[1]][1]
      }
      )
      
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
      
      
      
      
      #mat1 = read.table(textConnection(distlines[matstarts[1]:matends[1]])); colnames(mat1)=labels; rownames(mat1) = labels
      #mat2 = read.table(textConnection(distlines[matstarts[2]:matends[2]])); colnames(mat2)=labels; rownames(mat2) = labels
      if(sum(nchar(unique(distlines[matstarts[3]:matends[3]])),na.rm=T) > 0){
        mat3 = read.table(textConnection(distlines[matstarts[3]:matends[3]])); colnames(mat3)=labels; rownames(mat3) = labels
      } else {mat3=NULL}
      #mat4 = read.table(textConnection(distlines[matstarts[4]:matends[4]])); colnames(mat4)=labels; rownames(mat4) = labels
      #mat5 = read.table(textConnection(distlines[matstarts[5]:matends[5]])); colnames(mat5)=labels; rownames(mat5) = labels
      
      #tree1 =tryCatch(njs(as.matrix(mat1)),error=function(cond){return(NA)})
      #tree2 =tryCatch(njs(as.matrix(mat2)),error=function(cond){return(NA)})
      if(!(is.null(mat3))){tree3 =tryCatch(njs(as.matrix(mat3)),error=function(cond){return(NA)})}else{tree3=NA}
      
      #tree4 =tryCatch(njs(as.matrix(mat4)),error=function(cond){return(NA)})
      #tree5 =tryCatch(njs(as.matrix(mat5)),error=function(cond){return(NA)})
      
      #treelist = (list(tree1,tree2,tree3,tree4,tree5))
      
      ## write out tree3 to the distance matrix folder for GDMS
      
      if(!(is.null(mat3))){
        labels2 = paste("Ind_",seq(0,length(labels)-1),sep="")
        spp= toupper(strsplit(spp,as.character("_"))[[1]][2])
        colnames(mat3)=labels2; rownames(mat3) = labels2
        write.csv(mat3,paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",spp,"/",spp,"_distancematrix_",region,".csv",sep = ""))
      }
      
      # if (unique(is.na(treelist))==T) {
      #   print("FAILURE")
      # } else {
      #   class(treelist) = "multiPhylo"
      #   
      #   rfdist=RF.dist(treelist)
      #   rfdist
      #   
      #   
      #   cons = consensus(treelist,p=0.5)
      #   
      #   if(doPlots==T){
      #     png(outpng)
      #     par(mar=c(0,0,0,0))
      #     plot(cons,tip.color=color,type="radial")
      #     dev.off()
      #     
      #     png(threepng)
      #     par(mar=c(0,0,0,0))
      #     plot(tree3,tip.color=color,type="radial")
      #     dev.off()
      #     
      #     png(threepng)
      #     par(mar=c(0,0,0,0),
      #         bg = 'black', fg = 'white',
      #         col="white",
      #         font=2)
      #     plot(tree3,tip.color="black",type="radial",
      #          edge.color="#ffffbf",show.tip.label=T,edge.width=3)
      #     tiplabels(pch=16,col=color,cex=3)
      #     text(0,0,spp,col="darkgrey")
      #     dev.off()
      #   }
      if(!(is.na(tree3))) {write.tree(tree3,threenewick)}
      
    } 
  #}
  
  gzip(distfile,zippedfile,skip=T)  
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

## custom corrplot
draw_method_color <- function(coords, fg, bg) {
  symbols(coords, squares = rep(1, nrow(coords)), fg = fg, bg = bg,
          add = TRUE, inches = FALSE)
}
draw_grid <- function(coords, fg) {
  symbols(coords, add = TRUE, inches = FALSE, fg = fg, bg = NA,
          rectangles = matrix(1, nrow = nrow(coords), ncol = 2),
          lwd=0.25)
}

customcorrplot=function(corr, method = c("circle", "square", "ellipse", "number", 
                                         "shade", "color", "pie"), type = c("full", "lower", "upper"), 
                        add = FALSE, col = NULL, bg = "white", title = "", is.corr = TRUE, 
                        diag = TRUE, outline = FALSE, mar = c(0, 0, 0, 0), addgrid.col = NULL, 
                        addCoef.col = NULL, addCoefasPercent = FALSE, order = c("original","AOE", "FPC", "hclust", "alphabet"), 
                        hclust.method = c("complete", "ward", "ward.D", "ward.D2", "single", "average", "mcquitty", 
                                          "median", "centroid"), addrect = NULL, rect.col = "black", 
                        rect.lwd = 2, tl.pos = NULL, tl.cex = 1, tl.col = "red", 
                        tl.offset = 0.4, tl.srt = 90, cl.pos = NULL, cl.lim = NULL, 
                        cl.length = NULL, cl.cex = 0.8, cl.ratio = 0.15, cl.align.text = "c", 
                        cl.offset = 0.5, number.cex = 1, number.font = 2, number.digits = NULL, 
                        addshade = c("negative", "positive", "all"), shade.lwd = 1, 
                        shade.col = "white", p.mat = NULL, sig.level = 0.05, 
                        insig = c("pch", "p-value", "blank", "n", "label_sig"), pch = 4, pch.col = "black", 
                        pch.cex = 3, plotCI = c("n", "square", "circle", "rect"), 
                        lowCI.mat = NULL, uppCI.mat = NULL, na.label = "?", na.label.col = "black", 
                        win.asp = 1, ...) 
{
  method <- match.arg(method)
  type <- match.arg(type)
  order <- match.arg(order)
  hclust.method <- match.arg(hclust.method)
  addshade <- match.arg(addshade)
  insig <- match.arg(insig)
  plotCI <- match.arg(plotCI)
  if (win.asp != 1 && !(method %in% c("circle", "square"))) {
    stop("Parameter 'win.asp' is supported only for circle and square methods.")
  }
  asp_rescale_factor <- min(1, win.asp)/max(1, win.asp)
  stopifnot(asp_rescale_factor >= 0 && asp_rescale_factor <= 
              1)
  if (!is.matrix(corr) && !is.data.frame(corr)) {
    stop("Need a matrix or data frame!")
  }
  if (is.null(addgrid.col)) {
    addgrid.col <- switch(method, color = NA, shade = NA, 
                          "grey")
  }
  if (any(corr[!is.na(corr)] < cl.lim[1]) || any(corr[!is.na(corr)] > cl.lim[2])) {
    stop("color limits should cover matrix")
  }
  if (is.null(cl.lim)) {
    if (is.corr) {
      cl.lim <- c(-1, 1)
    }
    else {
      corr_tmp <- corr
      diag(corr_tmp) <- ifelse(diag, diag(corr_tmp), NA)
      cl.lim <- c(min(corr_tmp, na.rm = TRUE), max(corr_tmp, 
                                                   na.rm = TRUE))
    }
  }
  intercept <- 0
  zoom <- 1
  if (!is.corr) {
    c_max <- max(corr, na.rm = TRUE)
    c_min <- min(corr, na.rm = TRUE)
    if (c_max <= 0) {
      intercept <- -cl.lim[2]
      zoom <- 1/(diff(cl.lim))
    }
    else if (c_min >= 0) {
      intercept <- -cl.lim[1]
      zoom <- 1/(diff(cl.lim))
    }
    else {
      stopifnot(c_max * c_min < 0)
      stopifnot(c_min < 0 && c_max > 0)
      intercept <- 0
      zoom <- 1/max(abs(cl.lim))
    }
    if (zoom == Inf) {
      stopifnot(cl.lim[1] == 0 && cl.lim[2] == 0)
      zoom <- 0
    }
    corr <- (intercept + corr) * zoom
  }
  cl.lim2 <- (intercept + cl.lim) * zoom
  int <- intercept * zoom
  if (is.corr) {
    if (min(corr, na.rm = TRUE) < -1 - .Machine$double.eps^0.75 || 
        max(corr, na.rm = TRUE) > 1 + .Machine$double.eps^0.75) {
      stop("The matrix is not in [-1, 1]!")
    }
  }
  if (is.null(col)) {
    col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061"))(200)
  }
  n <- nrow(corr)
  m <- ncol(corr)
  min.nm <- min(n, m)
  ord <- seq_len(min.nm)
  if (order != "original") {
    ord <- corrplot::corrMatOrder(corr, order = order, hclust.method = hclust.method)
    corr <- corr[ord, ord]
  }
  if (is.null(rownames(corr))) {
    rownames(corr) <- seq_len(n)
  }
  if (is.null(colnames(corr))) {
    colnames(corr) <- seq_len(m)
  }
  apply_mat_filter <- function(mat) {
    x <- matrix(1:n * m, nrow = n, ncol = m)
    switch(type, upper = mat[row(x) > col(x)] <- Inf, lower = mat[row(x) < 
                                                                    col(x)] <- Inf)
    if (!diag) {
      diag(mat) <- Inf
    }
    return(mat)
  }
  getPos.Dat <- function(mat) {
    tmp <- apply_mat_filter(mat)
    Dat <- tmp[is.finite(tmp)]
    ind <- which(is.finite(tmp), arr.ind = TRUE)
    Pos <- ind
    Pos[, 1] <- ind[, 2]
    Pos[, 2] <- -ind[, 1] + 1 + n
    return(list(Pos, Dat))
  }
  getPos.NAs <- function(mat) {
    tmp <- apply_mat_filter(mat)
    ind <- which(is.na(tmp), arr.ind = TRUE)
    Pos <- ind
    Pos[, 1] <- ind[, 2]
    Pos[, 2] <- -ind[, 1] + 1 + n
    return(Pos)
  }
  Pos <- getPos.Dat(corr)[[1]]
  if (any(is.na(corr)) && is.character(na.label)) {
    PosNA <- getPos.NAs(corr)
  }
  else {
    PosNA <- NULL
  }
  AllCoords <- rbind(Pos, PosNA)
  n2 <- max(AllCoords[, 2])
  n1 <- min(AllCoords[, 2])
  nn <- n2 - n1
  m2 <- max(AllCoords[, 1])
  m1 <- min(AllCoords[, 1])
  mm <- max(1, m2 - m1)
  expand_expression <- function(s) {
    ifelse(grepl("^[:=$]", s), parse(text = substring(s, 
                                                      2)), s)
  }
  newrownames <- sapply(rownames(corr)[(n + 1 - n2):(n + 1 - 
                                                       n1)], expand_expression)
  newcolnames <- sapply(colnames(corr)[m1:m2], expand_expression)
  DAT <- getPos.Dat(corr)[[2]]
  len.DAT <- length(DAT)
  rm(expand_expression)
  assign.color <- function(dat = DAT, color = col) {
    newcorr <- (dat + 1)/2
    newcorr[newcorr <= 0] <- 0
    newcorr[newcorr >= 1] <- 1 - 1e-16
    color[floor(newcorr * length(color)) + 1]
  }
  col.fill <- assign.color()
  isFALSE <- function(x) identical(x, FALSE)
  isTRUE <- function(x) identical(x, TRUE)
  if (isFALSE(tl.pos)) {
    tl.pos <- "n"
  }
  if (is.null(tl.pos) || isTRUE(tl.pos)) {
    tl.pos <- switch(type, full = "lt", lower = "ld", upper = "td")
  }
  if (isFALSE(cl.pos)) {
    cl.pos <- "n"
  }
  if (is.null(cl.pos) || isTRUE(cl.pos)) {
    cl.pos <- switch(type, full = "r", lower = "b", upper = "r")
  }
  if (isFALSE(outline)) {
    col.border <- col.fill
  }
  else if (isTRUE(outline)) {
    col.border <- "black"
  }
  else if (is.character(outline)) {
    col.border <- outline
  }
  else {
    stop("Unsupported value type for parameter outline")
  }
  oldpar <- par(mar = mar, bg = "white")
  on.exit(par(oldpar), add = TRUE)
  if (!add) {
    plot.new()
    xlabwidth <- max(strwidth(newrownames, cex = tl.cex))
    ylabwidth <- max(strwidth(newcolnames, cex = tl.cex))
    laboffset <- strwidth("W", cex = tl.cex) * tl.offset
    for (i in 1:50) {
      xlim <- c(m1 - 0.5 - laboffset - xlabwidth * (grepl("l", 
                                                          tl.pos) | grepl("d", tl.pos)), m2 + 0.5 + mm * 
                  cl.ratio * (cl.pos == "r") + xlabwidth * abs(cos(tl.srt * 
                                                                     pi/180)) * grepl("d", tl.pos)) + c(-0.35, 0.15) + 
        c(-1, 0) * grepl("l", tl.pos)
      ylim <- c(n1 - 0.5 - nn * cl.ratio * (cl.pos == "b") - 
                  laboffset, n2 + 0.5 + laboffset + ylabwidth * 
                  abs(sin(tl.srt * pi/180)) * grepl("t", tl.pos)) + 
        c(-0.15, 0) + c(0, -1) * (type == "upper" && 
                                    tl.pos != "n") + c(0, 1) * grepl("d", tl.pos)
      plot.window(xlim, ylim, asp = 1, xaxs = "i", yaxs = "i")
      x.tmp <- max(strwidth(newrownames, cex = tl.cex))
      y.tmp <- max(strwidth(newcolnames, cex = tl.cex))
      laboffset.tmp <- strwidth("W", cex = tl.cex) * tl.offset
      if (max(x.tmp - xlabwidth, y.tmp - ylabwidth, laboffset.tmp - 
              laboffset) < 0.001) {
        break
      }
      xlabwidth <- x.tmp
      ylabwidth <- y.tmp
      laboffset <- laboffset.tmp
      if (i == 50) {
        warning(c("Not been able to calculate text margin, ", 
                  "please try again with a clean new empty window using ", 
                  "{plot.new(); dev.off()} or reduce tl.cex"))
      }
    }
    if (.Platform$OS.type == "windows") {
      grDevices::windows.options(width = 7, height = 7 * 
                                   diff(ylim)/diff(xlim))
    }
    plot.window(xlim = xlim, ylim = ylim, asp = win.asp, 
                xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  }
  laboffset <- strwidth("W", cex = tl.cex) * tl.offset
  symbols(Pos, add = TRUE, inches = FALSE, rectangles = matrix(1, 
                                                               len.DAT, 2), bg = bg, fg = bg)
  if (method == "circle" && plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, circles = asp_rescale_factor * 
              0.9 * abs(DAT)^0.5/2, fg = col.border, bg = col.fill)
  }
  if (method == "ellipse" && plotCI == "n") {
    ell.dat <- function(rho, length = 99) {
      k <- seq(0, 2 * pi, length = length)
      x <- cos(k + acos(rho)/2)/2
      y <- cos(k - acos(rho)/2)/2
      cbind(rbind(x, y), c(NA, NA))
    }
    ELL.dat <- lapply(DAT, ell.dat)
    ELL.dat2 <- 0.85 * matrix(unlist(ELL.dat), ncol = 2, 
                              byrow = TRUE)
    ELL.dat2 <- ELL.dat2 + Pos[rep(1:length(DAT), each = 100), 
                               ]
    polygon(ELL.dat2, border = col.border, col = col.fill)
  }
  if (is.null(number.digits)) {
    number.digits <- switch(addCoefasPercent + 1, 2, 0)
  }
  stopifnot(number.digits%%1 == 0)
  stopifnot(number.digits >= 0)
  if (method == "number" && plotCI == "n") {
    text(Pos[, 1], Pos[, 2], font = number.font, col = col.fill, 
         labels = round((DAT - int) * ifelse(addCoefasPercent, 
                                             100, 1)/zoom, number.digits), cex = number.cex)
  }
  NA_LABEL_MAX_CHARS <- 2
  if (is.matrix(PosNA) && nrow(PosNA) > 0) {
    stopifnot(is.matrix(PosNA))
    if (na.label == "square") {
      symbols(PosNA, add = TRUE, inches = FALSE, squares = rep(1, 
                                                               nrow(PosNA)), bg = na.label.col, fg = na.label.col)
    }
    else if (nchar(na.label) %in% 0:NA_LABEL_MAX_CHARS) {
      symbols(PosNA, add = TRUE, inches = FALSE, squares = rep(1, 
                                                               nrow(PosNA)), fg = bg, bg = bg)
      text(PosNA[, 1], PosNA[, 2], font = number.font, 
           col = na.label.col, labels = na.label, cex = number.cex, 
           ...)
    }
    else {
      stop(paste("Maximum number of characters for NA label is:", 
                 NA_LABEL_MAX_CHARS))
    }
  }
  if (method == "pie" && plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, circles = rep(0.5, 
                                                           len.DAT) * 0.85, fg = col.border)
    pie.dat <- function(theta, length = 100) {
      k <- seq(pi/2, pi/2 - theta, length = 0.5 * length * 
                 abs(theta)/pi)
      x <- c(0, cos(k)/2, 0)
      y <- c(0, sin(k)/2, 0)
      cbind(rbind(x, y), c(NA, NA))
    }
    PIE.dat <- lapply(DAT * 2 * pi, pie.dat)
    len.pie <- unlist(lapply(PIE.dat, length))/2
    PIE.dat2 <- 0.85 * matrix(unlist(PIE.dat), ncol = 2, 
                              byrow = TRUE)
    PIE.dat2 <- PIE.dat2 + Pos[rep(1:length(DAT), len.pie), 
                               ]
    polygon(PIE.dat2, border = "black", col = col.fill)
  }
  if (method == "shade" && plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, squares = rep(1, 
                                                           len.DAT), bg = col.fill, fg = addgrid.col)
    shade.dat <- function(w) {
      x <- w[1]
      y <- w[2]
      rho <- w[3]
      x1 <- x - 0.5
      x2 <- x + 0.5
      y1 <- y - 0.5
      y2 <- y + 0.5
      dat <- NA
      if ((addshade == "positive" || addshade == "all") && 
          rho > 0) {
        dat <- cbind(c(x1, x1, x), c(y, y1, y1), c(x, 
                                                   x2, x2), c(y2, y2, y))
      }
      if ((addshade == "negative" || addshade == "all") && 
          rho < 0) {
        dat <- cbind(c(x1, x1, x), c(y, y2, y2), c(x, 
                                                   x2, x2), c(y1, y1, y))
      }
      return(t(dat))
    }
    pos_corr <- rbind(cbind(Pos, DAT))
    pos_corr2 <- split(pos_corr, 1:nrow(pos_corr))
    SHADE.dat <- matrix(na.omit(unlist(lapply(pos_corr2, 
                                              shade.dat))), byrow = TRUE, ncol = 4)
    segments(SHADE.dat[, 1], SHADE.dat[, 2], SHADE.dat[, 
                                                       3], SHADE.dat[, 4], col = shade.col, lwd = shade.lwd)
  }
  if (method == "square" && plotCI == "n") {
    draw_method_square(Pos, DAT, asp_rescale_factor, col.border, 
                       col.fill)
  }
  if (method == "color" && plotCI == "n") {
    draw_method_color(Pos, col.border, col.fill)
  }
  draw_grid(AllCoords, addgrid.col)
  if (plotCI != "n") {
    if (is.null(lowCI.mat) || is.null(uppCI.mat)) {
      stop("Need lowCI.mat and uppCI.mat!")
    }
    if (order != "original") {
      lowCI.mat <- lowCI.mat[ord, ord]
      uppCI.mat <- uppCI.mat[ord, ord]
    }
    pos.lowNew <- getPos.Dat(lowCI.mat)[[1]]
    lowNew <- getPos.Dat(lowCI.mat)[[2]]
    pos.uppNew <- getPos.Dat(uppCI.mat)[[1]]
    uppNew <- getPos.Dat(uppCI.mat)[[2]]
    if (!method %in% c("circle", "square")) {
      stop("Method shoud be circle or square if drawing confidence intervals.")
    }
    k1 <- (abs(uppNew) > abs(lowNew))
    bigabs <- uppNew
    bigabs[which(!k1)] <- lowNew[!k1]
    smallabs <- lowNew
    smallabs[which(!k1)] <- uppNew[!k1]
    sig <- sign(uppNew * lowNew)
    color_bigabs <- col[ceiling((bigabs + 1) * length(col)/2)]
    color_smallabs <- col[ceiling((smallabs + 1) * length(col)/2)]
    if (plotCI == "circle") {
      symbols(pos.uppNew[, 1], pos.uppNew[, 2], add = TRUE, 
              inches = FALSE, circles = 0.95 * abs(bigabs)^0.5/2, 
              bg = ifelse(sig > 0, col.fill, color_bigabs), 
              fg = ifelse(sig > 0, col.fill, color_bigabs))
      symbols(pos.lowNew[, 1], pos.lowNew[, 2], add = TRUE, 
              inches = FALSE, circles = 0.95 * abs(smallabs)^0.5/2, 
              bg = ifelse(sig > 0, bg, color_smallabs), fg = ifelse(sig > 
                                                                      0, col.fill, color_smallabs))
    }
    if (plotCI == "square") {
      symbols(pos.uppNew[, 1], pos.uppNew[, 2], add = TRUE, 
              inches = FALSE, squares = abs(bigabs)^0.5, bg = ifelse(sig > 
                                                                       0, col.fill, color_bigabs), fg = ifelse(sig > 
                                                                                                                 0, col.fill, color_bigabs))
      symbols(pos.lowNew[, 1], pos.lowNew[, 2], add = TRUE, 
              inches = FALSE, squares = abs(smallabs)^0.5, 
              bg = ifelse(sig > 0, bg, color_smallabs), fg = ifelse(sig > 
                                                                      0, col.fill, color_smallabs))
    }
    if (plotCI == "rect") {
      rect.width <- 0.25
      rect(pos.uppNew[, 1] - rect.width, pos.uppNew[, 2] + 
             smallabs/2, pos.uppNew[, 1] + rect.width, pos.uppNew[, 
                                                                  2] + bigabs/2, col = col.fill, border = col.fill)
      segments(pos.lowNew[, 1] - rect.width, pos.lowNew[, 
                                                        2] + DAT/2, pos.lowNew[, 1] + rect.width, pos.lowNew[, 
                                                                                                             2] + DAT/2, col = "black", lwd = 1)
      segments(pos.uppNew[, 1] - rect.width, pos.uppNew[, 
                                                        2] + uppNew/2, pos.uppNew[, 1] + rect.width, 
               pos.uppNew[, 2] + uppNew/2, col = "black", lwd = 1)
      segments(pos.lowNew[, 1] - rect.width, pos.lowNew[, 
                                                        2] + lowNew/2, pos.lowNew[, 1] + rect.width, 
               pos.lowNew[, 2] + lowNew/2, col = "black", lwd = 1)
      segments(pos.lowNew[, 1] - 0.5, pos.lowNew[, 2], 
               pos.lowNew[, 1] + 0.5, pos.lowNew[, 2], col = "grey70", 
               lty = 3)
    }
  }
  if (!is.null(p.mat) && insig != "n") {
    if (order != "original") {
      p.mat <- p.mat[ord, ord]
    }
    pos.pNew <- getPos.Dat(p.mat)[[1]]
    pNew <- getPos.Dat(p.mat)[[2]]
    if (insig == "label_sig") {
      if (!is.character(pch)) 
        pch <- "*"
      place_points <- function(sig.locs, point) {
        text(pos.pNew[, 1][sig.locs], pos.pNew[, 2][sig.locs], 
             labels = point, col = pch.col, cex = pch.cex, 
             lwd = 2)
      }
      if (length(sig.level) == 1) {
        place_points(sig.locs = which(pNew < sig.level), 
                     point = pch)
      }
      else {
        l <- length(sig.level)
        for (i in seq_along(sig.level)) {
          iter <- l + 1 - i
          pchTmp <- paste(rep(pch, i), collapse = "")
          if (i == length(sig.level)) {
            locs <- which(pNew < sig.level[iter])
            if (length(locs)) {
              place_points(sig.locs = locs, point = pchTmp)
            }
          }
          else {
            locs <- which(pNew < sig.level[iter] & pNew > 
                            sig.level[iter - 1])
            if (length(locs)) {
              place_points(sig.locs = locs, point = pchTmp)
            }
          }
        }
      }
    }
    else {
      ind.p <- which(pNew > sig.level)
      p_inSig <- length(ind.p) > 0
      if (insig == "pch" && p_inSig) {
        points(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
               pch = pch, col = pch.col, cex = pch.cex, lwd = 2)
      }
      if (insig == "p-value" && p_inSig) {
        text(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
             round(pNew[ind.p], 2), col = pch.col)
      }
      if (insig == "blank" && p_inSig) {
        symbols(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
                inches = FALSE, squares = rep(1, length(pos.pNew[, 
                                                                 1][ind.p])), fg = addgrid.col, bg = bg, add = TRUE)
      }
    }
  }
  if (cl.pos != "n") {
    colRange <- assign.color(dat = cl.lim2)
    ind1 <- which(col == colRange[1])
    ind2 <- which(col == colRange[2])
    colbar <- col[ind1:ind2]
    if (is.null(cl.length)) {
      cl.length <- ifelse(length(colbar) > 20, 11, length(colbar) + 
                            1)
    }
    labels <- seq(cl.lim[1], cl.lim[2], length = cl.length)
    if (cl.pos == "r") {
      vertical <- TRUE
      xlim <- c(m2 + 0.5 + mm * 0.02, m2 + 0.5 + mm * cl.ratio)
      ylim <- c(n1 - 0.5, n2 + 0.5)
    }
    if (cl.pos == "b") {
      vertical <- FALSE
      xlim <- c(m1 - 0.5, m2 + 0.5)
      ylim <- c(n1 - 0.5 - nn * cl.ratio, n1 - 0.5 - nn * 
                  0.02)
    }
    corrplot::colorlegend(colbar = colbar, labels = round(labels, 2), 
                          offset = cl.offset, ratio.colbar = 0.3, cex = cl.cex, 
                          xlim = xlim, ylim = ylim, vertical = vertical, align = cl.align.text)
  }
  if (tl.pos != "n") {
    pos.xlabel <- cbind(m1:m2, n2 + 0.5 + laboffset)
    pos.ylabel <- cbind(m1 - 0.5, n2:n1)
    if (tl.pos == "td") {
      if (type != "upper") {
        stop("type should be \"upper\" if tl.pos is \"dt\".")
      }
      pos.ylabel <- cbind(m1:(m1 + nn) - 0.5, n2:n1)
    }
    if (tl.pos == "ld") {
      if (type != "lower") {
        stop("type should be \"lower\" if tl.pos is \"ld\".")
      }
      pos.xlabel <- cbind(m1:m2, n2:(n2 - mm) + 0.5 + laboffset)
    }
    if (tl.pos == "d") {
      pos.ylabel <- cbind(m1:(m1 + nn) - 0.5, n2:n1)
      pos.ylabel <- pos.ylabel[1:min(n, m), ]
      symbols(pos.ylabel[, 1] + 0.5, pos.ylabel[, 2], add = TRUE, 
              bg = bg, fg = addgrid.col, inches = FALSE, squares = rep(1, 
                                                                       length(pos.ylabel[, 1])))
      text(pos.ylabel[, 1] + 0.5, pos.ylabel[, 2], newcolnames[1:min(n, 
                                                                     m)], col = tl.col, cex = tl.cex, ...)
    }
    else {
      text(pos.xlabel[, 1], pos.xlabel[, 2], newcolnames, 
           srt = tl.srt, adj = ifelse(tl.srt == 0, c(0.5, 
                                                     0), c(0, 0)), col = tl.col, cex = tl.cex, offset = tl.offset, 
           ...)
      text(pos.ylabel[, 1], pos.ylabel[, 2], newrownames, 
           col = tl.col, cex = tl.cex, pos = 2, offset = tl.offset, 
           ...)
    }
  }
  title(title, ...)
  if (!is.null(addCoef.col) && method != "number") {
    text(Pos[, 1], Pos[, 2], col = addCoef.col, labels = round((DAT - 
                                                                  int) * ifelse(addCoefasPercent, 100, 1)/zoom, number.digits), 
         cex = number.cex, font = number.font)
  }
  if (type == "full" && plotCI == "n" && !is.null(addgrid.col)) {
    rect(m1 - 0.5, n1 - 0.5, m2 + 0.5, n2 + 0.5, border = addgrid.col)
  }
  if (!is.null(addrect) && order == "hclust" && type == "full") {
    corrRect.hclust(corr, k = addrect, method = hclust.method, 
                    col = rect.col, lwd = rect.lwd)
  }
  invisible(corr)
}



if(compareNewick==T) {
  print("yes")
  
  path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/NEWICK"
  setwd(path)
  
  allnewicks = c(1,"1A","1B",2:4,"4A",5:28,"mtDNA","LG2","LG5","LGE22","Z","all")
  
  newickfiles1 = list.files(path,pattern="newick$",recursive = T,full.names = T)
  newickfiles1 = newickfiles1[!grepl("window",newickfiles1)]

  if(doPlots==F){
    png("all_species_rfdist_corrplots.png",#height=600,width=1000,
        height=2.5,width=6.5,units="in",res=600)
    par(mfrow=c(2,5),mar=c(0,0,0,0))
  } 
  for (spp in specieslist) {
    print(spp)
    #outplot = paste(spp,"_corrplot.png",sep="")
    
    #newickfiles2 = list.files(path,pattern=spp,recursive = T,full.names = T)
    #newickfiles = intersect(newickfiles1,newickfiles2)
    ## need a better system for this
    
    newickfiles= newickfiles1[grepl(spp,newickfiles1)]

    #names=sapply(newickfiles,FUN=function(x){strsplit(x,"/")[[1]][11]}) ## before was 12
    names=basename(newickfiles)
    names=sapply(names,FUN=function(x){strsplit(x,"\\.")[[1]][4]})
    names=sapply(names,FUN=function(x){y=strsplit(x,"\\_")[[1]];return(y[length(y)])})
    names=as.character(names)
    #print(names)
    
    newickfiles=newickfiles[which(names %in% allnewicks)]
    names=names[which(names %in% allnewicks)]
    
    newickfiles=newickfiles[!duplicated(names)]
    names=names[!duplicated(names)]
    
    missingcolnames=allnewicks[which(!(allnewicks %in% names))]
    
    trees = lapply(newickfiles,FUN=read.tree)
    class(trees) = "multiPhylo"
    phydist=phytools::multiRF(trees)
    colnames(phydist) = names
    rownames(phydist) = names
    phydist=as.matrix(phydist)
    
    if(nrow(phydist)!=length(allnewicks)){
      rowtoadd=rep(NA,ncol(phydist))
      for(missing in missingcolnames){
        #print(missing)
        phydist=rbind(phydist,rowtoadd)
      }
      coltoadd=rep(NA,nrow(phydist))
      for(missing in missingcolnames){
        #print(missing)
        phydist=cbind(phydist,coltoadd)
      }
      
      names=c(names,missingcolnames)
      colnames(phydist) = names
      rownames(phydist) = names
      
    }
    
    
    #phydist=unique(phydist)
    #phydist=unique(t(phydist))
    #phydist=t(phydist)
    

    
    phydist = phydist[,gtools::mixedorder(colnames(phydist))]
    phydist = phydist[gtools::mixedorder(rownames(phydist)),]
    
    allcol=which(colnames(phydist)=="all")
    neworder=c(1:(allcol-1),(allcol+1):length(colnames(phydist)),allcol)
    phydist=phydist[neworder,neworder]
    
    spectralpal=colorRampPalette(c("#d7191c","#fdae61","#ffffbf","#abdda4","#2b83ba"))
    
    if(doPlots==T){
      png(outplot)
      corrplot::corrplot(phydist,is.corr=F,diag=T,order="original",method="color",
                         cl.lim=c(0,44),
                         #col=colorRampPalette(c(viridis::viridis(5)))(45),
                         col=c(rep(rgb(0,0,0,0),44),spectralpal(45)),
                         main=paste("\n\n",spp,sep=""),cl.length=5)
      dev.off()
    } else {
      customcorrplot(phydist,is.corr=F,diag=T,order="original",method="color",
                     cl.lim=c(0,44),
                     #col=colorRampPalette(c(viridis::viridis(5)))(45),
                     col=c(rep(rgb(0,0,0,0),44),spectralpal(45)),
                     main="",mar=c(0,0,0,0),na.label="",
                     #addgrid.col="lightgrey",
                     tl.offset=0.1,tl.srt=90,
                     tl.cex=0.25,cl.cex=0.5,cl.length=5)
    }
    #"Toxostoma-curvirostre-called.geno/DISTS/Toxostoma-curvirostre-called.geno.PseudoNC.all.vcf_distances_tree3.newick"
  }
  if(doPlots==F){
    dev.off()
  }
}

