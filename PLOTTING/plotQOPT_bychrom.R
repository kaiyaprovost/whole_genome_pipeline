library(plotrix)
library(rgdal)
library(RColorBrewer)

reloadEnv = F

## import color scheme
col = brewer.pal(12,"Paired")
red=col[6]
blue=col[2]
green=col[3]
yellow=col[7]
grey=col[9]

generatePiePlot = function(bg,spp,ourcols,structure_big=chrom_merge,k=k,chrom_lab=chrom_lab) {
  raster::plot(bg,
               col = "grey",
               colNA = "white",
               legend = F,
               xaxt="n",yaxt="n",
               main = paste(spp,k,chrom_lab)
  )
  
  plot(border, add = T)
  points(structure_big$JLONG, structure_big$JLAT, pch = 4)
  for (i in 1:nrow(structure_big)) {
    
    if(!(is.na(structure_big[i,ourcols[1]]))) {
      
      print(i)
      floating.pie(
        xpos = structure_big$JLONG[i],
        ypos = structure_big$JLAT[i],
        x = c(as.numeric(structure_big[i,c(ourcols)])),
        radius = 0.3,
        col = c("#1F78B4","#B2DF8A","#E31A1C","#FDBF6F","#CAB2D6")
      )
    }
  }
}
generateStrPlot = function(metadata_cols=1:12,structure=structure_all,spp="bellii",k=NULL,k_cols=NULL) {
  
  ## if k is set default to some columns
  if (is.null(k)) {
    k = length(k_cols)
  } else if (k == 2) {
    k_cols = c(42:43)
  } else if (k==3) {
    k_cols = c(44:46)
  } else if (k==4) {
    k_cols = c(47:50)
  } else if (k==5) {
    k_cols = c(51:55)
  } 
  
  structure = structure[order(structure$LONG,structure$LAT),]
  structure_big = structure[structure$SP==spp,]
  kdf = structure_big[,c(metadata_cols,k_cols,ncol(structure_big))]
  kdfm = t(structure_big[,c(k_cols)])
  
  if(ncol(kdfm) >= 1) {
  two = barplot(kdfm,col=c("#1F78B4","#B2DF8A","#E31A1C","#FDBF6F","#CAB2D6"),axisnames=F,axes=F,xpd=F)
  #mtext(text = kdf$LABEL, side = 1, at = two, line = 0.1,las=2,cex=0.5)
  } else {
    barplot(0,main="FAIL")
  }
}
generatePiePlot_shapefile = function(shp,spp,ourcols,structure_big=chrom_merge,k=k,chrom_lab=chrom_lab,main) {
  plot(shp,ylim=c(27,35),xlim=c(-117,-95),col="grey",main=main)
  points(structure_big$JLONG, structure_big$JLAT, pch = 4)
  for (i in 1:nrow(structure_big)) {
    
    if(!(is.na(structure_big[i,ourcols[1]]))) {
      
      print(i)
      floating.pie(
        xpos = structure_big$JLONG[i],
        ypos = structure_big$JLAT[i],
        x = c(as.numeric(structure_big[i,c(ourcols)])),
        radius = 0.3,
        col = c("#1F78B4","#B2DF8A","#E31A1C","#FDBF6F","#CAB2D6")
      )
    }
  }
}


if (reloadEnv == T) {
  Env = raster::stack('/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/bio_2-5m_bil/bio1.bil'
                      #list.files(
                      #path = '/Users/kprovost/Dropbox (AMNH)/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/',
                      #pattern = "\\.bil$",
                      #full.names = T
                      #)[1]
  )
  #ext = raster::extent(c(-125,-60, 10, 50)) ## make sure this will play nice with your points
  #ext = raster::extent(c(-115,-97,26,37))
  ext = raster::extent(c(-117, -98, 27, 35))
  Env = raster::crop(Env, ext)
  bg = Env[[1]] ## just for plotting
  
  border = rgdal::readOGR(
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/mexstates/mexico_and_states_subset.shp"
  )
  #plot(border)
  
}

path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/pca qopt files/chroms/clumpped/by_chrom/"
setwd(path)

structure_all = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/AllSpeciesMetadata_allK_9june2020.csv")
structure_all = structure_all[order(structure_all$RG),]

labelfilelist = list.files(pattern=".chromosomes.order.txt")
datafilelist = sub("chromosomes.order.txt","perm_datafile",labelfilelist)

for (j in 1:length(labelfilelist)) {
  
  labelfile = labelfilelist[j]
  datafile  = datafilelist[j]
  print(labelfile)
  dat    = readLines(datafile)
  labels = readLines(labelfile)
  splitfile = strsplit(labelfile,"\\.")[[1]]
  spp = splitfile[1]
  k = splitfile[2]
  numk = as.numeric(substr(k,2,10))
  structure_spp = structure_all[structure_all$SP==spp,1:12]
  
  breaks = which(dat == "")
  ends = breaks-1
  starts = breaks-(ends[1])
  
  numchrom = length(starts)
  
  ## get chromosome
  
  pdf(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/by_chrom_images/",spp,".",k,".allchrom.pdf",sep=""),
      width=11,height=7)
  par(mfcol=c(2,3),mar=c(1,1,1,1))
  for (i in 1:numchrom) {
    
    chrom_start = starts[i]
    chrom_end = ends[i]
    chrom_lab = strsplit(labels[i],"\\.")[[1]][2]
    
    chrom_data = dat[c(chrom_start:chrom_end)]
    chrom_df = read.table(textConnection(chrom_data))
    chrom_df = chrom_df[,c(-1:-5)]
    chrom_merge = cbind(chrom_df,structure_spp[,1:12])
    
    ## with all combined
    generatePiePlot(bg=bg,spp=spp,ourcols=1:numk,structure_big=chrom_merge,k=k,chrom_lab=chrom_lab)
    
    ## do the str plot 
    chrom_merge$LABEL = paste(chrom_merge$STATE,chrom_merge$CO,sep="")
    generateStrPlot(structure=chrom_merge,spp=spp,k_cols=1:numk)
    
    
    
  }
  dev.off()
  
  ## chromosome averages
  pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_AllC.pdf",
      width=9.5,height=4)
  par(mfrow=c(2,2),mar=c(0,0,1,0))
  for (s in 1:length(unique(structure_all$SP))) {
    
    spp = as.character(unique(structure_all$SP)[s])
    print(spp)
    structure_big = structure_all[structure_all$SP == spp, ]
    
    for (k in 2:5) {
      
      if(k == 2) {
        ourcols = c(42:43)
      } else if (k==3) {
        ourcols = c(44:46)
      } else if (k==4) {
        ourcols = c(47:50)
      } else if (k==5) {
        ourcols = c(51:55)
      } 
      
      
      
      raster::plot(bg,
                   col = "grey",
                   colNA = "white",
                   legend = F,
                   xaxt="n",yaxt="n",
                   main = spp#,xlim=c(-117, -98),ylim=c(27, 35)
                   
      )
      plot(border, add = T)
      points(structure_big$JLONG, structure_big$JLAT, pch = 4)
      for (i in 1:nrow(structure_big)) {
        
        if(!(is.na(structure_big[i,ourcols[1]]))) {
          
          print(i)
          floating.pie(
            xpos = structure_big$JLONG[i],
            ypos = structure_big$JLAT[i],
            x = c(as.numeric(structure_big[i,c(ourcols)])),
            radius = 0.3,
            col = c(blue,green,red,yellow,grey)
          )
        }
      }
      #dev.off()
      
      
    }
  }
  dev.off()
  
}


## doing this again one more time
path="/Users/kprovost/Downloads/by_chrom"
labelfilelist = list.files(path=path,pattern=".chromosomes.order.txt",full.names = T)
labelfilelist = labelfilelist[!(grepl("K4",labelfilelist))]
labelfilelist = labelfilelist[!(grepl("K5",labelfilelist))]
labelfilelist = labelfilelist[!(grepl("K6",labelfilelist))]

datafilelist = sub("chromosomes.order.txt","perm_datafile.reclump.txt",labelfilelist)


structure_all=read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/AllSpeciesMetadata_allK_9june2020.csv",
                       header=T,sep="\t")
structure_all = structure_all[order(structure_all$RG),]

shp1 = shapefile("/Users/kprovost/Downloads/enm_layers/mexstates/mexstates.shp")
shp2 = shapefile("/Users/kprovost/Downloads/enm_layers/cb_2016_us_state_500k/cb_2016_us_state_500k.shp")
shp=bind(shp1,shp2)
#plot(shp,ylim=c(25,35),xlim=c(-120,-90))

pdf("test.pdf",width=7,height=7)
par(mfrow=c(5,2),mar=c(0,0,0,0))
for (j in 1:length(labelfilelist)) {
  
  labelfile = labelfilelist[j]
  datafile  = datafilelist[j]
  print(labelfile)
  dat    = readLines(datafile)
  labels = readLines(labelfile)
  splitfile = strsplit(basename(labelfile),"\\.")[[1]]
  spp = splitfile[1]
  k = splitfile[2]
  numk = as.numeric(substr(k,2,10))
  structure_spp = structure_all[structure_all$SP==spp,1:12]
  
  breaks = which(trimws(dat) == "")
  ends = breaks-1
  starts = breaks-(ends[1])
  
  numchrom = length(starts)
  numinds = ends[1]
  
  #sppdf_raster = raster()
  sppdf = NULL
  
  for(start_index in 1:length(starts)){
    #print(paste("start:",start_index))
    this_start = starts[start_index]
    this_end = ends[start_index]
    chrom_data = dat[c(this_start:this_end)]
    which(trimws(chrom_data) == "")
    
    chrom_df = read.table(textConnection(chrom_data),sep = "\t")
    chrom_df = chrom_df[,c(-1:-5)]
    chrom_df = as.matrix(chrom_df)
    colnames(chrom_df) = 1:ncol(chrom_df)
    rownames(chrom_df) = 1:nrow(chrom_df)
    
    if(is.null(sppdf)){
      sppdf = chrom_df
    } else {
    sppdf = cbind(sppdf,chrom_df)
    }
    
    #sppdf_raster = stack(sppdf_raster,raster(chrom_df))
    
  }
  #par(mfrow=c(2,2))
  #plot(sppdf_raster[[1:3]])
  #plot(quan[[3]])
  
  finalquans = NULL
  
  for(name in as.character(1:numk)){
    thesecols = sppdf[,colnames(sppdf)==name]
    #thesequans = matrixStats::rowQuantiles(thesecols,probs=c(0.25,0.5,0.75))
    thesemeans= as.data.frame(rowMeans(thesecols,na.rm = T))
    colnames(thesemeans) = "Mean"
    #thesequans = cbind(thesequans,thesemeans)
    thesequans = thesemeans
    if(is.null(finalquans)) {
      finalquans = thesequans
    } else {
      finalquans = cbind(finalquans,thesequans)
    }
  }
  
  #plotquan=finalquans[,colnames(finalquans)=="Mean"]
  chrom_merge=cbind(finalquans,structure_spp[,c("JLAT","JLONG")])
  ourcols=1:ncol(finalquans)
  generatePiePlot_shapefile(shp,spp,ourcols,structure_big=chrom_merge,k=k,chrom_lab=chrom_lab,main=quanname)
  
  
  #png(paste(spp,"_",k,"_","chromquan.png",sep=""),height=800,width=1200)
  #par(mfrow=c(2,2))
  #for(quanname in unique(colnames(finalquans))){
  #  plotquan = finalquans[,colnames(finalquans)==quanname]
  #  chrom_merge=cbind(plotquan,structure_spp[,c("JLAT","JLONG")])
  #  ourcols=1:ncol(plotquan)
  #  generatePiePlot_shapefile(shp,spp,ourcols,structure_big=chrom_merge,k=k,chrom_lab=chrom_lab,main=quanname)
  #}
  #dev.off()
  
}
dev.off()

