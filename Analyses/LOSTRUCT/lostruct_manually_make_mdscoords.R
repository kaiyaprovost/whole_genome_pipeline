library(lostruct)
library(colorspace)
library(jsonlite)
library(RColorBrewer)

# distance matrix

knum=2
npc=2
maxrows=11000

## time estimate: seconds ~ 3.022e-05*maxrows^2 + -0.003461*maxrows +0.7488 
## O(n^2)

pca_list=list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/LOSTRUCT/lostruct_results/duplicated",
                    pattern=".pca.csv",full.names = T,recursive=T)
pca_list = pca_list[!grepl("fla",pca_list)]

regions_list = gsub(".pca.csv",".regions.csv",pca_list)
coords_list = gsub(".pca.csv",".coords.csv",pca_list)
distmat_list = gsub(".pca.csv",".distmat.png",pca_list)
bychrom_list = gsub(".pca.csv",".mds_bychrom.png",pca_list)
mdspng_list = gsub(".pca.csv",".mds.png",pca_list)
finalpng_list = gsub(".pca.csv",".final_mds.png",pca_list)

for (file_index in 1:length(pca_list)) {
  print(file_index)
  
  pca_file = pca_list[file_index]
  regions_file = regions_list[file_index]
  mds.file = coords_list[file_index]
  pca_distmat_file = distmat_list[file_index]
  bychrom_file = bychrom_list[file_index]
  mdspng = mdspng_list[file_index]
  finalpng = finalpng_list[file_index]
  
  all.pcas=read.table(pca_file,fill=T,header=T,sep=",")
  all.regions=read.table(regions_file,fill=T,header=T,sep=",")
  all.lengths = as.data.frame(table(all.regions[,1]))
  
  if(nrow(all.pcas[complete.cases(all.pcas),])<=0){
    print("NO DATA")
  } else {
  
  subset=1:maxrows
  
  #pc.distmat <- pc_dist(as.matrix(all.pcas[subset,]), npc=npc) ## 500 ~ 10 seconds?
  #na.inds <- is.na( all.pcas[subset,1] ) # there may be windows with missing data
  #plot(pc.distmat)
  #pc.distmat.nomissing=pc.distmat[!na.inds,!na.inds]
  #plot(pc.distmat.nomissing)
  
  pca.regions = cbind(all.regions,all.pcas)
  head(pca.regions)
  
  #pca.regions.nomissing = pca.regions[!na.inds,]
  
  complete.pca = pca.regions[complete.cases(pca.regions),]
  

    if (nrow(complete.pca) < maxrows){
      subset=1:nrow(complete.pca)
      cat("NEW SUBSET:",nrow(complete.pca),"\n")
    }
    
    system.time(pc.distmat2 <- pc_dist(as.matrix(complete.pca[subset,4:ncol(complete.pca)]), npc=npc)) ## long step
    
    ## seconds per rows
    # rows=c(100,200,250,300,400,500,
    #        600,700,800,1000,2000)
    # seconds=c(0.279,1.179,1.819,2.546,4.298,7.029,
    #           9.692,13.734,16.919,26.856,114.814)
    # plot(rows,seconds,pch=16)
    # mod=lm(seconds~rows)
    # abline(mod)
    # mod=lm(seconds~I(rows^2)+rows)
    # rowvalues <- seq(1, max(subset), 5)
    # predicted <- (predict(mod,list(rows=rowvalues)))
    # lines(rowvalues, predicted,lwd=2, col = "red", xlab = "Time (s)", ylab = "Counts")
    # 
    # predict(mod,list(rows=100000))
    ## an estimate for 100000 rows of data is 83 hours 
    
    #plot(pc.distmat)
    png(pca_distmat_file)
    plot(pc.distmat2)
    dev.off()
    
    #cmdscale_object = cmdscale(pc.distmat.nomissing, k=knum)[ ifelse( na.inds, NA, cumsum(!na.inds) ), ]
    #cmdscale_object = cmdscale(pc.distmat.nomissing, k=knum)
    cmdscale_object = cmdscale(pc.distmat2, k=knum) ## sort of long
    
    colnames(cmdscale_object) <- paste0("MDS",seq_len(knum))
    
    #mds.coords=cbind(pca.regions.nomissing[subset,1:4],cmdscale_object)
    mds.coords=cbind(complete.pca[subset,1:4],cmdscale_object)
    
    write.csv( mds.coords, mds.file, row.names=FALSE )
    
    png(mdspng)
    plot(mds.coords$MDS1,mds.coords$MDS2)
    dev.off()
    
    png(bychrom_file)
    par(mfrow=c(2,1))
    plot(mds.coords$MDS1,ylab="MDS1")
    plot(mds.coords$MDS2,ylab="MDS2")
    dev.off()
    
    
    
    # figure out where to plot things at
    chroms <- unique(complete.pca$chrom)
    chrom.starts <- tapply( complete.pca$start, complete.pca$chrom, min, na.rm=TRUE )
    chrom.ends <- tapply( complete.pca$end, complete.pca$chrom, max, na.rm=TRUE )
    chrom.spacing <- floor(.05*mean(chrom.ends))
    chrom.offsets <- c(0,cumsum(chrom.spacing+chrom.ends))
    names(chrom.offsets) <- c(chroms,"end")
    chrom.dividers <- c(0,chrom.offsets[-1])-chrom.spacing/2
    chrom.mids <- chrom.dividers[-1] - diff(chrom.dividers)/2
    complete.pca$pos <- chrom.offsets[match(complete.pca$chrom,chroms)]+(complete.pca$start+complete.pca$end)/2
    
    
    
    mds.corners <- corners( mds.coords[,c("MDS1","MDS2")], prop=.05 )
    corner.cols <- brewer.pal(3,"Dark2")
    corner.pch <- c(15,17,19)
    ccols <- rep("black",nrow(mds.coords))
    cpch <- rep(20,nrow(mds.coords))
    for (k in 1:ncol(mds.corners)) {
      ccols[ mds.corners[,k] ] <- corner.cols[k]
      cpch[ mds.corners[,k] ] <- corner.pch[k]
    }
    # centroids of the corners in MDS space
    corner.mds <- do.call(rbind, lapply(1:ncol(mds.corners), 
                                        function (ii){
                                          colMeans(mds.coords[mds.corners[,ii],-(c(1:2,ncol(mds.coords)))])
                                        } ) )
    mds.coords$ccols = ccols
    write.csv( mds.coords, mds.file, row.names=FALSE )
    
    mds.cols <- (1:ncol(mds.coords))[-(c(1:4,ncol(mds.coords)))]
    
    
    layout_heights <- function (k,dl=0,ncol=1) {
      # to set up layout without 'dl' lines between plots
      # use like layout(1:5,heights=layout_heights(5))
      if (k==1) return(1)
      layout(matrix(seq_len(k*ncol),ncol=ncol))  # this changes par("csi")
      ds <- dl*par("lheight")*par("csi")
      eps=par("mai")[c(1,3)]
      dh=(par("din")[2]-sum(eps)-(k-1)*ds)/k
      return(c(eps[2]+dh+ds/2,rep(dh+ds,k-2),eps[1]+dh+ds/2)/par("din")[2])
    }
    chrom.plot <- function (y,ylab='',main='',chrom.labels=TRUE,regions,...) {
      plot(0, type='n', xlim=range(chrom.offsets/1e6), ylim=range(y,finite=TRUE), 
           xlab='', xaxt='n', yaxt='n', ylab=ylab, main=main)
      if (length(chroms)>1) for (k in 1:floor(length(chroms)/2)) {
        rect( xleft=chrom.dividers[2*k-1]/1e6, xright=chrom.dividers[2*k]/1e6, 
              ybottom=par("usr")[3], ytop=par("usr")[4], 
              border=NA, col=adjustcolor("grey",0.25) )
      }
      abline( v=chrom.dividers/1e6, lty=3, col=adjustcolor("grey",0.5), lwd=2 )
      if (chrom.labels) axis( 1, at=chrom.mids/1e6, labels=paste("chromosome", chroms), las=0, tick=FALSE )
      points( regions$pos/1e6, y, ...)
    }
    
    png(finalpng)
    spacing <- 1
    opar <- par(mar=c(4,4,2,1)+.1,mgp=c(2.5,0.8,0))
    layout(matrix(c(rep(1,length(mds.cols)),1+seq_along(mds.cols)),ncol=2),
           widths=c(1,2), heights=layout_heights(length(mds.cols),dl=spacing,ncol=2))
    plot( mds.coords[,mds.cols[1:2]], pch=cpch, 
          col=adjustcolor(ccols,0.75),  asp=1,
          xaxt='n', yaxt='n',
          xlab="MDS coordinate 1", ylab="MDS coordinate 2" )
    points( corner.mds, pch=20, cex=5,
            col=adjustcolor(corner.cols,0.25))
    text( corner.mds, labels=seq_len(nrow(corner.mds)), 
          col=corner.cols, cex=2, lwd=2 )
    opar2 <- par(mar=c(par("mar"),spacing/2)[c(5,2,3,4)])
    for (k in mds.cols) {
      lastone <- (k==mds.cols[length(mds.cols)])
      if (lastone) { par(mar=c(par("mar"),opar2$mar[1])[c(5,2,3,4)]) }
      chrom.plot( y=mds.coords[,k], pch=20,regions=complete.pca[subset,], 
                  xlab=if (lastone) { "Position (Mb)"} else { "" }, # main=paste("MDS coordinate",match(k,mds.cols)),
                  chrom.labels=lastone,
                  ylab=colnames(mds.coords)[k],
                  col=adjustcolor(ccols,0.75) )
      # do this for all but first
      par(mar=c(par("mar"),spacing/2)[c(1,2,5,4)])
    }
    par(opar)
    dev.off()
    
    
    
  } 
}
