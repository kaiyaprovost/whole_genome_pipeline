setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/Distances/")

files = c(
  #"Amphispiza-bilineata_ngsdist",
  #"Melozone-fusca_ngsdist",
  #"Polioptila-melanura_ngsdist",
  #"Campylorhynchus-brunneicapillus_ngsdist",
  #"Cardinalis-sinuatus_ngsdist",
  #"Toxostoma-crissale_ngsdist",
  "Toxostoma-curvirostre_ngsDist"
  #"Vireo-bellii-NOWEIRD_ndsDist",
  #"Vireo-bellii-NOWEIRD_ndsdist",
  #"Phainopepla-nitens-NOWEIRD_ngsDist",
  #"Auriparus-flaviceps-NOWEIRD_ngsDist",
  #"Phainopepla-nitensDXY_ngsDist"
)

for (file in files) {
  
  #spp = substr(file,71,100)
  spp = strsplit(file,"_")[[1]][1]
  print(spp)
  
  ## the first distance matrix is the computed distance
  ## the other ones are bootstraps
  
  dists = readLines(file)
  dists = dists[2:length(dists)]
  numinds = as.numeric(dists[1])
  reps = numinds+2 
  #firstmatrix = dists[1:reps]
  #bootmatrix = dists[(reps+1):length(dists)]
  
  firstmatrix = read.csv(file,nrow=reps-1,header=F,sep="\t",row.names = 1)
  firstmatrix = firstmatrix[2:nrow(firstmatrix),]
  colnames(firstmatrix) = rownames(firstmatrix)
  
  write.csv(firstmatrix,paste(file,"_FIRSTMATRIX.csv"))
  
  firstdf = as.data.frame(firstmatrix)
  #colnames(firstdf) = NULL
  #rownames(firstdf) = NULL
  
  
  png(paste(file,"distancetree_corrplot_first.png",sep="_"))
  corrplot::corrplot(as.matrix(firstdf)*100,is.corr=F,method="color",order="alphabet",
                     main=spp)
  dev.off()
  
  #pairwisefirst = data.frame( t(combn(names(firstmatrix),2)), dist=t(firstmatrix)[lower.tri(firstmatrix)] )
  #newdist = as.dist(xtabs(pairwisefirst[, 3] ~ pairwisefirst[, 2] + pairwisefirst[, 1]),diag=T,upper=T)
  
  latermatrix = read.csv(file,header=F,skip=reps,sep="\t")
  latermatrix = latermatrix[latermatrix$V1!=numinds,]
  colnames(latermatrix)[2:ncol(latermatrix)] = as.character(unique(latermatrix$V1))
  laterdf = as.data.frame(latermatrix)
  

  agglater = aggregate(latermatrix[,-1],by=list(latermatrix$V1),FUN=function(x) {mean(x,na.rm=T)})
  rownames(agglater) = agglater$Group.1
  agglater = as.data.frame(agglater)
  agglater = agglater[,2:ncol(agglater)]
  #colnames(agglater) = sort(rownames(agglater))
  png(paste(file,"distancetree_corrplot_boots.png",sep="_"))
  
  corrplot::corrplot(as.matrix(agglater)*100,is.corr=F,method="color",order="alphabet",
                     main=spp)
  dev.off()
  
  write.csv(agglater,paste(file,"_BOOTS.csv"))
  
  
  library(ape)
  library(phangorn)
  
  UPGMA <- phangorn::upgma(firstmatrix)
  png(paste(file,"distancetree.png",sep="_"))
  par(mfrow=c(1,2),mar=c(0,1,1,0))
  plot(UPGMA, main=spp,type="unrooted",cex=0.5)
  plot(UPGMA, main=spp,type="phylogram")
  dev.off()
  
  UPGMAboot <- phangorn::upgma(agglater)
  png(paste(file,"distancetree_boots.png",sep="_"))
  par(mfrow=c(1,2),mar=c(0,1,1,0))
  plot(UPGMAboot, main=spp,type="unrooted",cex=0.5)
  plot(UPGMAboot, main=spp,type="phylogram")
  dev.off()
  
  
  
  
}
