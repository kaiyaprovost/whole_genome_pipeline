files=list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/",
                 pattern="imiss.imiss",full.names = T)

for(file in files){
  
  df = read.table(file,header=T,sep="\t")
  print(basename(file))
  print(mean(df$F_MISS,na.rm=T))
  
}

## generate window files for lmiss
## every 10000 windows of 100000 starting at 0-50000-100000

file = "/Users/kprovost/Downloads/MISSING/Vireo-bellii-called.geno.vcf.lmiss.lmiss"
df = data.table::fread(file,data.table=F)




"Error in colMeans(mattemp[temp$POS >= starts[i] & temp$POS < stops[i],  : 
  'x' must be an array of at least two dimensions
Called from: colMeans(mattemp[temp$POS >= starts[i] & temp$POS < stops[i], 
    ], na.rm = T)
Error during wrapup: unimplemented type (29) in 'eval'

Error: no more error handlers available (recursive errors?); invoking 'abort' restart
Error during wrapup: INTEGER() can only be applied to a 'integer', not a 'unknown type #29'
Error: no more error handlers available (recursive errors?); invoking 'abort' restart"
## got through mtdna, 1-3, 1a 1b, then failed on 4
full= lapply(unique(df$CHR),FUN=function(chrom) {
  print(chrom)
  temp = df[df$CHR==chrom,]
  end = max(temp$POS,na.rm=T)
  if(round(end,-4)<end){
    end=round(end,-4)+10000
  } else {
    end=round(end,-4)
  }
  if(end < 100000) {
    starts=0
    stops=100000
  } else {
    starts = seq(0,end,10000)
    stops=starts+100000
  }
  print(paste("LENGTH TO DO:",length(starts)))
  mattemp = as.matrix(temp[,c(3,5,6)])
  
  out = lapply(1:length(starts),FUN=function(i){
    a1=colMeans(mattemp[temp$POS>=starts[i] & temp$POS<stops[i],],na.rm=T)
    a2=colSums(mattemp[temp$POS>=starts[i] & temp$POS<stops[i],],na.rm=T)
    a3=matrixStats::colSds(mattemp[temp$POS>=starts[i] & temp$POS<stops[i],],na.rm=T)
    names(a1) = paste("mean",colnames(mattemp),sep="_")
    names(a2) = paste("sum",colnames(mattemp),sep="_")
    names(a3) = paste("sd",colnames(mattemp),sep="_")
    a4=matrix(cbind(a1,a2,a3),nrow=1)
    colnames(a4) = c(names(a1),names(a2),names(a3))
    a4 = as.data.frame(a4)
    a4$sum_N_DATA_div_sum_N_MISS = a4$sum_N_MISS / a4$sum_N_DATA
    a4$start = starts[i]
    a4$stop = stops[i]
    a4$midPos = (a4$stop+a4$start)/2
    a4$chr= chrom
    return(a4)
  })
  outdf = do.call(rbind,out)
  return(outdf)
})
