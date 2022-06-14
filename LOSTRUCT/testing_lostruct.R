library(lostruct)
library(vcfR)

## bcf.file="Vireo-bellii-called.geno.PseudoNC_007897.1_Tgut_mtDNA.fixedchroms.converted.vcf.gz"

## vcf files must not have any spaces in the indvidual names for this to work 



size=100000
type="bp"
npc=2
weights=1

all.pcas <- numeric(0)       # will be a numeric matrix of eigen values/vectors
all.lengths <- numeric(0)    # will be a numeric vector of numbers of windows per chromosome
all.regions <- data.frame()  # will be a data.frame of the chromsome, start, stop for each window

for(bcf.file in listfiles) {

sites=vcf_positions(bcf.file)
samples=vcf_samples(bcf.file)

win.fn <- vcf_windower(bcf.file, size=size,type=tolower(type),
                            sites=sites,samples=samples)

#win.fn <- vcf_windower_bp_2(bcf.file, size=size,
#                       sites=sites,samples=samples)
these.regions <- region(win.fn)()

pca.stuff <- eigen_windows(data=win.fn, k=npc, w=weights) 

all.pcas <- rbind( all.pcas, pca.stuff )
all.lengths <- c(all.lengths, nrow(pca.stuff))
all.regions <- rbind( all.regions, these.regions )

}


# as in POPRES_SNPdata_recode12.R
#chr22 <- read_tped("POPRES_Genotypes_QC2_v2_TXT.tped.gz", chrom=22)
#templost=read.table("/Users/kprovost/temp.temp.vcf",na.strings=".")
#windows=eigen_windows(data=templost[,-1],k=3,win=10)
#dist=pc_dist(windows)

snps <- read.vcfR("/Users/kprovost/Downloads/cardinalis vcf/cardcard16_filtered.vcf")
#snps2 <- pegas::read.vcf("/Users/kprovost/Downloads/cardinalis vcf/cardcard16_filtered.vcf")

vcf=snps@gt[,-1]
#vcf=data.table::fread("/Users/kprovost/Downloads/cardinalis vcf/cardcard16_filtered_TEST.vcf")
vcf=as.data.frame(vcf)
#vcf=vcf[,10:ncol(vcf)]
val_het=c("0|1","1|0")
val_maj="0|0"
val_min="1|1"

for(colnum in 1:ncol(vcf)){
  print(colnum)
  vcf[is.na(vcf[,colnum]),colnum] = 0
  vcf[vcf[,colnum] %in% c(val_maj,val_min),colnum] = 1
  vcf[vcf[,colnum] %in% val_het,colnum] = 2
  vcf[,colnum] = as.numeric(vcf[,colnum])
}

head(vcf)
write.table(vcf,"/Users/kprovost/Downloads/cardinalis vcf/cardcard16_filtered.vcf.lostruct")
windows=eigen_windows(data=vcf,k=2,win=1000) ## takes long time, longer if window smaller, or if more K
write.table(windows,"/Users/kprovost/Downloads/cardinalis vcf/cardcard16_filtered.vcf.lostruct.windows")
dist=pc_dist(windows,npc=2)
keep=complete.cases(dist[,1])
dist_full=dist[keep,keep]

fit2d=cmdscale(dist_full,eig=TRUE, k=2 )
plot( fit2d$points, xlab="Coordinate 1", 
      ylab="Coordinate 2" #,col=rainbow(1.2*nrow(dist_full)) 
      )


snps <- read_vcf("/Users/kprovost/Downloads/cardinalis vcf/cardcard16_filtered.vcf")
pcs <- eigen_windows(snps,k=2)
pcdist <- pc_dist(pcs,npc=2)