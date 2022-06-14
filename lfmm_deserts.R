library(lfmm)
library(LEA)
library('raster')
library('RStoolbox') ## not working on huxley
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("LEA")
#devtools::install_github("https://github.com/bcm-uga/lfmm")

#https://bcm-uga.github.io/lfmm/articles/lfmm

setwd("/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/CURVIROSTRE/GENOME/")
#setwd("~/Dropbox (AMNH)/")

reloadEnv = F

## make a pca of the rasters
if(file.exists("worldclim_pca.tif")) {
  pca1 = stack("worldclim_pca.tif")
} else {
  
  if (reloadEnv == T) {
    Env = raster::stack(#'/Users/kprovost/Dropbox (AMNH)/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/bio1.bil'
      list.files(
        path = '/Users/kprovost/Dropbox (AMNH)/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/',
        pattern = "\\.bil$",
        full.names = T
      )#[1]
    )
    ext = raster::extent(c(-117, -98, 27, 35))
    Env = raster::crop(Env, ext)
    bg = Env[[1]] ## just for plotting
    
    border = rgdal::readOGR(
      "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/mexstates/mexico_and_states_subset.shp"
    )
    #plot(border)
    
  } else {
    Env <- getData("worldclim",var="bio",res=10)
    Env = crop(Env,extent(-130,-70,20,50))
    
  }
  
  
  
  
pca1 <- RStoolbox::rasterPCA(Env,filename="worldclim_pca.tif",format="GTiff")
plot(pca1$model$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained") ## k=3?
write.table(summary(pca1$model)$loadings,"worldclim_pca_loadings.txt",sep="\t")
plot(pca1$map)
plot(pca1$map,1)
library(ggplot2)
RStoolbox::ggRGB(pca1$map,1,2,3, stretch="lin", q=0)
}

latlongs = read.csv("~/locations_genome_dissertation.csv",header=T)
#latlongs = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/locations_genome_dissertation.csv",header=T)
latlongs_extract = extract(pca1,latlongs[,c("LONG","LAT")])
latlongs = cbind(latlongs,latlongs_extract)

## first convert vcf to lfmm
vcffile="Toxostoma-curvirostre-called.geno.PseudoNC.all.fixedchroms.converted.vcf"
lfmmfile=paste(tools::file_path_sans_ext(vcffile),".lfmm",sep="")
projectfile=sub("lfmm","snmfProject",lfmmfile)
if(!file.exists(lfmmfile)){
fullfile = readr::read_file(vcffile)
fullfile=sub(" ","",fullfile)
readr::write_file(fullfile,vcffile)
## this will crash R if the header contains spaces in the names, or any other kind of text
## try reading in the first few lines and seeing if there are spaces?
lfmmfile = LEA::vcf2lfmm(vcffile)
}
genofile = sub("lfmm","geno",lfmmfile)

vcfspecies = stringr::str_to_upper(strsplit(vcffile,"-",fixed=F)[[1]][2])
envfile=paste(vcfspecies,"_env_pc1_lfmm.txt",sep="") ## want to only do this once per species
if(file.exists(envfile)){
  envalues = read.table(envfile,sep="\t",header=F)
} else {
envalues = latlongs[latlongs$SPECIES==vcfspecies,colnames(latlongs_extract)]
envalues = envalues[,1] ## only keep pc1
## get the environments for each individual  
#thesevalues = extract(pca1,latlongs_species)
write.table(envalues,envfile,row.names = F,sep="\t",col.names =F)
}

# main options
# K = number of ancestral populations
# entropy = TRUE: computes the cross-entropy criterion,
# CPU = 4 the number of CPUs.
project = NULL
if(file.exists(projectfile)){ 
  project = load.snmfProject(projectfile) 
} else {
project = snmf(genofile, ## takes a while
               K = 1:3,
               entropy = TRUE,
               repetitions = 10,
               project = "new")
}
par(mfrow=c(1,1))
plot(project, col = "blue", pch = 19, cex = 1.2)
## k=2 is  what we will go with for curvirostre? sometimes k=1, sometimes  k=2
best_1 = which.min(cross.entropy(project, K = 1))
best_2 = which.min(cross.entropy(project, K = 2))
best_3 = which.min(cross.entropy(project, K = 3))
my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold")
png(paste(projectfile,"_K2-3_ancestry.png",sep=""))
par(mfrow=c(2,1))
barchart(project, K = 2, run = best_2,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp_2
axis(1, at = 1:length(bp_2$order),
     labels = bp_2$order, las=1,
     cex.axis = .4)
barchart(project, K = 3, run = best_3,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp_3
axis(1, at = 1:length(bp_3$order),
     labels = bp_3$order, las=1,
     cex.axis = .4)
dev.off()

# Population differentiation tests
if(file.exists(paste(projectfile,"_K2_pvalues.txt",sep=""))){
  pvalues = read.table(paste(projectfile,"_K2_pvalues.txt",sep=""))
} else {
p = snmf.pvalues(project, ## also takes a while because goes through every snp
                 entropy = TRUE,
                 ploidy = 2,
                 K = 2)
pvalues = p$pvalues
write.table(pvalues,paste(projectfile,"_K2_pvalues.txt",sep=""))
}
if(!file.exists(paste(projectfile,"_K2_pvalues.png",sep=""))){
png(paste(projectfile,"_K2_pvalues.png",sep=""))
par(mfrow = c(2,1))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)
dev.off()
}

if(file.exists(paste(projectfile,"_K3_pvalues.txt",sep=""))){
  pvalues3 = read.table(paste(projectfile,"_K3_pvalues.txt",sep=""))
} else {
  p3 = snmf.pvalues(project, ## also takes a while because goes through every snp
                   entropy = TRUE,
                   ploidy = 2,
                   K = 3)
  pvalues3 = p3$pvalues
  write.table(pvalues3,paste(projectfile,"_K3_pvalues.txt",sep=""))
}
if(!file.exists(paste(projectfile,"_K3_pvalues.png",sep=""))){
  png(paste(projectfile,"_K3_pvalues.png",sep=""))
  par(mfrow = c(2,1))
  hist(pvalues3, col = "orange")
  plot(-log10(pvalues3), pch = 19, col = "blue", cex = .5)
  dev.off()
}

## TODO: implement this imputed for K1, K2, K3? right now overwrites with whichever is assumed

## make this an if statement as well 
imputelfmm = paste(lfmmfile,"_imputed.lfmm",sep="")
if(!file.exists(imputelfmm)){ 
  project.missing = project
  impute(project.missing, lfmmfile,
         method = 'mode', K = 2, run = best_2)
}
if(!file.exists(paste(lfmmfile,"_imputed.geno",sep=""))){ 
  imputegeno = lfmm2geno(imputelfmm)
}

# project.impute = snmf(imputelfmm, K = 1:3,
#                        entropy = TRUE, repetitions = 10, ## maybe do 100
#                        project = "new")
# plot(project.impute, col = "red", pch = 19, cex = 1.2) ## still k=2 -- CIRCULAR

project_env = NULL
project_env = lfmm(lfmmfile, ## takes a while
               envfile,
               K = 2,
               repetitions = 5,
               project = "new")

project.impute_env = NULL
project.impute_env = lfmm(imputelfmm, ## takes a while 
               envfile,
               K = 2,
               repetitions = 5,
               project = "new")

## this needs to be imputed to work i think 
Y = data.table::fread(imputelfmm)
X = data.table::fread(envfile)
mod.lfmm <- lfmm_ridge(Y = Y, 
                       X = X, 
                       K = 2)  # K = 2 genetic clusters
pp <- lfmm_test(Y = Y, 
                X = X, 
                lfmm = mod.lfmm, 
                calibrate = "gif" ## median+MAD, gif, or NULL
                )
pvalues <- pp$calibrated.pvalue
pvalues

Y2 = data.table::fread(lfmmfile)
mod.lfmm2 <- lfmm_ridge(Y = Y2, 
                       X = X, 
                       K = 2)  # K = 2 genetic clusters
pp2 <- lfmm_test(Y = Y2, 
                X = X, 
                lfmm = mod.lfmm2, 
                calibrate = "gif" ## median+MAD, gif, or NULL
)
pvalues2 <- pp2$calibrated.pvalue
pvalues2
par(mfrow=c(1,1))
plot(pvalues,pvalues2)
abline(a=0,b=1)
abline(h=0.05,col="red")
abline(v=0.05,col="red")

par(mfrow=c(1,2))
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)
qqplot(rexp(length(pvalues2), rate = log(10)), ## hm, this looks better
       -log10(pvalues2), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)


## crashes
pc = pca(lfmmfile, scale = TRUE)

# Perfom Tracy-Widom tests on all eigenvalues.
tw = tracy.widom(pc)

# display p-values for the Tracy-Widom tests (first 5 pcs).
tw$pvalues[1:5]

# plot the percentage of variance explained by each component
plot(tw$percentage)



project.missing = snmf(lfmmfile, K = 4,
                       entropy = TRUE, repetitions = 10,
                       project = "new")


{
#### LEA EXAMPLES
# creation of a directory for LEA analyses
dir.create("LEA_analyses")
# set the created directory as the working directory
setwd("LEA_analyses")

## TODO: look up how to convert this data

library(LEA)
# Creation a the genotypic file: "genotypes.lfmm"
# The data include 400 SNPs for 50 individuals.
data("tutorial")
# Write genotypes in the lfmm format
write.lfmm(tutorial.R, "genotypes.lfmm")
# Write genotypes in the geno format
write.geno(tutorial.R, "genotypes.geno")
# creation of an environment gradient file: gradient.env.
# The .env file contains a single ecological variable
# for each individual.
write.env(tutorial.C, "gradients.env")
}
# run of pca
# Available options, K (the number of PCs),
# center and scale.
# Create files: genotypes.eigenvalues - eigenvalues,
# genotypes.eigenvectors - eigenvectors,
# genotypes.sdev - standard deviations,
# genotypes.projections - projections,
# Create a pcaProject object: pc.
pc = pca("genotypes.lfmm", scale = TRUE)

# Perfom Tracy-Widom tests on all eigenvalues.
# create file: tuto.tracyWidom - tracy-widom test information.
tw = tracy.widom(pc)

# display p-values for the Tracy-Widom tests (first 5 pcs).
tw$pvalues[1:5]

# plot the percentage of variance explained by each component
plot(tw$percentage)
{
# main options
# K = number of ancestral populations
# entropy = TRUE: computes the cross-entropy criterion,
# CPU = 4 the number of CPUs.
project = NULL
project = snmf("genotypes.geno", ## takes a while
               K = 1:10,
               entropy = TRUE,
               repetitions = 10,
               project = "new")

# plot cross-entropy criterion for all runs in the snmf project
plot(project, col = "blue", pch = 19, cex = 1.2)

# select the best run for K = 4
best = which.min(cross.entropy(project, K = 4))
my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold")
barchart(project, K = 4, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4)

# Population differentiation tests
p = snmf.pvalues(project,
                 entropy = TRUE,
                 ploidy = 2,
                 K = 4)
pvalues = p$pvalues
par(mfrow = c(2,1))
hist(pvalues, col = "orange")
plot(-log10(pvalues), pch = 19, col = "blue", cex = .5)


# creation of a genotypic matrix with missing genotypes
dat = as.numeric(tutorial.R)
dat[sample(1:length(dat), 100)] <- 9 ## 9 is missing data 
dat <- matrix(dat, nrow = 50, ncol = 400)
write.lfmm(dat, "genoM.lfmm")

project.missing = snmf("genoM.lfmm", K = 4,
                       entropy = TRUE, repetitions = 10,
                       project = "new")

# select the run with the lowest cross-entropy value
best = which.min(cross.entropy(project.missing, K = 4))
# Impute the missing genotypes
impute(project.missing, "genoM.lfmm",
       method = 'mode', K = 4, run = best)
## Missing genotype imputation for K = 4
## Missing genotype imputation for run = 3
## Results are written in the file: genoM.lfmm_imputed.lfmm
# Proportion of correct imputation results
dat.imp = read.lfmm("genoM.lfmm_imputed.lfmm")
mean( tutorial.R[dat == 9] == dat.imp[dat == 9] )
}

# main options:
# K the number of latent factors
# Runs with K = 6 and 5 repetitions.
project = NULL
project = lfmm("genotypes.lfmm",
               "gradients.env",
               K = 6,
               repetitions = 5,
               project = "new")
## LFMM uses a very naive imputation method which has low power when
## genotypes are missing: See impute() for a better imputation method.

# compute adjusted p-values
p = lfmm.pvalues(project, K = 6)
pvalues = p$pvalues

for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("Expected FDR:", alpha))
  L = length(pvalues)
  # return a list of candidates with expected FDR alpha.
  # Benjamini-Hochberg's algorithm:
  w = which(sort(pvalues) < alpha * (1:L) / L)
  candidates = order(pvalues)[w]
  # estimated FDR and True Positive Rate
  Lc = length(candidates)
  estimated.FDR = sum(candidates <= 350)/Lc
  print(paste("Observed FDR:",
              round(estimated.FDR, digits = 2)))
  
  # GWAS significance test
  par(mfrow = c(2,1))
  hist(pvalues, col = "lightblue")
  plot(-log10(pvalues), pch = 19, col = "blue", cex = .7)
  
  estimated.TPR = sum(candidates > 350)/50
  print(paste("Estimated TPR:",
              round(estimated.TPR, digits = 2)))
}

# Simulate non-null effect sizes for 10 target loci
#individuals
n = 100
#loci
L = 1000
# Environmental variable
X = as.matrix(rnorm(n))
# effect sizes
B = rep(0, L)
target = sample(1:L, 10)
B[target] = runif(10, -10, 10)

# Create 3 hidden factors and their loadings
U = t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X)) +
  matrix(rnorm(3*n), ncol = 3)
V <- matrix(rnorm(3*L), ncol = 3)

# Simulate a matrix containing haploid genotypes
Y <- tcrossprod(as.matrix(X), B) +
  tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)

# Fitting an LFMM with K = 3 factors
mod <- lfmm2(input = Y, env = X, K = 3)

# Computing P-values and plotting their minus log10 values
pv <- lfmm2.test(object = mod,
                 input = Y,
                 env = X,
                 linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .6, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")



###### EXAMPLES


##### STARTING

## Simulated phenotypes for Arabidopsis thaliana SNP data
data("example.data")
## Simulated (and real) methylation levels for sun-exposed tissue sampled
data("skin.exposure")

Y <- example.data$genotype
pc <- prcomp(Y)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(6,pc$sdev[6]^2, type = "h", lwd = 3, col = "blue")

Y <- skin.exposure$beta.value
pc <- prcomp(Y)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(2,pc$sdev[2]^2, type = "h", lwd = 3, col = "blue")

Y <- example.data$genotype
X <- example.data$phenotype #scaled phenotype

## Fit an LFMM, i.e, compute B, U, V estimates
mod.lfmm <- lfmm_ridge(Y = Y, 
                       X = X, 
                       K = 6)

## performs association testing using the fitted model:
pv <- lfmm_test(Y = Y, 
                X = X, 
                lfmm = mod.lfmm, 
                calibrate = "gif")

pvalues <- pv$calibrated.pvalue 
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)

## Manhattan plot
plot(-log10(pvalues), 
     pch = 19, 
     cex = .2, 
     xlab = "SNP", ylab = "-Log P",
     col = "grey")
points(example.data$causal.set, 
       -log10(pvalues)[example.data$causal.set], 
       type = "h", 
       col = "blue")

##### RIDGE ESTIMATES AND EWAS TEST

Y <- scale(skin.exposure$beta.value)
X <- scale(as.numeric(skin.exposure$exposure))

## Fit and LFMM, i.e, compute B, U, V estimates
mod.lfmm <- lfmm_ridge(Y = Y, 
                       X = X, 
                       K = 2)

## Perform association testing using the fitted model:
pv <- lfmm_test(Y = Y, 
                X = X, 
                lfmm = mod.lfmm, 
                calibrate = "gif")

## Manhattan plot
plot(-log10(pv$calibrated.pvalue), 
     pch = 19, 
     cex = .3,
     xlab = "Probe", ylab = "-Log P",
     col = "grey")
causal.set <- seq(11, 1496, by = 80)
points(causal.set, 
       -log10(pv$calibrated.pvalue)[causal.set], 
       col = "blue")


## LASSO ESTIMATES AND GWAS 

Y <- example.data$genotype
X <- example.data$phenotype #scaled phenotype

## Fit an LFMM, i.e, compute B, U, V estimates
mod.lfmm <- lfmm_lasso(Y = Y, 
                       X = X, 
                       K = 6,
                       nozero.prop = 0.01)

## performs association testing using the fitted model:
pv <- lfmm_test(Y = Y, 
                X = X, 
                lfmm = mod.lfmm, 
                calibrate = "gif")

pvalues <- pv$calibrated.pvalue 
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)

## Manhattan plot
plot(-log10(pvalues), 
     pch = 19, 
     cex = .2, 
     xlab = "SNP", ylab = "-Log P",
     col = "grey")
points(example.data$causal.set, 
       -log10(pvalues)[example.data$causal.set], 
       type = "h", 
       col = "blue")


#### PREDICTION OF PHENOTYPES 

## Simulation of 1000 genotypes for 100 individuals (y)
u <- matrix(rnorm(300, sd = 1), nrow = 100, ncol = 2) 
v <- matrix(rnorm(3000, sd = 2), nrow = 2, ncol = 1000)
y <- matrix(rbinom(100000, size = 2, 
                   prob = 1/(1 + exp(-0.3*(u%*%v 
                                           + rnorm(100000, sd = 2))))),
            nrow = 100,
            ncol = 1000)

## Simulation of 1000 phenotypes (x)
## Only the last 10 genotypes have significant effect sizes (b)
b <- matrix(c(rep(0, 990), rep(6000, 10)))
x <- y%*%b + rnorm(100, sd = 100)

mod <- lfmm_ridge(Y = y, 
                  X = x,
                  K = 2)

candidates <- 991:1000 #causal loci
b.values <- effect_size(Y = y, X = x, lfmm.object = mod) 
x.pred <- scale(y[,candidates], scale = F)%*% matrix(b.values[candidates])


##Compare simulated and predicted/fitted phenotypes
plot(x - mean(x), x.pred, 
     pch = 19, col = "grey", 
     xlab = "Observed phenotypes (centered)", 
     ylab = "Predicted from PRS")
abline(0,1)
abline(lm(x.pred ~ scale(x, scale = FALSE)), col = 2)


pred <- predict_lfmm(Y = y, 
                     X = x,
                     fdr.level = 0.25, 
                     mod)

##Compare simulated and predicted/fitted phenotypes
plot(x - mean(x), pred$pred, 
     pch = 19, col = "grey", 
     xlab = "Observed phenotypes (centered)", 
     ylab = "Predicted from PRS")
abline(0,1)
abline(lm(pred$pred ~ scale(x, scale = FALSE)), col = 2)



# http://membres-timc.imag.fr/Olivier.Francois/LEA/files/LEA_github.pdf


