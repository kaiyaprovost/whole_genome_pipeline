# Rscript "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/plotDXY.R"

specieslist = c(#"bel",
                #"bil","bru",
                #"cri",
                #"cur","fla",
                #"fus",
                #"mel",
                "sin",
                "nit"
                )
                
#specieslist = c("cri")
                
orig_num_to_keep = 100
windowsize = 100000
movesize = 10000

for (spp in specieslist) {
  
  print(spp)
  
  num_to_keep = orig_num_to_keep
  
  sonfile = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/raw/",spp,"_SON_Dxy_persite_chrfix.txt",sep="")

  print("READING CSV")
  sonb = read.csv(sonfile,sep="\t",header=T)
  
  print("READING IN SCAFFOLDS")
  num_on_scaff = rev(sort(table(sonb$chromo)))
  totalscaffs = length(num_on_scaff)
  print(paste("TOTAL SCAFFS READ:",totalscaffs))
  
  if (totalscaffs < num_to_keep) {
    num_to_keep=totalscaffs
  }
  
  outfile = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/",spp,"_SON_Dxy_WINDOWS_chrfix_",num_to_keep,".txt",sep="")
  pngfile = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/",spp,"_SON_DxyMEAN_WINDOWS_chrfix_",num_to_keep,".png",sep="")
  sumfile = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/",spp,"_SON_DxySUMS_WINDOWS_chrfix_",num_to_keep,".png",sep="")
  sdvfile = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/",spp,"_SON_DxySDVS_WINDOWS_chrfix_",num_to_keep,".png",sep="")
  
  
  print("SUBSETTING SCAFFOLDS")
  scaff_subset = names(num_on_scaff)[1:num_to_keep]
  
  son = sonb[sonb$chromo %in% scaff_subset,]
  matches = match(son$chromo,scaff_subset)
  son$chromnum = matches
  son = son[order(son[,"chromnum"], son[,"position"] ),]
  
  head(son)
  summary(son$dxy)
  summary(son$position)
  
  closest_snp = min(son$position)
  furthest_snp = max(son$position)
  
  ## son and chi are exactky the same! yay
  
  
  
  nscaf = (length(scaff_subset))
  
  print(paste("NUMBER OF SCAFFS TO OUTPUT:",nscaf))
  
  nums = as.numeric(num_on_scaff)
  
  scaf_list = c()
  start_list = c()
  mean_list = c()
  dev_list = c()
  snp_list = c()
  sum_list = c()
  
  print("BEGINNING SCAFFOLDS")
  #for (i in 1:100) {
    for (i in 1:length(scaff_subset)) {
    print(paste(i," of ",nscaf,spp))
    
    scaf = (scaff_subset[i])
    number = nums[i]
    this_scaf = son[son$chromo==scaf,]
    #closest = min(this_scaf$position)
    furthest = max(this_scaf$position)
    
    ## get windows
    window_ends = seq(windowsize,max(furthest,windowsize),movesize)
    window_starts = seq(1,furthest,movesize)
    window_starts = window_starts[1:length(window_ends)]
    windows = as.data.frame(cbind("starts"=window_starts,"ends"=window_ends))
    
    #print(nrow(windows))
    
    for (j in 1:nrow(windows)){
      window = windows[j,]
      #print(window)
      start = window$starts
      end = window$ends
      
      this_window = this_scaf[this_scaf$position>=start,]
      this_window = this_window[this_window$position<=end,]
      
      num_snps = nrow(this_window)
      average_dxy = mean(this_window$dxy)
      deviation_dxy = sd(this_window$dxy)
      sum_dxy = sum(this_window$dxy)
      
      scaf_list = c(scaf_list,scaf)
      start_list = c(start_list,start)
      mean_list = c(mean_list,average_dxy)
      dev_list = c(dev_list,deviation_dxy)
      snp_list = c(snp_list,num_snps)
      sum_list = c(sum_list,sum_dxy)
      
    }
    
    output_dataframe = cbind(scaf_list,
                             as.numeric(start_list),
                             as.numeric(mean_list),
                             as.numeric(dev_list),
                             as.numeric(snp_list),
                             as.numeric(sum_list))
    colnames(output_dataframe) = c("scafs","starts","means","stdvs","snps","sums")
    
  }
  head(output_dataframe)
  png(pngfile,width=700,height=350)
  #plot(output_dataframe[,"starts"],output_dataframe[,"means"])
  palette(c("red","orange","goldenrod","green","blue","purple",
            "cyan","black","brown","magenta"))
  plot(as.numeric(output_dataframe[,"means"]),col=as.numeric(as.factor(output_dataframe[,"scafs"])),cex=0.2,
       main=spp,xlab="Window (Scaffold)",ylab="Mean DXY",
       ylim=c(0,1))
  dev.off()
  
  png(sumfile,width=700,height=350)
  #plot(output_dataframe[,"starts"],output_dataframe[,"means"])
  palette(c("red","orange","goldenrod","green","blue","purple",
            "cyan","black","brown","magenta"))
  plot(as.numeric(output_dataframe[,"sums"]),col=as.numeric(as.factor(output_dataframe[,"scafs"])),cex=0.2,
       main=spp,xlab="Window (Scaffold)",ylab="Sum DXY")
  dev.off()
  
  write.table(output_dataframe,outfile)
  
  png(sdvfile,width=700,height=350)
  #plot(output_dataframe[,"starts"],output_dataframe[,"means"])
  palette(c("red","orange","goldenrod","green","blue","purple",
            "cyan","black","brown","magenta"))
  plot(as.numeric(output_dataframe[,"stdvs"]),col=as.numeric(as.factor(output_dataframe[,"scafs"])),cex=0.2,
       main=spp,xlab="Window (Scaffold)",ylab="STDV DXY")
  dev.off()
  
  #write.table(output_dataframe,outfile)
  
}
