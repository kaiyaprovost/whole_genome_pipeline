# Rscript "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/plotDXY.R"

## for i in *chrfix.txt; do echo $i; head -n 1 $i > $i.temp.txt; cat $i | grep "PseudoNC" | cut -d '_' -f 4 >> $i.temp.txt; done;

specieslist = rev(c("bel",
                    "bil",
                    "bru",
                    "cri",
                    "cur",
                    "fla",
                    "fus",
                    "mel",
                    "sin",
                    "nit"
))
chromlist=c("10","11","12","13","14","15","16","17","18","19",
            "1A","1B","2","20","21","22","23","24","25","26","27","28","3",
            "4","4A","5","6","7","8","9","LG2","LG5","LGE22",
            "mtDNA","Z")

#specieslist = c("cri")

orig_num_to_keep = 40
windowsize = 100000
movesize = 10000
printstatus = 10

listsonfiles = list.files(path="/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY",
                          pattern="???_*_???_4_Dxy_persite",full.names = T)
x <- file.info(listsonfiles)
listsonfiles <- listsonfiles[order(x$size)]

dataset=4

for (spp in (specieslist)) {
  
  print(spp)
  
  num_to_keep = orig_num_to_keep
  
  #sonfile = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/raw/",spp,"_SON_Dxy_persite_chrfix.txt",sep="")
  sppsonfiles = listsonfiles[grep(spp,listsonfiles)]
  
  for(sonfile in sppsonfiles) {
    
    chrom=strsplit(strsplit(sonfile,"-")[[1]][3],"\\.")[[1]][1]
    chrom=strsplit(strsplit(chrom,"/")[[1]][6],"_")[[1]][2]
    #chrom=strsplit(strsplit(sonfile,"-")[[1]][2],"\\.")[[1]][1]
    #dataset=strsplit(sonfile,"_")[[1]][4]
    
    print("READING CSV")
    sonb = read.csv(sonfile,sep="\t",header=T)
    print(chrom)
    if(is.null(chrom)) {
      chrom="full"
    }
    
    if(nrow(sonb) <=0){
      print("NO DATA")
      next
    } else {
      
      
      
      ## how many rows are just zero? 
      justzero=sum(sonb$dxy==0)
      percentzero=justzero/nrow(sonb)
      print(paste("Proportion that is zero:",round(percentzero,2)))
      
      print("READING IN SCAFFOLDS")
      num_on_scaff = rev(sort(table(sonb$chromo)))
      totalscaffs = length(num_on_scaff)
      print(paste("TOTAL SCAFFS READ:",totalscaffs))
      
      if (totalscaffs != num_to_keep) {
        num_to_keep=totalscaffs
      }
      
      outfile = paste("/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY/",spp,"_",dataset,"_SON_Dxy_WINDOWS_chrfix_",num_to_keep,"-",chrom,".txt",sep="")
      pngfile = paste("/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY/",spp,"_",dataset,"_SON_DxyMEAN_WINDOWS_chrfix_",num_to_keep,"-",chrom,".png",sep="")
      sumfile = paste("/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY/",spp,"_",dataset,"_SON_DxySUMS_WINDOWS_chrfix_",num_to_keep,"-",chrom,".png",sep="")
      sdvfile = paste("/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY/",spp,"_",dataset,"_SON_DxySDVS_WINDOWS_chrfix_",num_to_keep,"-",chrom,".png",sep="")
      rawfile = paste("/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY/",spp,"_",dataset,"_SON_DxyRAW_RAW_chrfix_",num_to_keep,"-",chrom,".png",sep="")
      
      
      #if(num_to_keep != 1){
      print("SUBSETTING SCAFFOLDS")
      scaff_subset = names(num_on_scaff)[1:num_to_keep]
      son = sonb[sonb$chromo %in% scaff_subset,]
      matches = match(son$chromo,scaff_subset)
      son$chromnum = matches
      son = son[order(son[,"chromnum"], son[,"position"] ),]
      head(son)
      summary(son$dxy)
      summary(son$position)
      nscaf = (length(scaff_subset))
      #} else {
      #  nscaf=1
      #}
      
      closest_snp = min(son$position)
      furthest_snp = max(son$position)
      
      ## son and chi are exactky the same! yay
      
      print(paste("NUMBER OF SCAFFS TO OUTPUT:",nscaf))
      
      nums = as.numeric(num_on_scaff)
      
      scaf_list = c()
      start_list = c()
      mean_list = c()
      dev_list = c()
      snp_list = c()
      sum_list = c()
      
      print("writing rawfile")
      png(rawfile,width=700,height=350)
      #plot(output_dataframe[,"starts"],output_dataframe[,"means"])
      palette(c("red","orange","goldenrod","green","blue","purple",
                "cyan","black","brown","magenta"))
      plot(x=as.numeric(son[,"position"]),y=as.numeric(son[,"dxy"]),
           col=as.numeric(as.factor(son[,"chromo"])),cex=0.2,
           main=spp,xlab="Position",ylab="Raw DXY",
           ylim=c(0,1))
      text(x=mean(as.numeric(son[,"position"]),na.rm=T),
           y=0.9,labels=as.character(round(as.numeric(mean(son[,"dxy"],na.rm=T)),2)))
      dev.off()
      
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
        
        totwindows=(nrow(windows))
        print(totwindows)
        
        for (j in 1:nrow(windows)){
          window = windows[j,]
          if(j %% printstatus == 0) {
            print(paste(j,"/",totwindows))
          }
          start = window$starts
          end = window$ends
          
          this_window = this_scaf[this_scaf$position>=start,]
          this_window = this_window[this_window$position<=end,]
          
          nonzero = this_window[which(this_window$dxy!=0),]
          numzero = sum(this_window$dxy==0)
          
          num_snps = nrow(this_window)
          sum_dxy = sum(nonzero$dxy,na.rm=T)
          #average_dxy = mean(this_window$dxy,na.rm=T)
          average_dxy = sum_dxy/num_snps
          #deviation_dxy = sd(this_window$dxy,na.rm=T)
          deviation_dxy = sqrt((sum((this_window$dxy-average_dxy)^2))/num_snps)
          
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
      png(pngfile,width=700,height=700)
      par(mfrow=c(2,1))
      #plot(output_dataframe[,"starts"],output_dataframe[,"means"])
      palette(c("red","orange","goldenrod","green","blue","purple",
                "cyan","black","brown","magenta"))
      plot(as.numeric(output_dataframe[,"means"]),col=as.numeric(as.factor(output_dataframe[,"scafs"])),cex=0.2,
           main=spp,xlab="Window (Scaffold)",ylab="Mean DXY",
           ylim=c(0,1))
      plot(as.numeric(output_dataframe[,"means"]),col=as.numeric(as.factor(output_dataframe[,"scafs"])),cex=0.2,
           main=spp,xlab="Window (Scaffold)",ylab="Mean DXY")
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
      
      write.table(output_dataframe,outfile)
    }
  }
}
