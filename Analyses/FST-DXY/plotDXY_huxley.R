# Rscript "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/plotDXY.R"

## for i in *chrfix.txt; do echo $i; head -n 1 $i > $i.temp.txt; cat $i | grep "PseudoNC" | cut -d '_' -f 4 >> $i.temp.txt; done;

listson_number=1

overwrite=F
rawonly=F
plotMeans=F

specieslist = sort(c("bel","bil","bru",
                     "cri","cur","fla","fus",
                     "mel","nit","sin"))
chromlist=sort(c("1","1A","1B","2","3","4","4A","5",
                 "6","7","8","9","10","11","12","13",
                 "14","15","16","17","18","19","20",
                 "21","22","23","24","25","26","27",
                 "28","LG2","LG5","LGE22","mtDNA","Z"))

orig_num_to_keep = 40
windowsize = 100000
movesize = 10000
printstatus = 10

listsonfiles = list.files(path="/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY",
                          pattern="???_*_???_4_Dxy_persite.txt",full.names = T,recursive = T)
#x <- file.info(listsonfiles)
#listsonfiles <- listsonfiles[order(x$size)] ## smallest first

dataset=4

num_to_keep = orig_num_to_keep

for(chrom in chromlist) {
  
  chromfiles = listsonfiles[grep(chrom,listsonfiles)]
  
  if(length(chromfiles)>0){
    for(spp in specieslist) {
      
      sppsonfiles = chromfiles[grep(spp,chromfiles)]
      
      if(length(sppsonfiles)>0){
        
        for(sonfile in sppsonfiles) {
          
          print(sonfile)
          if(file.exists(sonfile)){
            
            print("READING CSV")
            #sonb = read.csv(sonfile,sep="\t",header=T,stringsAsFactors = F)
            sonb=data.table::fread(sonfile,sep="\t",header=T,stringsAsFactors=F,
                                   colClasses=c("character","numeric","numeric"),
                                   showProgress=T)
            
            chrom=strsplit(basename(sonfile),"_")[[1]][2]
            dataset=strsplit(basename(sonfile),"_")[[1]][4]
            
            sonb$position = as.numeric(as.character(sonb$position))
            sonb$dxy = as.numeric(as.character(sonb$dxy))
            sonb = sonb[complete.cases(sonb),]
            
            
            if(nrow(sonb) <=0){
              print("NO DATA")
              next
            } else {
              
              ## how many rows are just zero? 
              justzero=sum(sonb$dxy==0)
              percentzero=justzero/nrow(sonb)
              print(paste("Proportion that is zero:",round(percentzero,2)))
              
              only_zero_data = sonb[sonb$dxy==0,]
              only_full_data = sonb[sonb$dxy!=0,]
              
              print("READING IN SCAFFOLDS")
              num_on_scaff = rev(sort(table(sonb$chromo)))
              totalscaffs = length(num_on_scaff)
              print(paste("TOTAL SCAFFS READ:",totalscaffs))
              
              if (totalscaffs != num_to_keep) {
                num_to_keep=totalscaffs
              }
              
              outfile = paste("/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY/",spp,"_",dataset,"_SON_Dxy_WINDOWS_chrfix_",num_to_keep,"-",chrom,".txt",sep="")
              outfile_TEMP = paste("/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY/",spp,"_",dataset,"_SON_Dxy_WINDOWS_chrfix_",num_to_keep,"-",chrom,".temp",sep="")
              pngfile = paste("/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY/",spp,"_",dataset,"_SON_DxyMEAN_WINDOWS_chrfix_",num_to_keep,"-",chrom,".png",sep="")
              sumfile = paste("/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY/",spp,"_",dataset,"_SON_DxySUMS_WINDOWS_chrfix_",num_to_keep,"-",chrom,".png",sep="")
              sdvfile = paste("/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY/",spp,"_",dataset,"_SON_DxySDVS_WINDOWS_chrfix_",num_to_keep,"-",chrom,".png",sep="")
              rawfile = paste("/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY/",spp,"_",dataset,"_SON_DxyRAW_RAW_chrfix_",num_to_keep,"-",chrom,".png",sep="")
              
              
              if(file.exists(pngfile) && overwrite==F) {
                print("SKIPPING")
              } else {
                
                if(num_to_keep != 1){
                  print("SUBSETTING SCAFFOLDS")
                  scaff_subset = names(num_on_scaff)[1:num_to_keep]
                  
                  son_z = only_zero_data[only_zero_data$chromo %in% scaff_subset,]
                  son_f = only_full_data[only_full_data$chromo %in% scaff_subset,]
                  
                  matches = match(son_f$chromo,scaff_subset)
                  son_f$chromnum = matches
                  son_f = son_f[order(son_f[,"chromnum"], son_f[,"position"] ),]
                  head(son_f)
                  summary(son_f$dxy)
                  summary(son_f$position)
                  nscaf = (length(scaff_subset))
                  
                } else {
                  nscaf=1
                  son_z = only_zero_data
                  son_f = only_full_data
                  scaff_subset = son_f$chromo[1]
                }
                
                closest_snp = min(min(son_f$position),min(son_z$position))
                furthest_snp = max(max(son_f$position),max(son_z$position))
                
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
                plot(x=as.numeric(son_f$position),y=as.numeric(son_f$dxy),
                     col=as.numeric(as.factor(son_f$chromo)),cex=0.2,
                     main=spp,xlab="Position",ylab="Raw DXY NO ZERO",
                     ylim=c(0,0.5))
                text(x=mean(as.numeric(son_f$position),na.rm=T),
                     y=0.9,labels=as.character(round(as.numeric(mean(son_f$dxy,na.rm=T)),2)))
                dev.off()
                
                if(rawonly==F){
                  
                  print("BEGINNING SCAFFOLDS")
                  #for (i in 1:100) {
                  for (i in 1:length(scaff_subset)) {
                    print(paste(i," of ",nscaf,spp))
                    
                    scaf = (scaff_subset[i])
                    number = nums[i]
                    
                    if(nscaf==1){
                      this_scaf_f = son_f
                      this_scaf_z = son_z
                    } else {
                      this_scaf_f = son_f[son_f$chromo==scaf,]
                      this_scaf_z = son_z[son_z$chromo==scaf,]
                    }
                    
                    furthest = max(max(this_scaf_f$position),max(this_scaf_z$position))
                    
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
                      
                      #this_window = this_scaf_f[this_scaf$position>=start,]
                      #this_window = this_window[this_window$position<=end,]
                      
                      nonzero = this_scaf_f[this_scaf_f$position>=start & this_scaf_f$position<=end,]
                      numzero = sum(this_scaf_z$position>=start & this_scaf_z$position<=end)
                      
                      num_snps = numzero + nrow(nonzero)
                      sum_dxy = sum(nonzero$dxy,na.rm=T)
                      #average_dxy = mean(this_window$dxy,na.rm=T)
                      average_dxy = sum_dxy/num_snps
                      #deviation_dxy = sd(this_window$dxy,na.rm=T)
                      deviation_dxy = sqrt((sum((c(nonzero$dxy,rep(0,numzero))-average_dxy)^2))/num_snps)
                      
                      # scaf_list = c(scaf_list,scaf)
                      # start_list = c(start_list,start)
                      # mean_list = c(mean_list,average_dxy)
                      # dev_list = c(dev_list,deviation_dxy)
                      # snp_list = c(snp_list,num_snps)
                      # sum_list = c(sum_list,sum_dxy)
                      
                      
                      
                      if(j==1 && i==1){
                        write.table(cbind("scafs","starts","means","stdvs","snps","sums"),outfile_TEMP,
                                    row.names = F,col.names = F,append=T)
                      }
                      
                      write.table(cbind(scaf,start,average_dxy,deviation_dxy,num_snps,sum_dxy),
                                  outfile_TEMP,
                                  append=T,row.names = F,col.names = F)
                      
                      # output_dataframe = cbind(scaf_list,
                      #                          as.numeric(start_list),
                      #                          as.numeric(mean_list),
                      #                          as.numeric(dev_list),
                      #                          as.numeric(snp_list),
                      #                          as.numeric(sum_list))
                      # colnames(output_dataframe) = c("scafs","starts","means","stdvs","snps","sums")
                    }
                  }
                  #head(output_dataframe)
                  output_dataframe = read.table(outfile_TEMP,header = T,stringsAsFactors = F)
                  output_dataframe=unique(output_dataframe)
                  
                  if(plotMeans==T){
                    
                    
                    png(pngfile,width=700,height=700)
                    par(mfrow=c(2,1))
                    #plot(output_dataframe[,"starts"],output_dataframe[,"means"])
                    palette(c("red","orange","goldenrod","green","blue","purple",
                              "cyan","black","brown","magenta"))
                    plot(as.numeric(output_dataframe$means),
                         col=as.numeric(as.factor(output_dataframe$scafs)),cex=0.2,
                         main=spp,xlab="Window (Scaffold)",ylab="Mean DXY",
                         ylim=c(0,0.5))
                    plot(as.numeric(output_dataframe$means),
                         col=as.numeric(as.factor(output_dataframe$scafs)),cex=0.2,
                         main=spp,xlab="Window (Scaffold)",ylab="Mean DXY")
                    dev.off()
                    
                    png(sumfile,width=700,height=350)
                    #plot(output_dataframe[,"starts"],output_dataframe[,"means"])
                    palette(c("red","orange","goldenrod","green","blue","purple",
                              "cyan","black","brown","magenta"))
                    plot(as.numeric(output_dataframe$sums),
                         col=as.numeric(as.factor(output_dataframe$scafs)),cex=0.2,
                         main=spp,xlab="Window (Scaffold)",ylab="Sum DXY")
                    dev.off()
                    
                    #write.table(output_dataframe,outfile)
                    
                    png(sdvfile,width=700,height=350)
                    #plot(output_dataframe[,"starts"],output_dataframe[,"means"])
                    palette(c("red","orange","goldenrod","green","blue","purple",
                              "cyan","black","brown","magenta"))
                    plot(as.numeric(output_dataframe$stdvs),
                         col=as.numeric(as.factor(output_dataframe$scafs)),cex=0.2,
                         main=spp,xlab="Window (Scaffold)",ylab="STDV DXY")
                    dev.off()
                    
                  }
                  write.table(output_dataframe,outfile,append=T)
                  
                }
              }
            }
            
            
            
            file.rename(sonfile,
                        paste("/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/DXY/DONE/",basename(sonfile),sep=""))
            
          } else {print("file does not exist, skipping")}
        }
      }
    }
  }
}