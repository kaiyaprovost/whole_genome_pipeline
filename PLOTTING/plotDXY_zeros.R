# Rscript "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/plotDXY_zeros.R"

## for i in *chrfix.txt; do echo $i; head -n 1 $i > $i.temp.txt; cat $i | grep "PseudoNC" | cut -d '_' -f 4 >> $i.temp.txt; done;

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("A folder must be input for script to work", call.=FALSE)
} else {
  folder = args[1]
}
# folder = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/"

chromlist=c("mtDNA",1,"1A","1B",2:4,"4A",5:28,"Z","LG2","LG5","LGE22")
chromLengths = c(76753,132631505,98848917,1686487,156780119,134338477,73199573,21100058,
                 65970967,38010788,40056312,26845016,26829339,21327751,22200599,21367209,
                 19355181,16997682,16163187,74801,12219782,12615983,11868285,15549540,8239338,
                 5000089,8114938,8747811,3133669,6990484,6667436,6132206,78297685,204071,3391,1340002)

orig_num_to_keep = 40
windowsize = 100000
movesize = 10000
printstatus = 10


listsonfiles = list.files(path=folder,
                          pattern="Dxy_persite",full.names = T)
listsonfiles=listsonfiles[!(grepl("windows.dxy.txt",listsonfiles))]
#x <- file.info(listsonfiles)
#listsonfiles <- listsonfiles[order(x$size)]


for(sonfile in listsonfiles) {
  num_to_keep = orig_num_to_keep
  outfile = paste(sonfile,".windows.dxy.txt",sep="")
  
  
  if(file.exists(outfile)){
    print(sonfile)
    print("SKIPPING")
  } else {
    print("READING CSV")
    
    sonb=data.table::fread(sonfile,sep="\t",
                           colClasses=c("character","numeric","numeric"),
                           showProgress=T,data.table = F)
    
    ## check if the names match?
    overlapchr=intersect(chromlist,unique(sonb$chromo))
    
    if(length(overlapchr)==0) {
      ## get rid of any that start with NW, PseudoNW, or SS
      prefixes=substr(sonb$chromo,1,2)
      sonb=sonb[prefixes!="NW",]
      sonb=sonb[prefixes!="SS",]
      sonb=sonb[!(is.na(prefixes)),]
      prefixes=substr(sonb$chromo,1,8)
      sonb=sonb[prefixes!="PseudoNW",]
      
      ## now get rid of the prefix altogether
      oldchromo = unique(sonb$chromo)
      newchromo=sapply(oldchromo,FUN=function(x){
        y=strsplit(x,"_")[[1]]
        z=y[length(y)]
        return(z)
      })
      for(i in 1:length(newchromo)){
        print(i)
        oldchr = oldchromo[i]
        newchr = newchromo[i]
        sonb$chromo[sonb$chromo==oldchr] = as.character(newchr)
      }
      write.table(sonb,paste(sonfile,"_RFIXED.txt",sep=""),
                  row.names = F,col.names = T,append=F,quote = F)
    }
    
    
    
    ## how many rows are just zero? 
    justzero=sum(sonb$dxy==0)
    percentzero=justzero/nrow(sonb)
    print(paste("Proportion that is zero:",round(percentzero,2)))
    if(percentzero==0) {
      print("WARNING: NEED TO ADD ZEROES WHEN CALCULATING")
      ## need length of chromosomes
      addZeroes=T
    } else {
      addZeroes=F
    }
    
    print("READING IN SCAFFOLDS")
    num_on_scaff = rev(sort(table(sonb$chromo)))
    totalscaffs = length(num_on_scaff)
    print(paste("TOTAL SCAFFS READ:",totalscaffs))
    
    if (totalscaffs != num_to_keep) {
      num_to_keep=totalscaffs
    }
    
    ## for each chromosome?
    
    for(chr_i in 1:length(chromlist)){
      chrom = chromlist[chr_i]
      chromlength = chromLengths[chr_i]
      subset = sonb[sonb$chromo==chrom,]
      
      if(nrow(subset)>0) {
        
        snps_present = nrow(subset)
        snps_sum = sum(subset$dxy)
        if(addZeroes==T){
          
          chrom_average = snps_sum / chromlength
          
        } else {
          
          chrom_average = snps_sum / snps_present
        }
        
        ## start windows
        window_ends = seq(windowsize,max(chromlength,windowsize),movesize)
        window_starts = seq(1,chromlength,movesize)
        window_starts = window_starts[1:length(window_ends)]
        windows = as.data.frame(cbind("starts"=window_starts,"ends"=window_ends))
        number_windows = nrow(windows)
        
        ## loop over windows
        for(win_i in 1:nrow(windows)){
          window = windows[win_i,]
          if(win_i %% printstatus == 0) { print(paste(chrom,win_i,"/",number_windows)) }
          start = window$starts
          end = window$ends
          windows_subset = subset[subset$position>=start & subset$position<=end,]
          window_snps_sum = sum(windows_subset$dxy,na.rm=T)
          windows_snps_present = nrow(windows_subset)
          
          if(addZeroes==T){
            window_average = window_snps_sum / windowsize
            dxylist = c(windows_subset$dxy,rep(0,windowsize-windows_snps_present))
            dif = dxylist-window_average
            difsq = dif^2
            sumdifsq = sum(difsq)
            window_stdev = sqrt(sumdifsq/windowsize)
          } else {
            window_average = window_snps_sum / windows_snps_present
            dif = windows_subset$dxy-window_average
            difsq = dif^2
            sumdifsq = sum(difsq)
            window_stdev = sqrt(sumdifsq/windows_snps_present)
          }
          
          
          
          
          if(chr_i==1 && win_i==1){
            write.table(cbind("scafs","starts","means","stdvs","snps","sums"),outfile,
                        row.names = F,col.names = F,append=T)
          }
          
          write.table(cbind(chrom,start,window_average,window_stdev,windows_snps_present,window_snps_sum),
                      outfile,
                      append=T,row.names = F,col.names = F)
          
        }
      } 
    }
  }
}

