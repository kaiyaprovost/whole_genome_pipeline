## TODO: REORGANIZE

library(AICcmodavg)

steps=c(7,8) ## up to 13

setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/")

specieslist = sort(c("bil","bel",
  "fla",
  "fus","cur","bru","cri",
  "mel","nit","sin"))
size = 100

makedxyplot = function(speciesname,variable1num=4,variable2num=10,colvarnum=19,dat=bigplot,plotsize="SMALL",
                       lowquan1=0.05,highquan1=0.95,lowquan2=0.05,highquan2=0.95){
  pngname = paste(speciesname,
                  names(dat)[variable1num],
                  names(dat)[variable2num],
                  plotsize,".may2020.png",
                  sep="_")
  
  palette(    c(      "black",      "red",
                      "darkred",      "blue",      "purple",      "magenta",
                      "darkblue",      "orange",      "cyan"    )  )
  
  dat = dat[dat$species==speciesname,]
  var1 = dat[,variable1num]
  var2 = dat[,variable2num]
  cols = dat[,colvarnum]
  
  
  quan1 = quantile(var1, lowquan1)
  quan2 = quantile(var1, highquan1)
  quan3 = quantile(var2, lowquan2)
  quan4 = quantile(var2, highquan2)
  
  png(pngname)
  if (plotsize=="SMALL") {
    plot(var1, var2, col = as.numeric(cols),main=speciesname,
         xlab=paste(names(dat)[variable1num],lowquan1,highquan1),
         ylab=paste(names(dat)[variable2num],lowquan2,highquan2))
    text(
      x = mean(var1),
      y = mean(var2),
      labels = paste(
        "Mean",names(dat)[variable1num],
        round(mean(var1), 2),
        "\nMean",names(dat)[variable2num],
        round(mean(var2), 2),
        sep = " "
      ),
      col = "white"
    )
  } else {
    plot(var1, var2, col = as.numeric(cols),main=speciesname,
         xlab=paste(names(dat)[variable1num],lowquan1,highquan1),
         ylab=paste(names(dat)[variable2num],lowquan2,highquan2),
         xlim=c(0,1),ylim=c(0,1))
    
    text(
      x = 0.8,
      y = 0.8,
      labels = paste(
        "Mean",names(dat)[variable1num],
        round(mean(var1), 2),
        "\nMean",names(dat)[variable2num],
        round(mean(var2), 2),
        sep = " "
      ),
      col = "black"
    )
    
  }
  abline(v = quan1, col = "grey")
  abline(v = quan2, col = "grey")
  abline(h = quan3, col = "grey")
  abline(h = quan4, col = "grey")
  
  
  
  
  dev.off()
}
makedxyplotSTDEV = function(speciesname,variable1num=4,variable2num=10,colvarnum=19,dat=bigplot,plotsize="SMALL",numsd=5){
  pngname = paste(speciesname,
                  names(dat)[variable1num],
                  names(dat)[variable2num],
                  plotsize,"STDEV",numsd,".may2020.png",
                  sep="_")
  
  palette(    c(      "black",      "red",
                      "darkred",      "blue",      "purple",      "magenta",
                      "darkblue",      "orange",      "cyan"    )  )
  
  dat = dat[dat$species==speciesname,]
  var1 = dat[,variable1num]
  var2 = dat[,variable2num]
  cols = dat[,colvarnum]
  
  mean1 = mean(var1,na.rm = T)
  mean2 = mean(var2,na.rm = T)
  sd1 = sd(var1,na.rm = T)
  sd2 = sd(var2,na.rm = T)
  
  quan1 = mean1-(numsd*sd1)
  quan2 = mean1+(numsd*sd1)
  quan3 = mean2-(numsd*sd2)
  quan4 = mean2+(numsd*sd2)
  
  png(pngname)
  if (plotsize=="SMALL") {
    plot(var1, var2, col = as.numeric(cols),main=speciesname,
         xlab=paste(names(dat)[variable1num],"STDEVS:",numsd),
         ylab=paste(names(dat)[variable2num],"STDEVS:",numsd))
    text(
      x = mean(var1),
      y = mean(var2),
      labels = paste(
        "Mean",names(dat)[variable1num],
        round(mean(var1), 2),
        "\nMean",names(dat)[variable2num],
        round(mean(var2), 2),
        sep = " "
      ),
      col = "white"
    )
  } else {
    plot(var1, var2, col = as.numeric(cols),main=speciesname,
         xlab=paste(names(dat)[variable1num],"STDEVS:",numsd),
         ylab=paste(names(dat)[variable2num],"STDEVS:",numsd),
         xlim=c(0,1),ylim=c(0,1))
    
    text(
      x = 0.8,
      y = 0.8,
      labels = paste(
        "Mean",names(dat)[variable1num],
        round(mean(var1), 2),
        "\nMean",names(dat)[variable2num],
        round(mean(var2), 2),
        sep = " "
      ),
      col = "black"
    )
    
  }
  abline(v = quan1, col = "grey")
  abline(v = quan2, col = "grey")
  abline(h = quan3, col = "grey")
  abline(h = quan4, col = "grey")
  
  
  
  
  dev.off()
}


if(1 %in% steps) {
  print("STEP 1 dxy div fst"); {
    for (spp in sort(specieslist)) {
      print("step 1")
      print(spp)
      
      if (spp == "bel") {
        fullspecies = "Vireo_bellii-NOWEIRD"
      } else if (spp == "bil") {
        fullspecies = "Amphispiza_bilineata"
      } else if (spp == "bru") {
        fullspecies = "Campylorhynchus_brunneicapillus"
      } else if (spp == "cri") {
        fullspecies = "Toxostoma_crissale"
      } else if (spp == "cur") {
        fullspecies = "Toxostoma_curvirostre"
      } else if (spp == "fla") {
        fullspecies = "Auriparus_flaviceps-NOWEIRDMIN1C"
      } else if (spp == "fus") {
        fullspecies = "Melozone_fusca"
      } else if (spp == "mel") {
        fullspecies = "Polioptila_melanura"
      } else if (spp == "nit") {
        fullspecies = "Phainopepla_nitens"
      } else if (spp == "sin") {
        fullspecies = "Cardinalis_sinuatus"
      }
      print(fullspecies)
      
      dxyfile = paste(
        #"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/scaff_",
        #size,"/",spp,"_SON_Dxy_WINDOWS_",size,".txt",
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/scaff_chromosomes/",spp,"/",spp,"_4_SON_Dxy_WINDOWS_chrfix_1-ALL.txt",
        sep = "")
      fstfile = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/FST/3.slidingWindowJobs/SON_CHI_",
        fullspecies,
        "_FST_slidingwindow_chrfix.fst",
        sep = ""
      )
      
      if(!(file.exists(fstfile))) {
        fstfile = paste(
          "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/FST/3.slidingWindowJobs/SON_CHI_",
          fullspecies,
          "_FST_slidingwindow.fst",
          sep = ""
        )
      }
      
      outfile = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/scaff_",
        spp,"_bothFST_DXY_ALL.may2020.txt",
        sep = ""
      )
      pngfile0 = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/images/scaff_",
        spp,"_corDXYFST_ALL.may2020.png",
        sep = ""
      )
      pngfile1 = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/images/scaff_",
        spp,"_DXYdivFST_ALL.may2020.png",
        sep = ""
      )
      pngfile2 = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/images/scaff_",
        spp,"_FSTdivDXY_ALL.may2020.png",
        sep = ""
      )
      pngfile3 = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/",
        spp,"_panelFSTDXY_ALL.may2020.png",
        sep = ""
      )
      pngfile4 = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/images/scaff_",
        spp,"_panelFSTDXYsums_ALL.may2020.png",
        sep = ""
      )
      
      print("reading")
      dxy = read.csv(dxyfile, sep = " ",row.names = NULL)
      dxy=dxy[,c("scafs","starts","means","stdvs","snps","sums")]
      
      #newscafs=sapply(as.character(dxy$scafs),FUN=function(x){strsplit(x,"_")[[1]][4]},simplify = T)
      #dxy$scafs=newscafs    
      
      lines = readLines(fstfile)
      end = substr(lines[1],
                   start = nchar(lines[1]) - 2,
                   stop = nchar(lines[1]))
      
      if (end != "Fst") {
        newline1 = paste(lines[1], "\tFst", sep = "")
        lines[1] = newline1
        writeLines(lines, fstfile)
      }
      
      fst = read.csv(fstfile, sep = "\t")
      
      dxy$midPos = dxy$starts + 24999 ## won't work now? 
      names(dxy)[1] = "chr"
      names(dxy)[3] = "dxymeans"
      
      print("calculating")
      
      if (length(intersect(dxy$chr,fst$chr)) == 0) {
        print("CHROMS DO NOT MATCH")
        print(unique(dxy$chr))
        print(unique(fst$chr))
      } else {
        
        
        
        
        both_dxyfst = merge(dxy, fst, by = c("chr", "midPos"))
        if (nrow(both_dxyfst) <= 0) {
          dxy$midPos = dxy$starts + 59999 ## won't work now? 
          both_dxyfst = merge(dxy, fst, by = c("chr", "midPos"),all=T)
        } else {
          both_dxyfst = merge(dxy, fst, by = c("chr", "midPos"),all=T)
        }
        both_dxyfst$dxy_div_fst = both_dxyfst$dxymeans / both_dxyfst$Fst
        both_dxyfst$fst_div_dxy = both_dxyfst$Fst / both_dxyfst$dxymeans ## USE THIS?
        
        both_dxyfst$fst_div_dxysum = both_dxyfst$Fst / both_dxyfst$sums ## USE THIS?
        both_dxyfst$dxysum_div_fst = both_dxyfst$sums / both_dxyfst$Fst
        
        
        print("sorting")
        
        orderneeded = names(sort(table(both_dxyfst$chr),decreasing =T))
        both_dxyfst$chr = factor(both_dxyfst$chr , levels = orderneeded)
        both_dxyfst = (both_dxyfst[order(both_dxyfst$chr, both_dxyfst$midPos),])
        
        print("png0"); {
          png(pngfile0, width = 700, height = 700)
          par(
            bg = NA,
            col.axis = "white",
            fg = "white",
            col.lab = "white",
            col.main = "white"
          )
          plot(
            both_dxyfst$Fst,
            both_dxyfst$dxymeans,
            ylim = c(0, 1),
            xlim = c(0, 1),
            main = paste(spp, "DXY vs FST")
          )
          abline(a = 0, b = 1, col = "red")
          dev.off() }
        
        print("png1");
        if (sum(!(is.na(as.numeric(both_dxyfst$dxy_div_fst)))) > 0) {
          png(pngfile1, width = 700, height = 350)
          par(
            bg = NA,
            col.axis = "white",
            fg = "white",
            col.lab = "white",
            col.main = "white"
          )
          palette(
            c(
              "red",
              "cyan",
              "goldenrod",
              "green",
              "blue",
              "purple",
              "blue",
              "black",
              "brown",
              "magenta"
            )
          )
          plot(
            as.numeric(both_dxyfst$dxy_div_fst),
            col = as.numeric(as.factor(both_dxyfst$chr)),
            cex = 0.2,
            main = spp,
            xlab = "Window (Scaffold)",
            ylab = "DXY div FST"
          )
          dev.off()
        }
        
        print("png2"); 
        if (sum(!(is.na(as.numeric(both_dxyfst$fst_div_dxy)))) > 0) {
          
          png(pngfile2, width = 700, height = 350)
          par(
            bg = NA,
            col.axis = "white",
            fg = "white",
            col.lab = "white",
            col.main = "white"
          )
          palette(
            c(
              "red",
              "cyan",
              "goldenrod",
              "green",
              "blue",
              "purple",
              "blue",
              "black",
              "brown",
              "magenta"
            )
          )
          plot(
            as.numeric(both_dxyfst$fst_div_dxy),
            col = as.numeric(as.factor(both_dxyfst$chr)),
            cex = 0.2,
            main = spp,
            xlab = "Window (Scaffold)",
            ylab = "FST div DXY"
          )
          dev.off()
        }
        
        print("png3"); {
          png(pngfile3, width = 700, height = 700)
          par(
            bg = NA,
            col.axis = "white",
            fg = "white",
            col.lab = "white",
            col.main = "white"
          )
          
          if (sum(!(is.na(as.numeric(both_dxyfst$fst_div_dxy)))) > 0) {
            par(mfrow = c(4, 1)) 
          } else {
            par(mfrow = c(2, 1)) 
          }
          palette(
            c(
              "red",
              "cyan",
              "goldenrod",
              "green",
              "blue",
              "purple",
              "blue",
              "black",
              "brown",
              "magenta"
            )
          )
          plot(
            as.numeric(both_dxyfst$dxymeans),
            col = as.numeric(as.factor(both_dxyfst$chr)),
            cex = 0.2,
            main = spp,
            xlab = "Window (Scaffold)",
            #ylim = c(0, 1),
            ylab = "DXY only"
            
          )
          plot(
            as.numeric(both_dxyfst$Fst),
            col = as.numeric(as.factor(both_dxyfst$chr)),
            cex = 0.2,
            main = spp,
            xlab = "Window (Scaffold)",
            #ylim = c(0, 1),
            ylab = "FST only"
          )
          if (sum(!(is.na(as.numeric(both_dxyfst$fst_div_dxy)))) > 0) {
            plot(
              as.numeric(both_dxyfst$dxy_div_fst),
              col = as.numeric(as.factor(both_dxyfst$chr)),
              cex = 0.2,
              main = spp,
              xlab = "Window (Scaffold)",
              ylab = "DXY div FST"
            )
            plot(
              as.numeric(both_dxyfst$fst_div_dxy),
              col = as.numeric(as.factor(both_dxyfst$chr)),
              cex = 0.2,
              main = spp,
              xlab = "Window (Scaffold)",
              ylab = "FST div DXY"
            )
          }
          dev.off()
        }
        
        print("png4"); {
          png(pngfile4, width = 700, height = 700)
          if (sum(!(is.na(as.numeric(both_dxyfst$fst_div_dxysum)))) > 0) {
            par(mfrow = c(4, 1)) 
          } else {
            par(mfrow = c(2, 1)) 
          }
          par(
            bg = NA,
            col.axis = "white",
            fg = "white",
            col.lab = "white",
            col.main = "white"
          )
          palette(
            c(
              "red",
              "cyan",
              "goldenrod",
              "green",
              "blue",
              "purple",
              "blue",
              "black",
              "brown",
              "magenta"
            )
          )
          plot(
            as.numeric(both_dxyfst$sums),
            col = as.numeric(as.factor(both_dxyfst$chr)),
            cex = 0.2,
            main = spp,
            xlab = "Window (Scaffold)",
            ylab = "DXYsums only"
          )
          plot(
            as.numeric(both_dxyfst$Fst),
            col = as.numeric(as.factor(both_dxyfst$chr)),
            cex = 0.2,
            main = spp,
            xlab = "Window (Scaffold)",
            #ylim = c(0, 1),
            ylab = "FST only"
          )
          if (sum(!(is.na(as.numeric(both_dxyfst$fst_div_dxysum)))) > 0) {
            plot(
              as.numeric(both_dxyfst$dxysum_div_fst),
              col = as.numeric(as.factor(both_dxyfst$chr)),
              cex = 0.2,
              main = spp,
              xlab = "Window (Scaffold)",
              ylab = "DXYsums div FST"
            )
            plot(
              as.numeric(both_dxyfst$fst_div_dxysum),
              col = as.numeric(as.factor(both_dxyfst$chr)),
              cex = 0.2,
              main = spp,
              xlab = "Window (Scaffold)",
              ylab = "FST div DXYsums"
            )
          }
        }
        dev.off()
        
        #plot(both_dxyfst$dxymeans,ylim=c(0,1))
        #points(both_dxyfst$Fst,col="red")
        
        write.table(both_dxyfst, outfile)
      }
      
    }
  }
}

if(2 %in% steps) {
  print("STEP 2 quantiles"); {
    for (spp in specieslist) {
      print("step 2")
      print(spp)
      testfile = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/scaff_",
        spp,"_bothFST_DXY_ALL.may2020.txt",
        sep = ""
      )
      #testfile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/sin_bothFST_DXY_1.txt"
      test = read.csv(testfile, sep = " ")
      #test = test[complete.cases(test), ]
      summary(test)
      quan1 = quantile(test$Fst, 0.05,na.rm=T)
      quan2 = quantile(test$Fst, 0.95,na.rm=T)
      test$quantileFST = 1 ## middle
      test$quantileFST[test$Fst <= quan1] = 2 ## bottom
      test$quantileFST[test$Fst >= quan2] = 4 ## top
      #plot(test$Fst,col=as.numeric(test$quantileFST),pch=as.numeric(test$quantileFST))
      
      quan3 = quantile(test$dxymeans, 0.05,na.rm=T)
      quan4 = quantile(test$dxymeans, 0.95,na.rm=T)
      test$quantiledxymeans = 8
      test$quantiledxymeans[test$dxymeans <= quan3] = 16
      test$quantiledxymeans[test$dxymeans >= quan4] = 32
      #plot(test$dxymeans,col=as.numeric(test$quantiledxymeans/8),pch=as.numeric(test$quantiledxymeans/8))
      
      
      test$sumquantile = test$quantiledxymeans + test$quantileFST
      ## 9  = 1+8  = fst not outlier, dxy not outlier
      ## 10 = 2+8  = fst low,         dxy not outlier
      ## 12 = 4+8  = fst high,        dxy not outlier
      ## 17 = 1+16 = fst not outlier, dxy low
      ## 18 = 2+16 = fst low,         dxy low
      ## 20 = 4+16 = fst high,        dxy low -- SWEEP
      ## 33 = 1+32 = fst not outlier, dxy high
      ## 34 = 2+32 = fst low,         dxy high
      ## 36 = 4+32 = fst high,        dxy high -- ISLAND
      
      test$quancolor = test$sumquantile
      test$quancolor[test$quancolor == 9] = 1
      test$quancolor[test$quancolor == 10] = 2
      test$quancolor[test$quancolor == 12] = 3
      test$quancolor[test$quancolor == 17] = 4
      test$quancolor[test$quancolor == 18] = 5
      test$quancolor[test$quancolor == 20] = 6
      test$quancolor[test$quancolor == 33] = 7
      test$quancolor[test$quancolor == 34] = 8
      test$quancolor[test$quancolor == 36] = 9
      
      png(
        paste(
          "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/",
          spp,
          "_fstVSdxy_quantilemap_ALL.may2020.png",
          sep = ""
        ), height=700,width=700
      )
      palette(  c("black",  "red",  "darkred",  "blue",   "purple",    "magenta",  "darkblue", "orange",    "cyan"))
      par(
        bg = NA,
        col.axis = "white",
        fg = "white",
        col.lab = "white",
        col.main = "white"
      )
      plot(
        test$dxymeans,
        test$Fst,
        col = test$quancolor,
        type = "n",
        main = spp
      )
      abline(h = quan1, col = "grey")
      abline(h = quan2, col = "grey")
      abline(v = quan3, col = "grey")
      abline(v = quan4, col = "grey")
      points(
        test$dxymeans,
        test$Fst,
        col = test$quancolor,
        pch = as.numeric(test$quancolor)
      )
      
      text(
        x = mean(test$dxymeans),
        y = mean(test$Fst),
        labels = paste(
          "Mean FST:",
          round(mean(test$Fst), 2),
          "\nMean DXY:",
          round(mean(test$dxymeans), 2),
          sep = " "
        ),
        col = "white"
      )
      
      dev.off()
      
      png(
        paste(
          "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/",
          spp,
          "_fstVSdxy_panelquantile_ALL.may2020.png",
          sep = ""
        ),
        width = 700,
        height = 700
      )
      par(mfrow = c(2, 1))
      par(
        bg = NA,
        col.axis = "white",
        fg = "white",
        col.lab = "white",
        col.main = "white"
      )
      plot(
        test$dxymeans,
        col = test$quancolor,
        pch = as.numeric(test$quancolor),
        main = spp
      )
      plot(
        test$Fst,
        col = test$quancolor,
        pch = as.numeric(test$quancolor),
        main = spp
      )
      #plot(test$dxy_div_fst,col=test$quancolor,pch=as.numeric(test$quancolor),main=spp)
      dev.off()
      
      
      
      png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/",
                spp,"_fstVSdxy_panelquantile_legible.may2020.png",sep = ""),
          width = 700,height = 700)
      palette(c("black",
                "red",          "darkred",
                "blue",          "purple",
                "magenta",          "darkblue",
                "orange",          "cyan"        )      )
      par(mfrow = c(2, 1))
      par(mar=c(3,4,0,0))
      plot(test$dxymeans,col = as.numeric(as.factor(test$chr)),
           pch = as.numeric(test$quancolor),main = "",
           ylab="DXY",xlab="")
      plot(test$Fst, col = as.numeric(as.factor(test$chr)),
           pch = as.numeric(test$quancolor),main = "",ylab="FST",xlab="Window\n")
      dev.off()
      
      
      
      freq = as.data.frame(table(test$sumquantile))
      # freq$VarNames = c(
      #   "NONE",
      #   "LOW_FST",
      #   "HIGH_FST",
      #   "LOW_DXY",
      #   "BOTH_LOW",
      #   "HIGH_FST_LOW_DXY",
      #   "HIGH_DXY",
      #   "LOW_FST_HIGH_DXY",
      #   "BOTH_HIGH"
      # )
      freq$VarNames=as.character(freq$Var1)
      freq$VarNames[freq$Var1==9]="NONE"
      freq$VarNames[freq$Var1==10]="LOW_FST"
      freq$VarNames[freq$Var1==12]="HIGH_FST"
      freq$VarNames[freq$Var1==17]="LOW_DXY"
      freq$VarNames[freq$Var1==18]="BOTH_LOW"
      freq$VarNames[freq$Var1==20]="HIGH_FST_LOW_DXY"
      freq$VarNames[freq$Var1==33]="HIGH_DXY"
      freq$VarNames[freq$Var1==34]="LOW_FST_HIGH_DXY"
      freq$VarNames[freq$Var1==36]="BOTH_HIGH"
      
      ## 9  = 1+8  = fst not outlier, dxy not outlier
      ## 10 = 2+8  = fst low,         dxy not outlier
      ## 12 = 4+8  = fst high,        dxy not outlier
      ## 17 = 1+16 = fst not outlier, dxy low
      ## 18 = 2+16 = fst low,         dxy low
      ## 20 = 4+16 = fst high,        dxy low
      ## 33 = 1+32 = fst not outlier, dxy high
      ## 34 = 2+32 = fst low,         dxy high
      ## 36 = 4+32 = fst high,        dxy high
      
      # Basic Barplot
      png(
        paste(
          "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/",
          spp,
          "_fstVSdxy_lognumoutliers_ALL.may2020.png",
          sep = ""
        )
      )
      par(mar = c(10, 4, 1, 1))
      par(
        bg = NA,
        col.axis = "white",
        fg = "white",
        col.lab = "white",
        col.main = "white"
      )
      my_bar = barplot(
        freq$Freq,
        names.arg = freq$VarNames,
        log = "y",
        las = 2,
        main = spp,
        ylim = c(1, 50000),
        col = c(
          "grey",
          "red",
          "darkred",
          "blue",
          "purple",
          "magenta",
          "darkblue",
          "orange",
          "cyan"
        )
      )
      # Add the text
      text(my_bar, 100 , paste(freq$Freq, sep = "") , cex = 1)
      dev.off()
      
      
      freq2 = as.data.frame(table(test$quantiledxymeans, test$quantileFST))
      
      
      
      
      
      #hist(test$sumquantile)
      #hist(test$quantiledxymeans)
      #hist(test$quantileFST)
      
      ## 9  = 1+8  = fst not outlier, dxy not outlier
      ## 10 = 2+8  = fst low,         dxy not outlier
      ## 12 = 4+8  = fst high,        dxy not outlier
      ## 17 = 1+16 = fst not outlier, dxy low
      ## 18 = 2+16 = fst low,         dxy low
      ## 20 = 4+16 = fst high,        dxy low
      ## 33 = 1+32 = fst not outlier, dxy high
      ## 34 = 2+32 = fst low,         dxy high
      ## 36 = 4+32 = fst high,        dxy high
      
      dxyhigh = test[test$quantiledxymeans == 32, ]
      dxylow = test[test$quantiledxymeans == 16, ]
      dxymid = test[test$quantiledxymeans == 8, ]
      
      highD_lowF = dxyhigh[dxyhigh$quantileFST == 2, ]
      ISLANDS = dxyhigh[dxyhigh$quantileFST == 4, ]
      highD_midF = dxyhigh[dxyhigh$quantileFST == 1, ]
      
      lowD_lowF = dxylow[dxylow$quantileFST == 2, ]
      SWEEPS = dxylow[dxylow$quantileFST == 4, ]
      lowD_midF = dxylow[dxylow$quantileFST == 1, ]
      
      midD_lowF = dxymid[dxymid$quantileFST == 2, ]
      midD_highF = dxymid[dxymid$quantileFST == 4, ]
      midD_midF = dxymid[dxymid$quantileFST == 1, ]
      
      notISLANDS = rbind(dxymid, dxylow,
                         highD_lowF, highD_midF)
      
      notSWEEPS = rbind(dxymid, dxyhigh,
                        lowD_lowF, lowD_midF)
      
      
      head(ISLANDS)
      #par(ask=T)
      png(
        paste(
          "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/",
          spp,
          "_boxplotsIslandsSweeps_ALL.may2020.png",
          sep = ""
        )
      )
      par(mfrow = c(1, 2))
      par(
        bg = NA,
        col.axis = "white",
        fg = "white",
        col.lab = "white",
        col.main = "white"
      )
      boxplot(
        ISLANDS$dxymeans,
        notISLANDS$dxymeans,
        SWEEPS$dxymeans,
        notSWEEPS$dxymeans,
        main = spp,
        ylab = "mean DXY",
        ylim = c(0, 0.5),
        names = c("ISLANDS", "NOT ISLANDS",
                  "SWEEPS", "NOT SWEEPS"),
        las = 2
      )
      
      boxplot(
        ISLANDS$Fst,
        notISLANDS$Fst,
        SWEEPS$Fst,
        notSWEEPS$Fst,
        main = spp,
        ylab = "mean FST",
        ylim = c(0, 0.4),
        names = c("ISLANDS", "NOT ISLANDS",
                  "SWEEPS", "NOT SWEEPS"),
        las = 2
      )
      #par(ask=F)
      dev.off()
      
      
      
      
      sweepfile = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/",
        spp,
        "_sweeps_ALL.may2020.png",
        sep = ""
      )
      
      png(sweepfile)
      par(
        bg = NA,
        col.axis = "white",
        fg = "white",
        col.lab = "white",
        col.main = "white"
      )
      par(mar = c(4.5, 4, 0, 0),
          mfrow = c(2, 1))
      
      palette(c("grey","darkgrey"))
      plot(
        as.numeric(test[, "Fst"]),
        col = as.numeric(as.factor(test[,"chr"])),
        cex = 0.2,
        ylab = "Fst",
        type = "p"
      )
      palette(c(rgb(0, 0, 0, 0), "cyan"))
      points(as.numeric(test[, "Fst"]),
             col = 1 + as.numeric(test[, "sumquantile"] == 36),
             pch = 16)
      
      palette(c("grey","darkgrey"))
      plot(
        as.numeric(test[, "Fst"]),
        col = as.numeric(as.factor(test[,"chr"])),
        cex = 0.2,
        ylab = "Fst",
        type = "p"
      )
      palette(c(rgb(0, 0, 0, 0), "magenta"))
      points(as.numeric(test[, "Fst"]),
             col = 1 + as.numeric(test[, "sumquantile"] == 20),
             pch = 16)
      
      dev.off()
      
      
      ## sweeps including the burgundy ones 
      sweepexpfile = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/",
        spp,
        "_sweeps_expanded_ALL.may2020.png",
        sep = ""
      )
      
      png(sweepexpfile)
      par(mar = c(4.5, 4, 0, 0),
          mfrow = c(3, 1))
      
      palette(c("grey","darkgrey"))
      plot(
        as.numeric(test[, "Fst"]),
        col = as.numeric(as.factor(test[,"chr"])),
        cex = 0.2,
        ylab = "Fst",
        type = "p"
      )
      palette(c(rgb(0, 0, 0, 0), "cyan"))
      points(as.numeric(test[, "Fst"]),
             col = 1 + as.numeric(test[, "sumquantile"] == 36),
             pch = 16)
      
      palette(c("grey","darkgrey"))
      plot(
        as.numeric(test[, "Fst"]),
        col = as.numeric(as.factor(test[,"chr"])),
        cex = 0.2,
        ylab = "Fst",
        type = "p"
      )
      palette(c(rgb(0, 0, 0, 0), "magenta"))
      points(as.numeric(test[, "Fst"]),
             col = 1 + as.numeric(test[, "sumquantile"] == 20),
             pch = 16)
      
      palette(c("grey","darkgrey"))
      plot(
        as.numeric(test[, "Fst"]),
        col = as.numeric(as.factor(test[,"chr"])),
        cex = 0.2,
        ylab = "Fst",
        type = "p"
      )
      palette(c(rgb(0, 0, 0, 0), "darkred"))
      points(as.numeric(test[, "Fst"]),
             col = 1 + as.numeric(test[, "sumquantile"] == 12),
             pch = 16)
      
      dev.off()
      
      
      ## sweeps with dxy
      sweepdxyfile = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/",
        spp,
        "_sweeps_dxy_ALL.may2020.png",
        sep = ""
      )
      
      png(sweepdxyfile)
      par(mar = c(4.5, 4, 0, 0),
          mfrow = c(3, 1))
      
      palette(c("grey","darkgrey"))
      plot(
        as.numeric(test[, "dxymeans"]),
        col = as.numeric(as.factor(test[,"chr"])),
        cex = 0.2,
        ylab = "DXY",
        type = "p"
      )
      palette(c(rgb(0, 0, 0, 0), "cyan"))
      points(as.numeric(test[, "dxymeans"]),
             col = 1 + as.numeric(test[, "sumquantile"] == 36),
             pch = 16)
      
      palette(c("grey","darkgrey"))
      plot(
        as.numeric(test[, "dxymeans"]),
        col = as.numeric(as.factor(test[,"chr"])),
        cex = 0.2,
        ylab = "Dxy",
        type = "p"
      )
      palette(c(rgb(0, 0, 0, 0), "magenta"))
      points(as.numeric(test[, "dxymeans"]),
             col = 1 + as.numeric(test[, "sumquantile"] == 20),
             pch = 16)
      
      palette(c("grey","darkgrey"))
      plot(
        as.numeric(test[, "dxymeans"]),
        col = as.numeric(as.factor(test[,"chr"])),
        cex = 0.2,
        ylab = "Dxy",
        type = "p"
      )
      palette(c(rgb(0, 0, 0, 0), "darkred"))
      points(as.numeric(test[, "dxymeans"]),
             col = 1 + as.numeric(test[, "sumquantile"] == 12),
             pch = 16)
      
      dev.off()
      
      
      
    }
  }
}

if(3 %in% steps) {
  print("STEP 3 larger panels"); {
    for (spp in specieslist) {
      print(spp)
      print("step 3")
      infile = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/scaff_",
        spp,
        "_bothFST_DXY_ALL.may2020.txt",
        sep = ""
      )
      
      outfile = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/",
        spp,
        "_bothFST_DXY_ALL.may2020.png",
        sep = ""
      )
      
      
      textfiles = read.csv(infile, sep = " ")
      
      palette(c(
        "red",
        "cyan",
        "goldenrod",
        "green",
        "blue",
        "purple",
        "magenta"
      ))
      
      png(outfile)
      par(
        bg = NA,
        col.axis = "white",
        fg = "white",
        col.lab = "white",
        col.main = "white"
      )
      par(mfrow = c(2, 1),
          mar = c(4.5, 4, 0.3, 0.3))
      plot(
        as.numeric(textfiles[, "Fst"]),
        col = as.numeric(as.factor(textfiles[, "chr"])),
        cex = 0.2,
        ylab = "Fst",
        ylim = c(0, 1),
        type = "p"
      )
      plot(
        as.numeric(textfiles[, "dxymeans"]),
        col = as.numeric(as.factor(textfiles[, "chr"])),
        cex = 0.2,
        ylab = "Dxy",
        ylim = c(0, 1),
        type = "p"
      )
      dev.off()
      
      
      
      
    }
  }
}

## GET PLOTS OF SPP ALL ON SAME AXES
print("BIGPLOT")

bigplot = data.frame()
if(4 %in% steps) {
  print("STEP 4 generate bigplot"); {
    for (i in 1:length(specieslist)) {
      spp = specieslist[i]
      print(spp)
      print("step 4")
      testfile = paste(
        "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/scaff_",
        spp,
        "_bothFST_DXY_ALL.may2020.txt",
        sep = ""
      )
      #testfile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/sin_bothFST_DXY_1.txt"
      test = read.csv(testfile, sep = " ")
      test = test[complete.cases(test), ]
      test$species = spp
      
      quan1 = quantile(test$Fst, 0.05)
      quan2 = quantile(test$Fst, 0.95)
      test$quantileFST = 1 ## middle
      test$quantileFST[test$Fst <= quan1] = 2 ## bottom
      test$quantileFST[test$Fst >= quan2] = 4 ## top
      #plot(test$Fst,col=as.numeric(test$quantileFST),pch=as.numeric(test$quantileFST))
      
      quan3 = quantile(test$dxymeans, 0.05)
      quan4 = quantile(test$dxymeans, 0.95)
      test$quantiledxymeans = 8
      test$quantiledxymeans[test$dxymeans <= quan3] = 16
      test$quantiledxymeans[test$dxymeans >= quan4] = 32
      #plot(test$dxymeans,col=as.numeric(test$quantiledxymeans/8),pch=as.numeric(test$quantiledxymeans/8))
      
      
      test$sumquantile = test$quantiledxymeans + test$quantileFST
      ## 9  = 1+8  = fst not outlier, dxy not outlier
      ## 10 = 2+8  = fst low,         dxy not outlier
      ## 12 = 4+8  = fst high,        dxy not outlier
      ## 17 = 1+16 = fst not outlier, dxy low
      ## 18 = 2+16 = fst low,         dxy low
      ## 20 = 4+16 = fst high,        dxy low
      ## 33 = 1+32 = fst not outlier, dxy high
      ## 34 = 2+32 = fst low,         dxy high
      ## 36 = 4+32 = fst high,        dxy high
      
      test$quancolor = test$sumquantile
      test$quancolor[test$quancolor == 9] = 1
      test$quancolor[test$quancolor == 10] = 2
      test$quancolor[test$quancolor == 12] = 3
      test$quancolor[test$quancolor == 17] = 4
      test$quancolor[test$quancolor == 18] = 5
      test$quancolor[test$quancolor == 20] = 6
      test$quancolor[test$quancolor == 33] = 7
      test$quancolor[test$quancolor == 34] = 8
      test$quancolor[test$quancolor == 36] = 9
      
      if (i == 1) {
        bigplot = test
      } else {
        bigplot = rbind(test, bigplot)
        
      }
      
      
    }
  }
}
if(5 %in% steps) {
  print("STEP 5 bigplot output"); {
    
    head(bigplot)
    summary(bigplot$species)
    unique(bigplot$species)
    
    summary(bigplot$dxymeans)
    
    quan1 = quantile(bigplot$Fst, 0.05,na.rm=T)
    quan2 = quantile(bigplot$Fst, 0.95,na.rm=T)
    quan3 = quantile(bigplot$dxymeans, 0.05,na.rm=T)
    quan4 = quantile(bigplot$dxymeans, 0.95,na.rm=T)
    
    
    morphcolors = rbind(
      c("bel", 1, 1, 1, 1, 1),
      c("bil", 0, 1, 1, 1, 0),
      c("bru", 0, 1, 1, 1, 0),
      c("cri", 1, 1, 1, 0, 1),
      c("cur", 0, 1, 1, 0, 1),
      c("fla", 0, 1, 1, 0, 1),
      c("fus", 0, 0, 0, 1, 1),
      c("mel", 0, 0, 1, 0, 0),
      #c("nit",0,1,1,0,0),
      c("sin", 0, 0, 0, 1, 0)
    )
    
    
    
    
    #nrow(bigplot)
    #nrow(unique(bigplot))
    png("bigplot_biglim.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    plot(bigplot$dxymeans, bigplot$Fst#,xlim=c(0,1),ylim=c(0,1)
    )
    dev.off()
    
    bigplot$globalquantileFST = 1 ## middle
    bigplot$globalquantileFST[bigplot$Fst <= quan1] = 2 ## bottom
    bigplot$globalquantileFST[bigplot$Fst >= quan2] = 4 ## top
    bigplot$globalquantiledxymeans = 8
    bigplot$globalquantiledxymeans[bigplot$dxymeans <= quan3] = 16
    bigplot$globalquantiledxymeans[bigplot$dxymeans >= quan4] = 32
    
    png("bigplot_color_biglim.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    plot(bigplot$dxymeans, bigplot$Fst, col = as.numeric(bigplot$quancolor)#,xlim=c(0,1),ylim=c(0,1)
    )
    dev.off()
    
    bigplot$globalsumquantile = bigplot$globalquantiledxymeans + bigplot$globalquantileFST
    bigplot$globalquancolor = bigplot$globalsumquantile
    bigplot$globalquancolor[bigplot$globalquancolor == 9] = 1
    bigplot$globalquancolor[bigplot$globalquancolor == 10] = 2
    bigplot$globalquancolor[bigplot$globalquancolor == 12] = 3
    bigplot$globalquancolor[bigplot$globalquancolor == 17] = 4
    bigplot$globalquancolor[bigplot$globalquancolor == 18] = 5
    bigplot$globalquancolor[bigplot$globalquancolor == 20] = 6
    bigplot$globalquancolor[bigplot$globalquancolor == 33] = 7
    bigplot$globalquancolor[bigplot$globalquancolor == 34] = 8
    bigplot$globalquancolor[bigplot$globalquancolor == 36] = 9
    
    png("bigplot_color_biglim.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    plot(bigplot$dxymeans, bigplot$Fst, col = as.numeric(bigplot$quancolor)#,xlim=c(0,1),ylim=c(0,1)
    )
    dev.off()
    
    png("bigplot_globalcolor_biglim.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    plot(bigplot$dxymeans,
         bigplot$Fst,
         col = as.numeric(bigplot$globalquancolor)#,xlim=c(0,1),ylim=c(0,1)
    )
    dev.off()
    
    ## 9  = 1+8  = fst not outlier, dxy not outlier
    ## 10 = 2+8  = fst low,         dxy not outlier
    ## 12 = 4+8  = fst high,        dxy not outlier
    ## 17 = 1+16 = fst not outlier, dxy low
    ## 18 = 2+16 = fst low,         dxy low
    ## 20 = 4+16 = fst high,        dxy low
    ## 33 = 1+32 = fst not outlier, dxy high
    ## 34 = 2+32 = fst low,         dxy high
    ## 36 = 4+32 = fst high,        dxy high
    
    png("bigplot_quanVSquan_biglim.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    plot(bigplot$sumquantile, bigplot$globalsumquantile#,xlim=c(0,1),ylim=c(0,1)
    )
    dev.off()
    
    palette(c("white", "pink"))
    png("bigplot_boxplotdxy_islands.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$dxymeans[bigplot$sumquantile == 36] ~ bigplot$species[bigplot$sumquantile ==
                                                                      36],
      ylab = "DXY (per species quantiles)",
      col = as.numeric(morphcolors[, 6]) + 1
    )
    dev.off()
    
    png("bigplot_boxplotdxy_notislands.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$dxymeans[bigplot$sumquantile != 36] ~ bigplot$species[bigplot$sumquantile !=
                                                                      36],
      ylab = "DXY (per species quantiles)",
      col = as.numeric(morphcolors[, 6]) + 1
    )
    dev.off()
    
    palette(
      c(
        "black",
        "red",
        "green",
        "blue",
        "purple",
        "goldenrod",
        "blue",
        "purple",
        "magenta"
      )
    )
    
    png("bigplot_boxplotfst_islands_majorpca.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$Fst[bigplot$sumquantile == 36] ~ bigplot$species[bigplot$sumquantile ==
                                                                 36],
      ylab = "FST (per species quantiles)",
      main = "SIG IN PC1-PC3",
      border = as.numeric(morphcolors[, 2]) + 1
    )
    dev.off()
    
    png("bigplot_boxplotfst_islands_allpca.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$Fst[bigplot$sumquantile == 36] ~ bigplot$species[bigplot$sumquantile ==
                                                                 36],
      ylab = "FST (per species quantiles)",
      main = "SIG IN ANY PC",
      border = as.numeric(morphcolors[, 3]) + 1
    )
    dev.off()
    png("bigplot_boxplotfst_islands_rawmetric.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$Fst[bigplot$sumquantile == 36] ~ bigplot$species[bigplot$sumquantile ==
                                                                 36],
      ylab = "FST (per species quantiles)",
      main = "SIG IN RAW METRICS",
      border = as.numeric(morphcolors[, 4]) + 1
    )
    dev.off()
    
    png("bigplot_boxplotfst_islands_residualmetric.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$Fst[bigplot$sumquantile == 36] ~ bigplot$species[bigplot$sumquantile ==
                                                                 36],
      ylab = "FST (per species quantiles)",
      main = "SIG IN METRICS AFTER CORRECT FOR SIZE",
      border = as.numeric(morphcolors[, 5]) + 1
    )
    dev.off()
    
    palette(c("white", "pink"))
    png("bigplot_boxplotfst_notislands.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$Fst[bigplot$sumquantile != 36] ~ bigplot$species[bigplot$sumquantile !=
                                                                 36],
      ylab = "FST (per species quantiles)",
      col = as.numeric(morphcolors[, 6]) + 1
    )
    dev.off()
    
    
    png("bigplot_boxplotDXY_islands_majorpca.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$dxymeans[bigplot$sumquantile == 36] ~ bigplot$species[bigplot$sumquantile ==
                                                                      36],
      ylab = "DXY (per species quantiles)",
      main = "SIG IN PC1-PC3",
      border = as.numeric(morphcolors[, 2]) + 1
    )
    dev.off()
    png("bigplot_boxplotDXY_islands_allpca.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$dxymeans[bigplot$sumquantile == 36] ~ bigplot$species[bigplot$sumquantile ==
                                                                      36],
      ylab = "DXY (per species quantiles)",
      main = "SIG IN ANY PC",
      border = as.numeric(morphcolors[, 3]) + 1
    )
    dev.off()
    png("bigplot_boxplotDXY_islands_rawmetric.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$dxymeans[bigplot$sumquantile == 36] ~ bigplot$species[bigplot$sumquantile ==
                                                                      36],
      ylab = "DXY (per species quantiles)",
      main = "SIG IN RAW METRICS",
      border = as.numeric(morphcolors[, 4]) + 1
    )
    dev.off()
    png("bigplot_boxplotDXY_islands_residualmetric.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$dxymeans[bigplot$sumquantile == 36] ~ bigplot$species[bigplot$sumquantile ==
                                                                      36],
      ylab = "DXY (per species quantiles)",
      main = "SIG IN METRICS AFTER CORRECT FOR SIZE",
      border = as.numeric(morphcolors[, 5]) + 1
    )
    dev.off()
    
    palette(c("white", "pink"))
    png("bigplot_boxplotDXY_notislands.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$dxymeans[bigplot$sumquantile != 36] ~ bigplot$species[bigplot$sumquantile !=
                                                                      36],
      ylab = "DXY (per species quantiles)",
      col = as.numeric(morphcolors[, 6]) + 1
    )
    dev.off()
    
    
    if(length(bigplot$dxymeans[bigplot$globalsumquantile == 36])>0){
    png("bigplot_boxplotdxy_globalislands.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$dxymeans[bigplot$globalsumquantile == 36] ~ bigplot$species[bigplot$globalsumquantile ==
                                                                            36],
      ylab = "DXY (global quantiles)",
      col = as.numeric(morphcolors[, 6]) + 1
    )
    dev.off()
    }
    
    png("bigplot_boxplotdxy_globalnotislands.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$dxymeans[bigplot$globalsumquantile != 36] ~ bigplot$species[bigplot$globalsumquantile !=
                                                                            36],
      ylab = "DXY (global quantiles)",
      col = as.numeric(morphcolors[, 6]) + 1
    )
    dev.off()
    
    if(length(bigplot$Fst[bigplot$globalsumquantile == 36])>0){
    png("bigplot_boxplotfst_globalislands.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$Fst[bigplot$globalsumquantile == 36] ~ bigplot$species[bigplot$globalsumquantile ==
                                                                       36],
      ylab = "FST (global quantiles)",
      col = as.numeric(morphcolors[, 6]) + 1
    )
    dev.off()
    }
    
    png("bigplot_boxplotfst_globalnotislands.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    boxplot(
      bigplot$Fst[bigplot$globalsumquantile != 36] ~ bigplot$species[bigplot$globalsumquantile !=
                                                                       36],
      ylab = "FST (global quantiles)",
      col = as.numeric(morphcolors[, 6]) + 1
    )
    dev.off()
    
    write.csv(bigplot, "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.may2020.txt", row.names = F)
    
    
    
  }
}
if(6 %in% steps) {
  print("STEP 6 bigplot pngs"); {
    
    bigplot = read.csv(
      "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.may2020.txt"
    )
    
    nummorph = c(10,
                 5,
                 5,
                 14,
                 3,
                 5,
                 9,
                 2,
                 #3,
                 4)
    
    tab = (table(bigplot$sumquantile, bigplot$species))
    sppsums = colSums(tab)
    div = t(t(tab) / sppsums * 100)
    rownames(div) = c(
      "Normal",
      "LowFST",
      "HighFST",
      "LowDXY",
      "LowBoth",
      "Sweep",
      "HighDxy",
      "LowFSTHighDXY",
      "Island"
    )
    colSums(div)
    div_no9 = div[2:9,]
    palette(
      c(
        "black",
        "red",
        "darkred",
        "blue",
        "purple",
        "magenta",
        "darkblue",
        "orange",
        "cyan"
      )
    )
    png("bigplot_sitesperspp_local_all.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    barplot(
      div,
      ylab = "Percent of Sites",
      xlab = "Species",
      col = c(
        "black",
        "red",
        "darkred",
        "blue",
        "purple",
        "magenta",
        "darkblue",
        "orange",
        "cyan"
      ),
      ylim = c(0, 125)
    )
    legend(
      "topleft",
      legend = rownames(div),
      fill = c(
        "black",
        "red",
        "darkred",
        "blue",
        "purple",
        "magenta",
        "darkblue",
        "orange",
        "cyan"
      ),
      ncol = 3,
      border = "white",
      bty = "n"
    )
    dev.off()
    
    png("bigplot_sitesperspp_local_nocommon.may2020.png",width=1200)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    barplot(
      div_no9,
      ylab = "Percent of Sites",
      xlab = "Species",
      col = c(
        "red",
        "darkred",
        "blue",
        "purple",
        "magenta",
        "darkblue",
        "orange",
        "cyan"
      ),
      ylim = c(0, 25)
    )
    legend(
      "topleft",
      legend = rownames(div_no9),
      fill = c(
        "red",
        "darkred",
        "blue",
        "purple",
        "magenta",
        "darkblue",
        "orange",
        "cyan"
      ),
      ncol = 3,
      border = "white",
      bty = "n"
    )
    dev.off()
    
    div_swpIsl = div[c(6, 9),]
    #div_swpIsl = div_swpIsl[, c("bel", "cri", "cur", "fla", "fus", "bil", "bru", "mel", "nit", "sin")]
    div_swpIsl = div_swpIsl[, c("cri", "cur", "fus", "bil", "bru", "mel", "nit", "sin")]
    
    div_swpIsl = as.matrix(div_swpIsl)
    png("bigplot_sitestack_per_spp.may2020.png",
        width = 700,
        height = 300)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    par(mar = c(4.5, 4, 0, 0))
    barplot((div_swpIsl),
            beside = F,
            col = c("magenta", "cyan"),
            ylab = "Percent of Sites",
            xlab = "Species"
    )
    legend(
      "topleft",
      legend = c(rownames(div_swpIsl)),
      fill = c("magenta", "cyan"),
      border = "white",
      bty = "n"
    )
    dev.off()
    
    tab2 = as.data.frame(cbind(tab, nummorph))
    
    plot(tab2)
    corrplot::corrplot(cor(tab2), diag = F, method = "ellipse")
    plot(tab2$nummorph, tab2$`12`)
    
    ## 9  = 1+8  = fst not outlier, dxy not outlier
    ## 10 = 2+8  = fst low,         dxy not outlier
    ## 12 = 4+8  = fst high,        dxy not outlier
    ## 17 = 1+16 = fst not outlier, dxy low
    ## 18 = 2+16 = fst low,         dxy low
    ## 20 = 4+16 = fst high,        dxy low
    ## 33 = 1+32 = fst not outlier, dxy high
    ## 34 = 2+32 = fst low,         dxy high
    ## 36 = 4+32 = fst high,        dxy high
    
    
    ## get number and size of regions with islans or islands
    ## sweep = 20
    ## islands = 36
  }
}
if(7 %in% steps) {
  print("STEP 7 bigplot island size calc"); {
    for (spp in unique(bigplot$species)) {
      #print(spp)
      #print("step 7")
      subbig = bigplot[bigplot$species == spp,]
      
      sequenced = (rle(as.character(subbig$sumquantile)))
      
      seqvec = sequenced$values
      names(seqvec) = sequenced$lengths
      seqvec = sort(seqvec)
      numislands = length(seqvec[seqvec == "36"])
      #print(paste(numislands,"islands for",spp))
      
      numislans = length(seqvec[seqvec == "20"])
      print(paste(numislands, "islands,", numislans, "islans for", spp))
      
      sizeislands = (names(seqvec[seqvec == "36"]))
      #hist(as.numeric(sizeislands))
      
      seqtable = t(table(sequenced))
      #barplot(seqtable)
      
      
      
      
      
    }
    
  }
}

## brif islands detour
if(8 %in% steps) {
  print("STEP 8 plot hardcoded islands"); {
    names = c("bel","cri","cur","fla","fus",
              "bil","bru","mel","nit","sin")
    names2 = c("V", "c", "T", "f", "M", "b", "A", "m", "n", "s")
    #sweep = c(25, 49, 14, 92, 212,
    #          257, 72, 123, 72, 78)
    sweep = c(105,172,230,12,6,
              66,27,2,14,3)
    #islan = c(127, 127, 200, 123, 72,
    #          55, 66, 44, 85, 69)
    islan = c(1,1,3,0,0,
              1,2,1,2,6)
    str = c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0)
    
    x = data.frame(names, islan, sweep)
    modelA = glm(str ~ islan, family = "binomial")
    AICcmodavg::AICc(modelA)
    modelB = glm(str ~ sweep, family = "binomial")
    AICcmodavg::AICc(modelB)
    modelC = glm(str ~ islan + sweep, family = "binomial")
    AICcmodavg::AICc(modelC)
    #plot(islan,sweep,pch=names2,col=str+3) # plot with body size on x-axis and survival (0 or 1) on y-axis
    png("sweeps_islands_structure.may2020.png",
        width = 700,
        height = 400)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    par(mar = c(4.5, 4, 0, 0))
    palette(c("red", "yellow"))
    plot(
      islan,
      sweep,
      #pch = str + 16,
      pch=names2,
      col = str + 1,
      cex = 2,
      xlab = "Number Sweeps",
      ylab = "Number Islands"#,
      #xlim=c(0,200),
      #ylim=c(0,270)
    ) # plot with body size on x-axis and survival (0 or 1) on y-axis
    #points(islan, predict(modelC, newdata=x, type="response"),pch=21,bg=as.numeric(str)+2)
    #abline(a = 0,
    #       b = 1,
    #       col = "grey60",
    #       lty = 2,lwd=2)
    legend(
      "topright",
      legend = c("No Structure", "Structure"),
      col = c("red", "yellow"),
      pch = c(16, 17),
      cex = 2,
      border = "white",
      bg = "grey50"
    )
    
    dev.off()
    
    
    y = as.matrix(data.frame(islan, islan))
    rownames(y) = names
    colnames(y) = c("islans", "Islands")
    png("bigplot_regions_per_spp.may2020.png",
        width = 700,
        height = 300)
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    par(mar = c(4.5, 4, 0, 0))
    barplot(
      t(y),
      beside = T,
      col = c("magenta", "cyan"),
      ylab = "Number of Regions",
      xlab = "Species"
    )
    legend(
      "topleft",
      legend = c(rownames(div_swpIsl)),
      fill = c("magenta", "cyan"),
      border = "white",
      bty = "n"
    )
    dev.off()
    
    
    aggmean = aggregate(
      cbind(
        dxymeans,
        stdvs,
        snps,
        sums,
        Nsites,
        Fst,
        dxy_div_fst,
        fst_div_dxy
      ) ~ species,
      data = bigplot,
      FUN = mean
    )
    names(aggmean) = c("species", paste(names(aggmean)[2:length(names(aggmean))], "mean", sep =
                                          "."))
    #aggquan = (aggregate(cbind(dxymeans,stdvs,snps,sums,Nsites,Fst,dxy_div_fst,fst_div_dxy)~species,data=bigplot,FUN=min))
    
    #merged = merge(aggmean,aggquan,by="species")
    #merged = merge(merged,data.frame(islan,islan,str,"species"=names),by="species")
    merged = merge(aggmean, data.frame(sweep, islan, str, "species" = names), by =
                     "species")
    
    
    mean_vals2$species = names
    mean_vals3$species = names
    
    merged = merge(merged, mean_vals2, by = "species")
    merged = merge(merged, mean_vals3, by = "species")
    
    #merged = merged[ , order(names(merged))]
    colms = colnames(merged)
    #merged = matrix(as.numeric(unlist(merged)),nrow=nrow(merged))
    colnames(merged) = colms
    
    corrs = cor(merged[, c(-1,-13,-42)], use = "pairwise.complete.obs")
    corrssq = corrs * corrs
    corrssmall = corrssq[c(1:11), c(12:24, 36:52)]
    
    hist(corrssmall)
    write.table(corrssq, "correlations_rsq.may2020.txt", sep = "\t")
    write.table(corrs, "correlationsfull.may2020.txt", sep = "\t")
    write.table(corrssmall, "correlations_rsqsmall.may2020.txt", sep = "\t")
    
    
    corrplot::corrplot(
      corrssq,
      main = "R squares",
      method = "ellipse",
      is.corr = T,
      tl.cex = 0.5,
      order = "hclust"
    )
    
    
    corrplot::corrplot(
      corrssmall,
      main = "R squares",
      method = "number",
      is.corr = F,
      tl.cex = 0.5,
      number.cex = 0.5
    )
    
    plot(merged[, "PC3"], merged[, "Fst.mean"])
    
    
    
    ## pgls
    library(ape)
    #install_github("bomeara/utilitree")
    library(utilitree)
    library(phylobase)
    library(caper)
    treefile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/subset_from_birdtree_for_pgls/output.nex"
    tree = read.nexus(treefile)
    con = consensus(tree,p=0.5)
    plot(con)
    nodelabels()
    con2 = consensusBrlen(con,tree,type="median_brlen",return.val = 0)
    con4 = consensusBrlen(con,tree,type="mean_brlen",return.val = 0)
    con3 = con
    con3$edge.length = con2[4,]
    plot(con3)
    con3$edge.length = con4[3,]
    plot(con3) ## use this one i guess 
    con3$tip.label = c("bel","mel","bru",
                       "cri","cur","nit",
                       "fus","bil","sin","fla")
    
    #shorebird <- comparative.data(shorebird.tree, shorebird.data, Species, vcv=TRUE, vcv.dim=3)
    birddata = comparative.data(con3,merged,species,vcv=TRUE,vcv.dim=3)
    model1 = lm(BILL.LENGTH ~ sweep,data=merged)
    model2 = pgls(BILL.LENGTH ~ sweep,data=birddata)
    model1
    model2
    summary(model1)$r.squared
    summary(model2)$r.squared
    plot(model1)
    plot(model2)
    rownames(merged) = (merged$species)
    
    library(nlme)
    vf<-diag(vcv(con3))
    
    palette(c("red", "yellow"))
    png("fst_vs_pcs_pgls.may2020.png",width=1200)
    par(mfrow = c(1, 3),
        cex = 2,
        mar = c(4, 4, 1, 1))
    model1 = pgls((PC1 + 6) ~ Fst.mean,data=birddata)
    pval = anova(nlme::gls(PC1+6~ Fst.mean,data=merged,correlation=corPagel(1,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3),"Pval:",round(pval,3))
    plot(
      merged$Fst.mean,
      (merged$PC1 + 6),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model1)
    model2 = pgls((PC2 + 2) ~ Fst.mean,data=birddata)
    pval = anova(nlme::gls(PC2 + 2~ Fst.mean,data=merged,correlation=corPagel(1,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    
    towrite = paste("Rsq:",round(summary(model2)$r.squared,3),"Pval:",round(pval,3))
    plot(
      merged$Fst.mean,
      (merged$PC2 + 2),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model2)
    model3 = pgls((PC3 + 1) ~ Fst.mean,data=birddata)
    pval = anova(nlme::gls(PC3 + 1~ Fst.mean,data=merged,correlation=corPagel(1,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    towrite = paste("Rsq:",round(summary(model3)$r.squared,3),"Pval:",round(pval,3))
    plot(
      merged$Fst.mean,
      (merged$PC3 + 1),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model3)
    dev.off()
    
    png("dxy_vs_pcs_pgls.may2020.png",width=1200)
    par(mfrow = c(1, 3),
        cex = 2,
        mar = c(4, 4, 1, 1))
    model1 = pgls((PC1 + 6) ~ dxymeans.mean,data=birddata)
    towrite = round(summary(model1)$r.squared,3)
    plot(
      merged$dxymeans.mean,
      (merged$PC1 + 6),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model1)
    model2 = pgls((PC2 + 2) ~ dxymeans.mean,data=birddata)
    towrite = round(summary(model2)$r.squared,3)
    plot(
      merged$dxymeans.mean,
      (merged$PC2 + 2),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model2)
    model3 = pgls((PC3 + 1) ~ dxymeans.mean,data=birddata)
    towrite = round(summary(model3)$r.squared,3)
    plot(
      merged$dxymeans.mean,
      (merged$PC3 + 1),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model3)
    dev.off()
    
    png("sweeps_vs_pcs_pgls.may2020.png",
        width = 1200)
    par(mfrow = c(1, 3),
        cex = 2,
        mar = c(4, 4, 1, 1))
    model1 = pgls((PC1 + 6) ~ sweep,data=birddata)
    towrite = round(summary(model1)$r.squared,3)
    plot(
      merged$sweep,
      (merged$PC1 + 6),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model1)
    model2 = pgls((PC2 + 2) ~ sweep,data=birddata)
    towrite = round(summary(model2)$r.squared,3)
    plot(
      merged$sweep,
      (merged$PC2 + 2),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model2)
    model3 = pgls((PC3 + 1) ~ sweep,data=birddata)
    towrite = round(summary(model3)$r.squared,3)
    plot(
      merged$sweep,
      (merged$PC3 + 1),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model3)
    dev.off()
    
    
    png("islands_vs_pcs_pgls.may2020.png",
        width = 1200)
    par(mfrow = c(1, 3),
        cex = 2,
        mar = c(4, 4, 1, 1))
    model1 = pgls((PC1 + 6) ~ islan,data=birddata)
    towrite = round(summary(model1)$r.squared,3)
    plot(
      merged$islan,
      (merged$PC1 + 6),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model1)
    model2 = pgls((PC2 + 2) ~ islan,data=birddata)
    towrite = round(summary(model2)$r.squared,3)
    plot(
      merged$islan,
      (merged$PC2 + 2),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model2)
    model3 = pgls((PC3 + 1) ~ islan,data=birddata)
    towrite = round(summary(model3)$r.squared,3)
    plot(
      merged$islan,
      (merged$PC3 + 1),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model3)
    dev.off()
    
    
    png("str_vs_pcs_pgls.may2020.png",width=1200)
    par(mfrow = c(1, 3),
        cex = 2,
        mar = c(4, 4, 1, 1))
    model1 = pgls((PC1 + 6) ~ str,data=birddata)
    towrite = round(summary(model1)$r.squared,3)
    boxplot(
      (merged$PC1 + 6) ~ merged$str,
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite,
      ylab = "PC1"
    )
    #abline(model1)
    model2 = pgls((PC2 + 2) ~ str,data=birddata)
    towrite = round(summary(model2)$r.squared,3)
    boxplot(
      (merged$PC2 + 2) ~ merged$str,
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite,
      ylab = "PC2"
    )
    #abline(model2)
    model3 = pgls((PC3 + 1) ~ str,data=birddata)
    towrite = round(summary(model3)$r.squared,3)
    boxplot(
      (merged$PC3 + 1) ~ merged$str,
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite,
      ylab = "PC3"
    )
    #abline(model3)
    dev.off()
    
    
    ## top 3 with fst
    ## top 3 with dxy
    ## top 3 with sweep
    ## top 3 with islan
    ## top 3 with str
    
    palette(c("red", "yellow"))
    png("fst_top3_pgls.may2020.png",width=1200)
    par(mfrow = c(1, 3),
        cex = 2,
        mar = c(4, 4, 1, 1))
    model1 = pgls((BILL.LENGTH) ~ Fst.mean,data=birddata)
    
    towrite = round(summary(model1)$r.squared,3)
    plot(
      merged$Fst.mean,
      (merged$BILL.LENGTH),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model1)
    model2 = pgls((TARSUS.LENGTH) ~ Fst.mean,data=birddata)
    
    towrite = round(summary(model2)$r.squared,3)
    plot(
      merged$Fst.mean,
      (merged$TARSUS.LENGTH),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model2)
    model3 = pgls((RES_TAIL + 16) ~ Fst.mean,data=birddata)
    
    towrite = round(summary(model3)$r.squared,3)
    plot(
      merged$Fst.mean,
      (merged$RES_TAIL + 16),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model3)
    dev.off()
    
    png("dxy_top3_pgls.may2020.png",width=1200)
    par(mfrow = c(1, 3),
        cex = 2,
        mar = c(4, 4, 1, 1))
    model1 = pgls((RES_TAIL + 16) ~ dxymeans.mean,data=birddata)
    towrite = round(summary(model1)$r.squared,3)
    plot(
      merged$dxymeans.mean,
      (merged$RES_TAIL + 16),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model1)
    model2 = pgls((RES_WING.LENGTH.TO.PRIMARIES + 27) ~ dxymeans.mean,data=birddata)
    towrite = round(summary(model2)$r.squared,3)
    plot(
      merged$dxymeans.mean,
      (merged$RES_WING.LENGTH.TO.PRIMARIES + 27),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model2)
    model3 = pgls((RES_WING.LENGTH.TO.SECONDARIES + 27) ~ dxymeans.mean,data=birddata)
    towrite = round(summary(model3)$r.squared,3)
    plot(
      merged$dxymeans.mean,
      (merged$RES_WING.LENGTH.TO.SECONDARIES + 27),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model3)
    dev.off()
    
    png("sweeps_top3_pgls.may2020.png",width=1200)
    par(mfrow = c(1, 3),
        cex = 2,
        mar = c(4, 4, 1, 1))
    model1 = pgls((RES_BILL.LENGTH + 8) ~ sweep,data=birddata)
    towrite = round(summary(model1)$r.squared,3)
    plot(
      merged$sweep,
      (merged$RES_BILL.LENGTH + 8),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model1)
    model2 = pgls((BILL.LENGTH) ~ sweep,data=birddata)
    towrite = round(summary(model2)$r.squared,3)
    plot(
      merged$sweep,
      (merged$BILL.LENGTH),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model2)
    model3 = pgls((PC3 + 1) ~ sweep,data=birddata)
    towrite = round(summary(model3)$r.squared,3)
    plot(
      merged$sweep,
      (merged$PC3 + 1),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model3)
    dev.off()
    
    
    png("islands_top3_pgls.may2020.png",width=1200)
    par(mfrow = c(1, 3),
        cex = 2,
        mar = c(4, 4, 1, 1))
    model1 = pgls((RES_BEAKLATERALSURFACE_CALC + 90) ~ islan,data=birddata)
    towrite = round(summary(model1)$r.squared,3)
    plot(
      merged$islan,
      (merged$RES_BEAKLATERALSURFACE_CALC + 90),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model1)
    model2 = pgls((BEAKTOTAREAVOLUME_CALC) ~ islan,data=birddata)
    towrite = round(summary(model2)$r.squared,3)
    plot(
      merged$islan,
      (merged$BEAKTOTAREAVOLUME_CALC),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model2)
    model3 = pgls((RES_BEAKAREAVOLUME_CALC + 1) ~ islan,data=birddata)
    towrite = round(summary(model3)$r.squared,3)
    plot(
      merged$islan,
      ((merged$RES_BEAKAREAVOLUME_CALC + 1)),
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite
    )
    abline(model3)
    dev.off()
    
    
    png("str_top3_pgls.may2020.png",width=1200)
    par(mfrow = c(1, 3),
        cex = 2,
        mar = c(4, 4, 1, 1))
    model1 = pgls((RES_BILL.LENGTH + 8) ~ str,data=birddata)
    towrite = round(summary(model1)$r.squared,3)
    boxplot(
      (merged$RES_BILL.LENGTH + 8) ~ merged$str,
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite,
      ylab = "RES_BILL.LENGTH"
    )
    #abline(model1)
    model2 = pgls((BILL.LENGTH) ~ str,data=birddata)
    towrite = round(summary(model2)$r.squared,3)
    boxplot(
      (merged$BILL.LENGTH) ~ merged$str,
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite,
      ylab = "BILL.LENGTH"
    )
    #abline(model2)
    model3 = pgls((PC3 + 1) ~ str,data=birddata)
    towrite = round(summary(model3)$r.squared,3)
    boxplot(
      (merged$PC3 + 1) ~ merged$str,
      col = as.numeric(merged$str) + 1,
      pch = as.numeric(merged$str) + 16,
      cex = 2,
      main = towrite,
      ylab = "PC3"
    )
    #abline(model3)
    dev.off()
    
  }
}
if(9 %in% steps) {
  print("STEP 9 correlations"); {
    
    head(corrssmall)
    corrssmini = corrssmall[c(1,6,9,10,11),c(11:13)]
    corrplot::corrplot(corrssmini,is.corr=F)
    
    x=data.frame(matrix(nrow=3,ncol=5))
    y=data.frame(matrix(nrow=3,ncol=5))
    rownames(x) = names(merged)[c(24:26)]
    colnames(x) = names(merged)[c(2,7,10:12)]
    rownames(y) = rownames(x)
    colnames(y) = colnames(x)
    
    for (g in c(2,7,10:12)) {
      for (m in c(24:26)) {
        name_g = names(merged)[g]
        name_m = names(merged)[m]
        print(paste(g,name_g,"vs",m,name_m))
        
        model = pgls(birddata$data[,name_m]~birddata$data[,name_g],data=birddata)
        rsq = summary(model)$r.squared  
        pval = summary(model)$coefficients[2,4]
        
        x[name_m,name_g] = rsq
        y[name_m,name_g] = pval
        
      }
    }
    #View(x)
    rsq_cor = as.matrix(x)
    pva_cor = as.matrix(y)
    corrplot::corrplot(rsq_cor,is.corr=F)
    corrplot::corrplot(pva_cor,is.corr=F)
    
    
    write.table(rsq_cor,"rsquared_pgls.may2020.txt")
    write.table(pva_cor,"pval_pgls.may2020.txt")
    
    
    barplot(rsq_cor[,"Fst.mean"])
    
    png("Barplot_morphology_vs_genetics_pcs_pgls.may2020.png")
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white"
    )
    barplot((rsq_cor),beside=T,
            col=c("red","cyan","yellow"),
            xlab=("Genetic Metric"),
            ylab="R-squared value",
            names=c("DXY","FST","#Sweeps",
                    "#Islands",
                    "Structure"),
            ylim=c(0,0.3))
    legend("top",legend=rownames(rsq_cor),
           fill=c("red","cyan","yellow"),
           title="Morphological Metric",
           bty="n")
    dev.off()
    
    palette(c("magenta","brown","red","orange","yellow","goldenrod",
              "darkgreen","green","cyan","lightblue","blue","darkblue",
              "purple"))
    
    png("Barplot_morphology_vs_genetics_raw_pgls.may2020.png")
    # par(
    #   bg = NA,
    #   col.axis = "white",
    #   fg = "white",
    #   col.lab = "white",
    #   col.main = "white"
    # )
    barplot((rsq_cor[c(1:10,25:29),c(1,6,9,10,11)]),beside=T,
            xlab=("Genetic Metric"),
            ylab="R-squared value",ylim=c(0,1),
            col=c("magenta","brown","red","orange","yellow","goldenrod",
                  "darkgreen","green","cyan","lightblue","blue","darkblue",
                  "purple","grey"),
            names=c("DXY","FST","#Sweeps",
                    "#Islands",
                    "Structure"))
    legend("top",ncol=2,
           legend=c(rownames(rsq_cor[c(1:10,25:29),c(1,6,9,10,11)])),
           cex=0.5,fill=c("magenta","brown","red","orange","yellow","goldenrod",
                          "darkgreen","green","cyan","lightblue","blue","darkblue",
                          "purple","grey"),
           title="Morphological Metric",
           bty="n")
    dev.off()
    
    png("Barplot_morphology_vs_genetics_res_pgls.may2020.png")
    # par(
    #   bg = NA,
    #   col.axis = "white",
    #   fg = "white",
    #   col.lab = "white",
    #   col.main = "white"
    # )
    barplot((rsq_cor[c(29:41),c(1,6,9,10,11)]),beside=T,
            xlab=("Genetic Metric"),
            ylab="R-squared value",ylim=c(0,1),
            col=c("magenta","brown","red","orange","yellow","goldenrod",
                  "darkgreen","green","cyan","blue","darkblue",
                  "purple","grey"),
            names=c("DXY","FST","#Sweeps",
                    "#Islands",
                    "Structure"))
    legend("top",ncol=2,
           legend=c(rownames(rsq_cor)[c(29:41)]),
           cex=0.5,fill=c("magenta","brown","red","orange","yellow","goldenrod",
                          "darkgreen","green","cyan","blue","darkblue",
                          "purple","grey"),
           title="Morphological Metric",
           bty="n")
    dev.off()
    
    
    
    
    
    png("Barplot_morphology_vs_fst_all_pgls.may2020.png")
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white",
      mar=c(4,10,0,1)
    )
    bp = barplot((rsq_cor[c(1:13,25:41),"Fst.mean"]),beside=T,
                 xlab="R-squared value with Fst",xlim=c(0,1),
                 ylab="",
                 col=c("grey"),las=2,cex.names=0.5,
                 horiz=T
                 #names=c(colnames(corrssmall)[c(18:30)]),
                 #names=rep("",13)
    )
    #plot(bp)
    # text(bp, c(0.05,0.175,0.05,0.15,0.05,
    #            0.05,0.2,0.05,0.075,0.075,0.075,0.075,0.075), c("Bill Height","Bill Length",
    #               "Bill Width","Tarsus Length",
    #               "Primaries","Secondaries","Tail",
    #               "Kipps Index","Beak Base Area",
    #               "Beak Volume",
    #               "Beak Lateral Surface Area",
    #               "Beak Total Surface Area",
    #               "Beak Surface Volume Ratio"),
    #      cex=0.5,pos=3,srt=90) 
    # mtext("Residual Morphological Variable",side=1)
    dev.off()
    
    
    
    png("Barplot_morphology_vs_dxy_all_pgls.may2020.png")
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white",
      mar=c(4,10,0,1)
    )
    bp = barplot((rsq_cor[c(1:13,25:41),"dxymeans.mean"]),beside=T,
                 xlab="R-squared value with Dxy",xlim=c(0,1),
                 ylab="",
                 col=c("grey"),las=2,cex.names=0.5,
                 horiz=T
                 #names=c(colnames(corrssmall)[c(18:30)]),
                 #names=rep("",13)
    )
    #plot(bp)
    # text(bp, c(0.05,0.175,0.05,0.15,0.05,
    #            0.05,0.2,0.05,0.075,0.075,0.075,0.075,0.075), c("Bill Height","Bill Length",
    #               "Bill Width","Tarsus Length",
    #               "Primaries","Secondaries","Tail",
    #               "Kipps Index","Beak Base Area",
    #               "Beak Volume",
    #               "Beak Lateral Surface Area",
    #               "Beak Total Surface Area",
    #               "Beak Surface Volume Ratio"),
    #      cex=0.5,pos=3,srt=90) 
    # mtext("Residual Morphological Variable",side=1)
    dev.off()
    
    
    png("Barplot_morphology_vs_sweep_all_pgls.may2020.png")
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white",
      mar=c(4,10,0,1)
    )
    bp = barplot((rsq_cor[c(1:13,25:41),"sweep"]),beside=T,
                 xlab="R-squared value with #Sweeps",xlim=c(0,1),
                 ylab="",
                 col=c("grey"),las=2,cex.names=0.5,
                 horiz=T
                 #names=c(colnames(corrssmall)[c(18:30)]),
                 #names=rep("",13)
    )
    #plot(bp)
    # text(bp, c(0.05,0.175,0.05,0.15,0.05,
    #            0.05,0.2,0.05,0.075,0.075,0.075,0.075,0.075), c("Bill Height","Bill Length",
    #               "Bill Width","Tarsus Length",
    #               "Primaries","Secondaries","Tail",
    #               "Kipps Index","Beak Base Area",
    #               "Beak Volume",
    #               "Beak Lateral Surface Area",
    #               "Beak Total Surface Area",
    #               "Beak Surface Volume Ratio"),
    #      cex=0.5,pos=3,srt=90) 
    # mtext("Residual Morphological Variable",side=1)
    dev.off()
    
    
    
    png("Barplot_morphology_vs_islands_all_pgls.may2020.png")
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white",
      mar=c(4,10,0,1)
    )
    bp = barplot((rsq_cor[c(1:13,25:41),"islan"]),beside=T,
                 xlab="R-squared value with #Islands",xlim=c(0,1),
                 ylab="",
                 col=c("grey"),las=2,cex.names=0.5,
                 horiz=T
                 #names=c(colnames(corrssmall)[c(18:30)]),
                 #names=rep("",13)
    )
    #plot(bp)
    # text(bp, c(0.05,0.175,0.05,0.15,0.05,
    #            0.05,0.2,0.05,0.075,0.075,0.075,0.075,0.075), c("Bill Height","Bill Length",
    #               "Bill Width","Tarsus Length",
    #               "Primaries","Secondaries","Tail",
    #               "Kipps Index","Beak Base Area",
    #               "Beak Volume",
    #               "Beak Lateral Surface Area",
    #               "Beak Total Surface Area",
    #               "Beak Surface Volume Ratio"),
    #      cex=0.5,pos=3,srt=90) 
    # mtext("Residual Morphological Variable",side=1)
    dev.off()
    
    
    png("Barplot_morphology_vs_str_all_pgls.may2020.png")
    par(
      bg = NA,
      col.axis = "white",
      fg = "white",
      col.lab = "white",
      col.main = "white",
      mar=c(4,10,0,1)
    )
    bp = barplot((rsq_cor[c(1:13,25:41),"str"]),beside=T,
                 xlab="R-squared value with Structure",xlim=c(0,1),
                 ylab="",
                 col=c("grey"),las=2,cex.names=0.5,
                 horiz=T
                 #names=c(colnames(corrssmall)[c(18:30)]),
                 #names=rep("",13)
    )
    #plot(bp)
    # text(bp, c(0.05,0.175,0.05,0.15,0.05,
    #            0.05,0.2,0.05,0.075,0.075,0.075,0.075,0.075), c("Bill Height","Bill Length",
    #               "Bill Width","Tarsus Length",
    #               "Primaries","Secondaries","Tail",
    #               "Kipps Index","Beak Base Area",
    #               "Beak Volume",
    #               "Beak Lateral Surface Area",
    #               "Beak Total Surface Area",
    #               "Beak Surface Volume Ratio"),
    #      cex=0.5,pos=3,srt=90) 
    # mtext("Residual Morphological Variable",side=1)
    dev.off()
    
    
    boxplot(merged$KIPPSINDEX_CALC~merged$str)
    
    
    ## just get the pgls pvals
    anova(nlme::gls(PC1~ Fst.mean,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    anova(nlme::gls(PC2~ Fst.mean,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    anova(nlme::gls(PC3~ Fst.mean,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    
    anova(nlme::gls(PC1~ dxymeans.mean,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    anova(nlme::gls(PC2~ dxymeans.mean,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    ## sig 
    anova(nlme::gls(PC3~ dxymeans.mean,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    
    anova(nlme::gls(PC1~ sweep,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    anova(nlme::gls(PC2~ sweep,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    ## sig 
    anova(nlme::gls(PC3~ sweep,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    ## sig 
    anova(nlme::gls(PC1~ islan,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    anova(nlme::gls(PC2~ islan,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    anova(nlme::gls(PC3~ islan,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    
    anova(nlme::gls(PC1~ str,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    anova(nlme::gls(PC2~ str,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    anova(nlme::gls(PC3~ str,data=merged,correlation=corPagel(0.5,con3),method="REML",weights=varFixed(~vf)))$`p-value`[2]
    ## sig 
    
    
    ## check if they match difs in morph
    absmatrix
    names(absmatrix)
    
    popmerge = merge(absmatrix,aggmean,by="species")
    popmerge = merge(popmerge,data.frame(sweep, islan, str, "species" = names),by="species")
    
    plot(popmerge[,c(2:ncol(popmerge))])
    
    
    birddata = comparative.data(con3,popmerge,species,vcv=TRUE,vcv.dim=3)
    
    x=data.frame(matrix(nrow=3,ncol=5))
    y=data.frame(matrix(nrow=3,ncol=5))
    rownames(x) = names(popmerge)[c(2:4)]
    colnames(x) = names(popmerge)[c(5,10,13:15)]
    rownames(y) = rownames(x)
    colnames(y) = colnames(x)
    
    for (g in c(5,10,13:15)) {
      for (m in c(2:4)) {
        name_g = names(popmerge)[g]
        name_m = names(popmerge)[m]
        print(paste(g,name_g,"vs",m,name_m))
        
        model = pgls(birddata$data[,name_m]~birddata$data[,name_g],data=birddata)
        rsq = summary(model)$r.squared  
        pval = summary(model)$coefficients[2,4]
        
        x[name_m,name_g] = rsq
        y[name_m,name_g] = pval
        
      }
    }
    #View(x)
    rsq_des = as.matrix(x)
    pva_des = as.matrix(y)
    corrplot::corrplot(rsq_des,is.corr=F)
    corrplot::corrplot(pva_des,is.corr=F)
    
    
    write.table(rsq_des,"rsquared_pgls_des.may2020.txt")
    write.table(pva_des,"pval_pgls_des.txt.may2020")
    
    
    png("Barplot_despca_vs_genomics_pgls.may2020.png")
    #par(  bg = NA,  col.axis = "white",  fg = "white",  col.lab = "white",  col.main = "white")
    barplot((rsq_des),beside=T,       xlab=("Genetic Metric"),
            ylab="R-squared value",col=c("red","cyan","yellow"),
            names=c("DXY","FST","#Sweep","#Island","Structure"),
            ylim=c(0,0.3))
    legend("top",legend=rownames(rsq_cor),
           fill=c("red","cyan","yellow"),
           title="Morphological Metric",
           bty="n")
    dev.off()
    
    png("correlations_populationlevel_genomemorpho_pc1dif.may2020.png",width=1000,height=300)
    #par(  bg = NA,  col.axis = "white",  fg = "white",  col.lab = "white",  col.main = "white")
    par(mfrow=c(1,5))
    model1 = pgls(birddata$data[,"pc1dif"]~birddata$data[,"dxymeans.mean"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    plot(popmerge[,c(5,2)], xlab="DXY",ylab="Dif pc1dif between pops",col=birddata$data[,"str"]+1,pch=birddata$data[,"str"]+16,cex=2,main=towrite)
    abline(model1)
    model1 = pgls(birddata$data[,"pc1dif"]~birddata$data[,"Fst.mean"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    plot(popmerge[,c(10,2)],xlab="FST",ylab="Dif pc1dif between pops",col=birddata$data[,"str"]+1,pch=birddata$data[,"str"]+16,cex=2,main=towrite)
    abline(model1)
    model1 = pgls(birddata$data[,"pc1dif"]~birddata$data[,"sweep"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    plot(popmerge[,c(13,2)],xlab="Sweeps",ylab="Dif pc1dif between pops",col=birddata$data[,"str"]+1,pch=birddata$data[,"str"]+16,cex=2,main=towrite)
    abline(model1)
    model1 = pgls(birddata$data[,"pc1dif"]~birddata$data[,"islan"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    plot(popmerge[,c(14,2)],xlab="Islands",ylab="Dif pc1dif between pops",col=birddata$data[,"str"]+1,pch=birddata$data[,"str"]+16,cex=2,main=towrite)
    abline(model1)
    model1 = pgls(birddata$data[,"pc1dif"]~birddata$data[,"str"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    boxplot(birddata$data[,"pc1dif"]~birddata$data[,"str"],xlab="Structure",ylab="Dif pc1dif between pops",col=c("red","cyan"),cex=2,main=towrite)
    dev.off()
    
    png("correlations_populationlevel_genomemorpho_pc2dif.may2020.png",width=1000,height=300)
    #par(  bg = NA,  col.axis = "white",  fg = "white",  col.lab = "white",  col.main = "white")
    par(mfrow=c(1,5))
    model1 = pgls(birddata$data[,"pc2dif"]~birddata$data[,"dxymeans.mean"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    plot(popmerge[,c(5,3)], xlab="DXY",ylab="Dif pc2dif between pops",col=birddata$data[,"str"]+1,pch=birddata$data[,"str"]+16,cex=2,main=towrite)
    abline(model1)
    model1 = pgls(birddata$data[,"pc2dif"]~birddata$data[,"Fst.mean"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    plot(popmerge[,c(10,3)],xlab="FST",ylab="Dif pc2dif between pops",col=birddata$data[,"str"]+1,pch=birddata$data[,"str"]+16,cex=2,main=towrite)
    abline(model1)
    model1 = pgls(birddata$data[,"pc2dif"]~birddata$data[,"sweep"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    plot(popmerge[,c(13,3)],xlab="Sweeps",ylab="Dif pc2dif between pops",col=birddata$data[,"str"]+1,pch=birddata$data[,"str"]+16,cex=2,main=towrite)
    abline(model1)
    model1 = pgls(birddata$data[,"pc2dif"]~birddata$data[,"islan"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    plot(popmerge[,c(14,3)],xlab="Islands",ylab="Dif pc2dif between pops",col=birddata$data[,"str"]+1,pch=birddata$data[,"str"]+16,cex=2,main=towrite)
    abline(model1)
    model1 = pgls(birddata$data[,"pc2dif"]~birddata$data[,"str"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    boxplot(birddata$data[,"pc2dif"]~birddata$data[,"str"],xlab="Structure",ylab="Dif pc2dif between pops",col=c("red","cyan"),cex=2,main=towrite)
    dev.off()
    
    png("correlations_populationlevel_genomemorpho_pc3dif.may2020.png",width=1000,height=300)
    #par(  bg = NA,  col.axis = "white",  fg = "white",  col.lab = "white",  col.main = "white")
    par(mfrow=c(1,5))
    model1 = pgls(birddata$data[,"pc3dif"]~birddata$data[,"dxymeans.mean"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    plot(popmerge[,c(5,4)], xlab="DXY",ylab="Dif pc3dif between pops",col=birddata$data[,"str"]+1,pch=birddata$data[,"str"]+16,cex=2,main=towrite)
    abline(model1)
    model1 = pgls(birddata$data[,"pc3dif"]~birddata$data[,"Fst.mean"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    plot(popmerge[,c(10,4)],xlab="FST",ylab="Dif pc3dif between pops",col=birddata$data[,"str"]+1,pch=birddata$data[,"str"]+16,cex=2,main=towrite)
    abline(model1)
    model1 = pgls(birddata$data[,"pc3dif"]~birddata$data[,"sweep"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    plot(popmerge[,c(13,4)],xlab="Sweeps",ylab="Dif pc3dif between pops",col=birddata$data[,"str"]+1,pch=birddata$data[,"str"]+16,cex=2,main=towrite)
    abline(model1)
    model1 = pgls(birddata$data[,"pc3dif"]~birddata$data[,"islan"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    plot(popmerge[,c(14,4)],xlab="Islands",ylab="Dif pc3dif between pops",col=birddata$data[,"str"]+1,pch=birddata$data[,"str"]+16,cex=2,main=towrite)
    abline(model1)
    model1 = pgls(birddata$data[,"pc3dif"]~birddata$data[,"str"],data=birddata)
    towrite = paste("Rsq:",round(summary(model1)$r.squared,3))
    boxplot(birddata$data[,"pc3dif"]~birddata$data[,"str"],xlab="Structure",ylab="Dif pc3dif between pops",col=c("red","cyan"),cex=2,main=towrite)
    dev.off()
    
    ## do the randomization tests
    
    smallmerged = merged[,c(2,7,10:12,24:26,1)]
    pcs = smallmerged[,c(6:8)]
    genes = smallmerged[,c(1:5,9)]
    
    dfs = lapply(1:10000,FUN=function(x) {pcs[sample(1:nrow(pcs)), ]})
    rsqs = lapply(1:length(dfs),FUN=function(i){
      print(i)
      df = dfs[[i]]
      dfmerge = cbind(genes,df)
      birddata = comparative.data(con3,dfmerge,species,vcv=TRUE,vcv.dim=3)
      x=data.frame(matrix(nrow=3,ncol=5))
      rownames(x) = names(dfmerge)[c(7:9)]
      colnames(x) = names(dfmerge)[c(1:5)]
      for (g in c(1:5)) {
        for (m in c(7:9)) {
          name_g = names(dfmerge)[g]
          name_m = names(dfmerge)[m]
          model = pgls(birddata$data[,name_m]~birddata$data[,name_g],data=birddata)
          rsq = summary(model)$r.squared  
          x[name_m,name_g] = rsq
        }
      }
      rsq_des = as.matrix(x)
      return(rsq_des)
    })
    
    one1 = c()
    one2 = c()
    one3 = c()
    one4 = c()
    one5 = c()
    two1 = c()
    two2 = c()
    two3 = c()
    two4 = c()
    two5 = c()
    three1 = c()
    three2 = c()
    three3 = c()
    three4 = c()
    three5 = c()
    
    
    for (i in 1:length(rsqs)) {
      todo = rsqs[[i]]
      one1 = c(one1,todo[1,1])
      one2 = c(one2,todo[1,2])
      one3 = c(one3,todo[1,3])
      one4 = c(one4,todo[1,4])
      one5 = c(one5,todo[1,5])
      
      two1 = c(two1,todo[2,1])
      two2 = c(two2,todo[2,2])
      two3 = c(two3,todo[2,3])
      two4 = c(two4,todo[2,4])
      two5 = c(two5,todo[2,5])
      
      three1 = c(three1,todo[3,1])
      three2 = c(three2,todo[3,2])
      three3 = c(three3,todo[3,3])
      three4 = c(three4,todo[3,4])
      three5 = c(three5,todo[3,5])
      
    }
    
    matrix = cbind(one1,one2,one3,one4,one5,
                   two1,two2,two3,two4,two5,
                   three1,three2,three3,three4,three5)
    
    colnames(matrix) = c("dxy-pc1","fst-pc1","sweep-pc1","islan-pc1","str-pc1",
                         "dxy-pc2","fst-pc2","sweep-pc2","islan-pc2","str-pc2",
                         "dxy-pc3","fst-pc3","sweep-pc3","islan-pc3","str-pc3")
    
    rsq_cor
    rsq_cor_line = c(rsq_cor[1,],rsq_cor[2,],rsq_cor[3,])
    names(rsq_cor_line) = colnames(matrix)
    
    png("Randomization_Tests_For_Across_Spp_PGLS.may2020.png",
        height=600,width=900)
    par(mfrow=c(3,5))
    for (i in 1:15) {
      data = matrix[,i]
      hist(data,breaks=50,xlab="R-squared",
           main=colnames(matrix)[i],
           xlim=c(0,1))
      emp = rsq_cor_line[i]
      abline(v=emp,col="red",lty=2,lwd=2)
      
      quan = sum(data < emp)/length(data)
      #print(quan)
      legend("topright",legend=quan,title="Quantile",
             bty="n")
      
    }
    dev.off()
    
    
    
    
    ## overlay the fst plots -- see if things line up
    ##Sturnus -- Toxostoma, Campylorhynchus, Polioptila, Phainopepla
    ##Zonotrichia -- Amphispiza, Pipilo
    
    zonbig = bigplot[bigplot$species %in% c("bil","fus"),]
    zonsmall = zonbig[,c(1,2,4,10,15,17)]
    zonsmall1 = zonsmall[zonsmall$species=="bil",]
    zonsmall2 = zonsmall[zonsmall$species=="fus",]
    zono = merge(zonsmall1,zonsmall2,by=c("chr","midPos"))
    
    zonomeansdxy = rowMeans(zono[,c(3,7)])
    zonomeansfst = rowMeans(zono[,c(4,8)])
    zonomeansqua = rowMeans(zono[,c(6,10)])
    zono = cbind(zono,zonomeansdxy,zonomeansfst,zonomeansqua)
    
    png("same_zono_fst.may2020.png",width=800,height=500)
    par(mfrow=c(2,1),mar=c(4,4,0,0))
    plot(zono[,4],col=as.numeric(as.factor(zono$chr)),xlab=unique(zono[,5]),ylab="Fst")
    plot(zono[,8],col=as.numeric(as.factor(zono$chr)),xlab=unique(zono[,9]),ylab="Fst")
    dev.off()
    
    png("same_zono_dxy.may2020.png",width=800,height=500)
    par(mfrow=c(2,1),mar=c(4,4,0,0))
    plot(zono[,3],col=as.numeric(as.factor(zono$chr)),xlab=unique(zono[,5]),ylab="Dxy")
    plot(zono[,7],col=as.numeric(as.factor(zono$chr)),xlab=unique(zono[,9]),ylab="Dxy")
    dev.off()
    
    colfunc <- colorRampPalette(c("red", "black","cyan"))
    palette(colfunc(17))
    png("corr_zono_dxy.may2020.png",width=500,height=500)
    plot(zono[,3],
         zono[,7],
         xlab=paste(unique(zono[,5]),names(zono)[3]),
         ylab=paste(unique(zono[,9]),names(zono)[7]),
         pch=16,
         col=as.numeric((zono$quancolor.x-zono$quancolor.y)+9))
    dev.off()
    
    png("corr_zono_fst.may2020.png",width=500,height=500)
    plot(zono[,4],
         zono[,8],
         xlab=paste(unique(zono[,5]),names(zono)[4]),
         ylab=paste(unique(zono[,9]),names(zono)[8]),
         pch=16,
         col=as.numeric((zono$quancolor.x-zono$quancolor.y)+9))
    dev.off()
    
    
    png("same_zono_dxyfstmeans.may2020.png",width=800,height=500)
    par(mfrow=c(2,1),mar=c(4,4,0,0))
    plot(zono$zonomeansfst,col=as.numeric(as.factor(zono$chr)),ylab="Fst mean across spp")
    plot(zono$zonomeansdxy,col=as.numeric(as.factor(zono$chr)),ylab="Dxy mean across spp")
    dev.off()
    
    png("same_zono_quantilemeansDXY.may2020.png",width=800,height=500)
    plot(zono$zonomeansqua,col=as.numeric(as.factor(zono$chr)),ylab="Quantile DXY mean across spp")
    dev.off()
    
    
    
    
    sturbig = bigplot[bigplot$species %in% c("bru","cri","cur","nit","mel"),]
    stursmall = sturbig[,c(1,2,4,10,15,18)]
    stursmall1 = stursmall[stursmall$species=="bru",]
    stursmall2 = stursmall[stursmall$species=="cri",]
    stur = merge(stursmall1,stursmall2,by=c("chr","midPos")) 
    stursmall3 = stursmall[stursmall$species=="cur",]
    stur = merge(stur,stursmall3,by=c("chr","midPos")) 
    stursmall4 = stursmall[stursmall$species=="nit",]
    stur = merge(stur,stursmall4,by=c("chr","midPos")) 
    stursmall5 = stursmall[stursmall$species=="mel",]
    stur = merge(stur,stursmall5,by=c("chr","midPos")) 
    
    
    png("same_stur_fst.may2020.png",width=800,height=1250)
    par(mfrow=c(5,1),mar=c(4,4,0,0))
    plot(stur[,4],col=as.numeric(as.factor(stur$chr)),xlab=unique(stur[,5]),ylab="Fst")
    plot(stur[,8],col=as.numeric(as.factor(stur$chr)),xlab=unique(stur[,9]),ylab="Fst")
    plot(stur[,12],col=as.numeric(as.factor(stur$chr)),xlab=unique(stur[,13]),ylab="Fst")
    plot(stur[,16],col=as.numeric(as.factor(stur$chr)),xlab=unique(stur[,17]),ylab="Fst")
    plot(stur[,20],col=as.numeric(as.factor(stur$chr)),xlab=unique(stur[,21]),ylab="Fst")
    dev.off()
    
    png("same_stur_dxy.may2020.png",width=800,height=1250)
    par(mfrow=c(5,1),mar=c(4,4,0,0))
    plot(stur[,3],col=as.numeric(as.factor(stur$chr)),xlab=unique(stur[,5]),ylab="Dxy")
    plot(stur[,7],col=as.numeric(as.factor(stur$chr)),xlab=unique(stur[,9]),ylab="Dxy")
    plot(stur[,11],col=as.numeric(as.factor(stur$chr)),xlab=unique(stur[,13]),ylab="Dxy")
    plot(stur[,15],col=as.numeric(as.factor(stur$chr)),xlab=unique(stur[,17]),ylab="Dxy")
    plot(stur[,19],col=as.numeric(as.factor(stur$chr)),xlab=unique(stur[,21]),ylab="Dxy")
    dev.off()
    
    corrs = cor(stur[,c(3,7,11,15,19,4,8,12,16,20)])
    corrplot::corrplot(corrs,method="number",type="upper",diag=F)
    
    ## calculate a mean value for each row
    sturmeansdxy = rowMeans(stur[,c(3,7,11,15,19)])
    sturmeansfst = rowMeans(stur[,c(4,8,12,16,20)])
    sturmeansqua = rowMeans(stur[,c(6,10,14,18,22)])
    stur = cbind(stur,sturmeansdxy,sturmeansfst,sturmeansqua)
    
    png("same_stur_dxyfstmeans.may2020.png",width=800,height=500)
    par(mfrow=c(2,1),mar=c(4,4,0,0))
    plot(stur$sturmeansfst,col=as.numeric(as.factor(stur$chr)),ylab="Fst mean across spp")
    plot(stur$sturmeansdxy,col=as.numeric(as.factor(stur$chr)),ylab="Dxy mean across spp")
    dev.off()
    
    png("same_stur_quantilemeans.may2020.png",width=800,height=500)
    plot(stur$sturmeansqua,col=as.numeric(as.factor(stur$chr)),ylab="Quantile mean across spp")
    dev.off()
    
    palette(    c(      "white",      "white",
                        "white",      "white",      "white",      "magenta",
                        "white",      "white",      "cyan"    )  )
    png("same_stur_quancol.may2020.png",width=800,height=500)
    #par(mfrow=c(5,1),mar=c(4,4,0,0))
    plot(log(stur[,6]+2),col=as.numeric(as.factor(stur[,6])),xlab=unique(stur[,5]),ylab="quantile",ylim=c(log(20-2),log(36+2)))
    points(log(stur[,10]+1),col=as.numeric(as.factor(stur[,10])),xlab=unique(stur[,9]),ylab="quantile",ylim=c(log(20-2),log(36+2)))
    points(log(stur[,14]),col=as.numeric(as.factor(stur[,14])),xlab=unique(stur[,13]),ylab="quantile",ylim=c(log(20-2),log(36+2)))
    points(log(stur[,18]-1),col=as.numeric(as.factor(stur[,18])),xlab=unique(stur[,17]),ylab="quantile",ylim=c(log(20-2),log(36+2)))
    points(log(stur[,22]-2),col=as.numeric(as.factor(stur[,22])),xlab=unique(stur[,21]),ylab="quantile",ylim=c(log(20-2),log(36+2)))
    dev.off()
    
    
    ## now are the sweeps and islands in the same spot? 
    ## sweep = 20
    ## islands = 36
    
    
  }
}
if(10 %in% steps) {
  print("STEP 10 -- 5 ST DEVS"); {
    
    #bigplot = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.may2020.txt")
    names(bigplot)
    
    var1num=which(names(bigplot)=="dxymeans")
    var2num=which(names(bigplot)=="Fst")
    varcol=which(names(bigplot)=="quancolor")
    
    for (spp in specieslist) {
      print(spp)
      
      makedxyplot(speciesname=spp,variable1num=var1num,variable2num=var2num,colvarnum=varcol,dat=bigplot,plotsize="SMALL",
                  lowquan1=0.05,highquan1=0.95,lowquan2=0.05,highquan2=0.95)
      for (i in 1:5) {
        print(i)
        makedxyplotSTDEV(speciesname=spp,variable1num=var1num,variable2num=var2num,colvarnum=varcol,dat=bigplot,plotsize="SMALL",numsd=i)
      }  
    }
    
  }
}
if(11 %in% steps) {
  print("STEP 11 -- DXYPLOTS"); {
    for (spp in specieslist) {
      
      var1num=which(names(bigplot)=="dxymeans")
      var2num=which(names(bigplot)=="Fst")
      varcol=which(names(bigplot)=="quancolor")
      
      print(spp)
      makedxyplot(speciesname=spp,variable1num=var1num,variable2num=var2num,colvarnum=varcol,dat=bigplot,plotsize="LARGE")
      makedxyplot(speciesname=spp,variable1num=var1num,variable2num=var2num,colvarnum=varcol,dat=bigplot,plotsize="SMALL")
    }
  }
}

## overlaying bigplots over one another
if(12 %in% steps) {
bigplot=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/bigplot.may2020.txt",
                   sep=",",header=T)
bigplot = bigplot[order(bigplot$chr,bigplot$starts),]
bigplot$chrpos = paste(bigplot$chr,bigplot$starts)
bigplot$chrpos = as.numeric(as.factor(bigplot$chrpos))
png("testbigplotoverlay.may2020.png",height=1000,width=800)
par(mfrow=c(10,1),mar=c(0,3,0,0))
for (i in 1:length(unique(bigplot$species))){
  spp=sort(unique(bigplot$species))[i]
  temp = bigplot[bigplot$species==spp,]
  temp$zfst = (temp$Fst-mean(temp$Fst,na.rm=T))/sd(temp$Fst,na.rm=T)
  print(max(temp$zfst,na.rm=T))
  print(min(temp$zfst,na.rm=T))
  plot(temp$chrpos,temp$zfst,col=as.numeric(as.factor(temp$chr)),pch=16,ylim=c(-15,75),xlim=c(0,90000),
       cex=0.5)
}
dev.off()

bigplot = bigplot[order(bigplot$chr,bigplot$starts),]
bigplot$chrpos = paste(bigplot$chr,bigplot$starts)
bigplot$chrpos = as.numeric(as.factor(bigplot$chrpos))
png("testbigplotoverlaydxy.may2020.png",height=1000,width=800)
par(mfrow=c(10,1),mar=c(0,3,0,0))
for (i in 1:length(unique(bigplot$species))){
  spp=sort(unique(bigplot$species))[i]
  temp = bigplot[bigplot$species==spp,]
  temp$zfst = (temp$dxymeans-mean(temp$dxymeans,na.rm=T))/sd(temp$dxymeans,na.rm=T)
  print(max(temp$zfst,na.rm=T))
  print(min(temp$zfst,na.rm=T))
  plot(temp$chrpos,temp$zfst,col=as.numeric(as.factor(temp$chr)),pch=16,ylim=c(-15,75),xlim=c(0,90000),
       cex=0.5)
}
dev.off()

## combine by species  -- fst, dxy, nsites
bigplot = bigplot[order(bigplot$chr,bigplot$starts),]
bigplot$chrpos = paste(bigplot$chr,bigplot$starts)
bigplot$chrpos = as.numeric(as.factor(bigplot$chrpos))
btemp=c()
for (i in 1:length(unique(bigplot$species))){
  spp=sort(unique(bigplot$species))[i]
  temp = bigplot[bigplot$species==spp,]
  temp=temp[,c("chr","starts","dxymeans","Fst","Nsites","chrpos")]
  colnames(temp)[3:5] = paste(colnames(temp)[3:5],spp,sep=".")
  if(is.null(btemp)){
    btemp=temp
  } else {
    btemp=merge(btemp,temp,by=c("chr","starts","chrpos"),all=T)
  }
}

dim(btemp)

dxys=c(seq(4,33,by=3))
fsts=c(seq(5,33,by=3))
nsit=c(seq(6,33,by=3))

dcor=cor(btemp[,dxys],use="pairwise.complete.obs")
fcor=cor(btemp[,fsts],use="pairwise.complete.obs")
ncor=cor(btemp[,nsit],use="pairwise.complete.obs")

listspp=unlist(strsplit(colnames(dcor),"\\."))[seq(2,length(colnames(dcor))*2,2)]
colnames(dcor)=listspp
rownames(dcor)=listspp
listspp=unlist(strsplit(colnames(fcor),"\\."))[seq(2,length(colnames(fcor))*2,2)]
colnames(fcor)=listspp
rownames(fcor)=listspp
listspp=unlist(strsplit(colnames(ncor),"\\."))[seq(2,length(colnames(ncor))*2,2)]
colnames(ncor)=listspp
rownames(ncor)=listspp

pdf("dxy_fst_nsites_corr_fullgenome.may2020.pdf",height=1.5,width=3)
par(mfrow=c(1,3))
corrplot::corrplot(dcor,
                   method="ellipse",diag=F,order="alphabet",
                   tl.cex=0.5,cl.cex=0.5,number.cex=0.5,pch.cex=0.5,
                   title="\ndxy",
                   type="upper",
                   mar=c(0,0,0,0))
corrplot::corrplot(fcor,
                   method="ellipse",diag=F,order="alphabet",
                   tl.cex=0.5,cl.cex=0.5,number.cex=0.5,pch.cex=0.5,
                   title="\nfsts",
                   type="upper",
                   mar=c(0,0,0,0))
corrplot::corrplot(ncor,
                   method="ellipse",diag=F,order="alphabet",
                   tl.cex=0.5,cl.cex=0.5,number.cex=0.5,pch.cex=0.5,
                   title="\nnsites",
                   type="upper",
                   mar=c(0,0,0,0))
dev.off()


corrplot::corrplot(cor(btemp[,dxys],use="pairwise.complete.obs"),
                   method="number",diag=T,order="hclust")
## sin and bel dxys are 100% correlated -- probably wrong
corrplot::corrplot(cor(btemp[,fsts],use="pairwise.complete.obs"),
                   method="number",diag=T,order="hclust")
corrplot::corrplot(cor(btemp[,nsit],use="pairwise.complete.obs"),
                   method="number",diag=T,order="hclust")


pdf("dxy_corrs_chroms.may2020.pdf",height=5,width=8)
par(mfrow=c(4,7))
for(i in 1:length(unique(btemp$chr))) {
  print(i)
  thischr=unique(btemp$chr)[i]
  thisset=btemp[btemp$chr==thischr,]
  thiscorr=cor(thisset[,dxys],use="pairwise.complete.obs")
  if(sum(is.na(thiscorr)) >= 100){
    print("BAD")
    frame()
    text(0.5,0.5,paste("\ndxy chr",thischr))
  } else {
    
    blank=which(rowSums(is.na(thiscorr))==10)
    if (length(blank) > 0){
      thiscorr=thiscorr[-blank,-blank]
    }
    
    listspp=unlist(strsplit(colnames(thiscorr),"\\."))[seq(2,length(colnames(thiscorr))*2,2)]
    colnames(thiscorr)=listspp
    rownames(thiscorr)=listspp
    
  corrplot::corrplot(thiscorr,
                     method="ellipse",diag=F,order="alphabet",
                     tl.cex=0.5,cl.cex=0.5,number.cex=0.5,pch.cex=0.5,
                     title=paste("\ndxy chr",thischr),
                     type="upper",
                     mar=c(0,0,0,0))
  }
}
dev.off()

pdf("fst_corrs_chroms.may2020.pdf",height=5,width=8)
par(mfrow=c(4,7))
for(i in 1:length(unique(btemp$chr))) {
  print(i)
  thischr=unique(btemp$chr)[i]
  thisset=btemp[btemp$chr==thischr,]
  thiscorr=cor(thisset[,fsts],use="pairwise.complete.obs")
  if(sum(is.na(thiscorr)) >= 100){
    print("BAD")
    frame()
    text(0.5,0.5,paste("\nfst chr",thischr))
  } else {
    
    blank=which(rowSums(is.na(thiscorr))==10)
    if (length(blank) > 0){
      thiscorr=thiscorr[-blank,-blank]
    }
    
    listspp=unlist(strsplit(colnames(thiscorr),"\\."))[seq(2,length(colnames(thiscorr))*2,2)]
    colnames(thiscorr)=listspp
    rownames(thiscorr)=listspp
    
    corrplot::corrplot(thiscorr,
                       method="ellipse",diag=F,order="alphabet",
                       tl.cex=0.5,cl.cex=0.5,number.cex=0.5,pch.cex=0.5,
                       title=paste("\nfst chr",thischr),
                       type="upper",
                       mar=c(0,0,0,0))
  }
}
dev.off()

pdf("nsites_corrs_chroms.may2020.pdf",height=5,width=8)
par(mfrow=c(4,7))
for(i in 1:length(unique(btemp$chr))) {
  print(i)
  thischr=unique(btemp$chr)[i]
  thisset=btemp[btemp$chr==thischr,]
  thiscorr=cor(thisset[,nsit],use="pairwise.complete.obs")
  if(sum(is.na(thiscorr)) >= 100){
    print("BAD")
    frame()
    text(0.5,0.5,paste("\nnsites chr",thischr))
  } else {
    
    blank=which(rowSums(is.na(thiscorr))==10)
    if (length(blank) > 0){
      thiscorr=thiscorr[-blank,-blank]
    }
    
    listspp=unlist(strsplit(colnames(thiscorr),"\\."))[seq(2,length(colnames(thiscorr))*2,2)]
    colnames(thiscorr)=listspp
    rownames(thiscorr)=listspp
    
    corrplot::corrplot(thiscorr,
                       method="ellipse",diag=F,order="alphabet",
                       tl.cex=0.5,cl.cex=0.5,number.cex=0.5,pch.cex=0.5,
                       title=paste("\nnsites chr",thischr),
                       type="upper",
                       mar=c(0,0,0,0))
  }
}
dev.off()

}

## combine by species -- relative dxy and fst
## combine by species 
if(13 %in% steps) {
bigplot = bigplot[order(bigplot$chr,bigplot$starts),]
bigplot$chrpos = paste(bigplot$chr,bigplot$starts)
bigplot$chrpos = as.numeric(as.factor(bigplot$chrpos))


qtab=table(bigplot[,c("species","sumquantile")])
ftab=table(bigplot[,c("species","quantileFST")])
dtab=table(bigplot[,c("species","quantiledxymeans")])

btemp=c()
for (i in 1:length(unique(bigplot$species))){
  spp=sort(unique(bigplot$species))[i]
  temp = bigplot[bigplot$species==spp,]
  temp=temp[,c("chr","starts","quantiledxymeans","quantileFST","sumquantile","chrpos")]
  colnames(temp)[3:5] = paste(colnames(temp)[3:5],spp,sep=".")
  if(is.null(btemp)){
    btemp=temp
  } else {
    btemp=merge(btemp,temp,by=c("chr","starts","chrpos"),all=T)
  }
}

dim(btemp)

dxys=c(seq(4,33,by=3))
fsts=c(seq(5,33,by=3))
nsit=c(seq(6,33,by=3))

dcor=cor(btemp[,dxys],use="pairwise.complete.obs")
fcor=cor(btemp[,fsts],use="pairwise.complete.obs")
ncor=cor(btemp[,nsit],use="pairwise.complete.obs")

listspp=unlist(strsplit(colnames(dcor),"\\."))[seq(2,length(colnames(dcor))*2,2)]
colnames(dcor)=listspp
rownames(dcor)=listspp
listspp=unlist(strsplit(colnames(fcor),"\\."))[seq(2,length(colnames(fcor))*2,2)]
colnames(fcor)=listspp
rownames(fcor)=listspp
listspp=unlist(strsplit(colnames(ncor),"\\."))[seq(2,length(colnames(ncor))*2,2)]
colnames(ncor)=listspp
rownames(ncor)=listspp

png("dxy_quan_corrplot_data.may2020.png")
plot(btemp[,dxys])
dev.off()

pdf("quans_corr_fullgenome.may2020.pdf",height=1.5,width=3)
par(mfrow=c(1,3))
corrplot::corrplot(dcor,
                   method="ellipse",diag=F,order="alphabet",
                   tl.cex=0.5,cl.cex=0.5,number.cex=0.5,pch.cex=0.5,
                   title="\ndxy",
                   type="upper",
                   mar=c(0,0,0,0))
corrplot::corrplot(fcor,
                   method="ellipse",diag=F,order="alphabet",
                   tl.cex=0.5,cl.cex=0.5,number.cex=0.5,pch.cex=0.5,
                   title="\nfsts",
                   type="upper",
                   mar=c(0,0,0,0))
corrplot::corrplot(ncor,
                   method="ellipse",diag=F,order="alphabet",
                   tl.cex=0.5,cl.cex=0.5,number.cex=0.5,pch.cex=0.5,
                   title="\nnsites",
                   type="upper",
                   mar=c(0,0,0,0))
dev.off()


chisq.test(qtab)
chisq.test(ftab)
chisq.test(dtab)
}
