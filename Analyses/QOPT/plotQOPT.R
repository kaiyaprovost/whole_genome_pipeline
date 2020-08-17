library(plotrix)
library(rgdal)

reloadEnv = F

## import color scheme
library(RColorBrewer)
col = brewer.pal(12,"Paired")
red=col[6]
blue=col[2]
green=col[3]
yellow=col[7]
grey=col[9]

if (reloadEnv == T) {
  Env = raster::stack('/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/bio_2-5m_bil/bio1.bil'
    #list.files(
      #path = '/Users/kprovost/Dropbox (AMNH)/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/',
      #pattern = "\\.bil$",
      #full.names = T
    #)[1]
  )
  #ext = raster::extent(c(-125,-60, 10, 50)) ## make sure this will play nice with your points
  #ext = raster::extent(c(-115,-97,26,37))
  ext = raster::extent(c(-117, -98, 27, 35))
  Env = raster::crop(Env, ext)
  bg = Env[[1]] ## just for plotting
  
  border = rgdal::readOGR(
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/mexstates/mexico_and_states_subset.shp"
  )
  #plot(border)
  
}

## with all combined

structure_all = read.csv(
  "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/AllSpeciesMetadata_allK_9june2020.csv",
  row.names = NULL
)

#structure_all = structure_all[!(is.na(structure_all$P2_A)), ]

## chromosome averages
pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_AllC.pdf",
    width=9.5,height=4)
par(mfrow=c(2,2),mar=c(0,0,1,0))
for (s in 1:length(unique(structure_all$SP))) {
  
  spp = as.character(unique(structure_all$SP)[s])
  print(spp)
  structure_big = structure_all[structure_all$SP == spp, ]
  
  for (k in 2:5) {
    
    if(k == 2) {
      ourcols = c(42:43)
    } else if (k==3) {
      ourcols = c(44:46)
    } else if (k==4) {
      ourcols = c(45:48)
    } else if (k==5) {
      ourcols = c(49:53)
    } 
    
    
    
    raster::plot(bg,
                 col = "grey",
                 colNA = "white",
                 legend = F,
                 xaxt="n",yaxt="n",
                 main = spp#,xlim=c(-117, -98),ylim=c(27, 35)
                 
    )
    plot(border, add = T)
    points(structure_big$JLONG, structure_big$JLAT, pch = 4)
    for (i in 1:nrow(structure_big)) {
      
      if(!(is.na(structure_big[i,ourcols[1]]))) {
        
        print(i)
        floating.pie(
          xpos = structure_big$JLONG[i],
          ypos = structure_big$JLAT[i],
          x = c(as.numeric(structure_big[i,c(ourcols)])),
          radius = 0.3,
          col = c(blue,green,red,yellow,grey)
        )
      }
    }
    #dev.off()
    
    
  }
}
dev.off()


pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_P2.pdf")
for (s in 1:length(unique(structure_all$SP))) {
  spp = as.character(unique(structure_all$SP)[s])
  print(spp)
  structure_big = structure_all[structure_all$SP == spp, ]
  
  #png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/",
  #          spp,"_Metadata_P2_AP2.png",
  #          sep=""),width=760,
  #    height=460)
  
  
  
  raster::plot(bg,
               col = "grey",
               colNA = "white",
               legend = F,
               main = spp#,xlim=c(-117, -98),ylim=c(27, 35)
               
  )
  plot(border, add = T)
  points(structure_big$JLONG, structure_big$JLAT, pch = 4)
  for (i in 1:nrow(structure_big)) {
    
    if(!(is.na(structure_big$P2_A[i]))) {
      
      print(i)
      floating.pie(
        xpos = structure_big$JLONG[i],
        ypos = structure_big$JLAT[i],
        x = c(structure_big$P2_A[i], structure_big$P2_B[i]),
        radius = 0.3,
        col = c(blue,grey)
      )
    }
  }
  #dev.off()
  
  
  
}
dev.off()

pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_P3.pdf")
for (s in 1:length(unique(structure_all$SP))) {
  spp = as.character(unique(structure_all$SP)[s])
  print(spp)
  structure_big = structure_all[structure_all$SP == spp, ]
  
  #png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/",
  #          spp,"_Metadata_P2_AP2.png",
  #          sep=""),width=760,
  #    height=460)
  
  
  
  raster::plot(bg,
               col = "grey",
               colNA = "white",
               legend = F,
               main = spp#,xlim=c(-117, -98),ylim=c(27, 35)
               
  )
  plot(border, add = T)
  points(structure_big$JLONG, structure_big$JLAT, pch = 4)
  for (i in 1:nrow(structure_big)) {
    
    if(!(is.na(structure_big$P3_A[i]))) {
      
      print(i)
      floating.pie(
        xpos = structure_big$JLONG[i],
        ypos = structure_big$JLAT[i],
        x = c(structure_big$P3_A[i], structure_big$P3_B[i], structure_big$P3_C[i]),
        radius = 0.3,
        col = c(blue,grey,red)
      )
    }
  }
  #dev.off()
  
  
  
}
dev.off()

pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_P4.pdf")
for (s in 1:length(unique(structure_all$SP))) {
  spp = as.character(unique(structure_all$SP)[s])
  print(spp)
  structure_big = structure_all[structure_all$SP == spp, ]
  
  #png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/",
  #          spp,"_Metadata_P2_AP2.png",
  #          sep=""),width=760,
  #    height=460)
  
  
  
  raster::plot(bg,
               col = "grey",
               colNA = "white",
               legend = F,
               main = spp#,xlim=c(-117, -98),ylim=c(27, 35)
               
  )
  plot(border, add = T)
  points(structure_big$JLONG, structure_big$JLAT, pch = 4)
  for (i in 1:nrow(structure_big)) {
    
    if(!(is.na(structure_big$P4_A[i]))) {
      
      print(i)
      floating.pie(
        xpos = structure_big$JLONG[i],
        ypos = structure_big$JLAT[i],
        x = c(structure_big$P4_A[i], structure_big$P4_B[i],structure_big$P4_C[i], structure_big$P4_D[i]),
        radius = 0.3,
        col = c(blue,grey,red,yellow)
      )
    }
  }
  #dev.off()
  
  
  
}
dev.off()

pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_P5.pdf")
for (s in 1:length(unique(structure_all$SP))) {
  spp = as.character(unique(structure_all$SP)[s])
  print(spp)
  structure_big = structure_all[structure_all$SP == spp, ]
  
  #png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/",
  #          spp,"_Metadata_P2_AP2.png",
  #          sep=""),width=760,
  #    height=460)
  
  
  
  raster::plot(bg,
               col = "grey",
               colNA = "white",
               legend = F,
               main = spp#,xlim=c(-117, -98),ylim=c(27, 35)
               
  )
  plot(border, add = T)
  points(structure_big$JLONG, structure_big$JLAT, pch = 4)
  for (i in 1:nrow(structure_big)) {
    
    if(!(is.na(structure_big$P5_A[i]))) {
      
      print(i)
      floating.pie(
        xpos = structure_big$JLONG[i],
        ypos = structure_big$JLAT[i],
        x = c(structure_big$P5_A[i], structure_big$P5_B[i], structure_big$P5_C[i], structure_big$P5_D[i], structure_big$P5_E[i]),
        radius = 0.3,
        col = c(blue, green,red,yellow,grey)
      )
    }
  }
  #dev.off()
  
  
  
}
dev.off()

## COMBINED PCA

pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_AllP.pdf",
    width=9.5,height=4)
par(mfrow=c(2,2),mar=c(0,0,1,0))
for (s in 1:length(unique(structure_all$SP))) {
  
  spp = as.character(unique(structure_all$SP)[s])
  print(spp)
  structure_big = structure_all[structure_all$SP == spp, ]
  
  for (k in 2:5) {
  
    if(k == 2) {
      ourcols = c(27:28)
    } else if (k==3) {
      ourcols = c(30:32)
    } else if (k==4) {
      ourcols = c(33:36)
    } else if (k==5) {
      ourcols = c(37:41)
    } 
    
  

  raster::plot(bg,
               col = "grey",
               colNA = "white",
               legend = F,
               xaxt="n",yaxt="n",
               main = spp#,xlim=c(-117, -98),ylim=c(27, 35)
               
  )
  plot(border, add = T)
  points(structure_big$JLONG, structure_big$JLAT, pch = 4)
  for (i in 1:nrow(structure_big)) {
    
    if(!(is.na(structure_big[i,ourcols[1]]))) {
      
      print(i)
      floating.pie(
        xpos = structure_big$JLONG[i],
        ypos = structure_big$JLAT[i],
        x = c(as.numeric(structure_big[i,c(ourcols)])),
        radius = 0.3,
        col = c(blue,green,red,yellow,grey)
      )
    }
  }
  #dev.off()
  
  
  }
}
dev.off()

##

pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_K2.pdf")
for (s in 1:length(unique(structure_all$SP))) {
  spp = as.character(unique(structure_all$SP)[s])
  print(spp)
  structure_big = structure_all[structure_all$SP == spp, ]
  
  #png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/",
  #          spp,"_Metadata_K2.png",
  #          sep=""),width=760,
  #    height=460)
  
  
  
  raster::plot(bg,
               col = "grey",
               colNA = "white",
               legend = F,
               main = spp#,xlim=c(-117, -98),ylim=c(27, 35)
               
  )
  plot(border, add = T)
  points(structure_big$JLONG, structure_big$JLAT, pch = 4)
  for (i in 1:nrow(structure_big)) {
    
    if(!(is.na(structure_big$K2_A[i]))) {
      
      print(i)
      floating.pie(
        xpos = structure_big$JLONG[i],
        ypos = structure_big$JLAT[i],
        x = c(structure_big$K2_A[i], structure_big$K2_B[i]),
        radius = 0.3,
        col = c(blue,green)
      )
    }
  }
  #dev.off()
  
  
  
}
dev.off()

pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_K3.pdf")
for (s in 1:length(unique(structure_all$SP))) {
  spp = as.character(unique(structure_all$SP)[s])
  print(spp)
  structure_big = structure_all[structure_all$SP == spp, ]
  
  #png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/",
  #          spp,"_Metadata_K3.png",
  #          sep=""),width=760,
  #    height=460)
  
  
  
  raster::plot(bg,
               col = "grey",
               colNA = "white",
               legend = F,
               main = spp#,xlim=c(-117, -98),ylim=c(27, 35)
               
  )
  plot(border, add = T)
  points(structure_big$JLONG, structure_big$JLAT, pch = 4)
  for (i in 1:nrow(structure_big)) {
    
    if(!(is.na(structure_big$K3_A[i]))) {
      
      print(i)
      floating.pie(
        xpos = structure_big$JLONG[i],
        ypos = structure_big$JLAT[i],
        x = c(structure_big$K3_A[i], structure_big$K3_B[i], structure_big$K3_C[i]),
        radius = 0.3,
        col = c(blue,green, red)
      )
    }
  }
  #dev.off()
  
  
  
}
dev.off()

pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_K4.pdf")
for (s in 1:length(unique(structure_all$SP))) {
  spp = as.character(unique(structure_all$SP)[s])
  print(spp)
  structure_big = structure_all[structure_all$SP == spp, ]
  
  #png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/",
  #          spp,"_Metadata_K4.png",
  #          sep=""),width=760,
  #    height=460)
  
  
  
  raster::plot(bg,
               col = "grey",
               colNA = "white",
               legend = F,
               main = spp#,xlim=c(-117, -98),ylim=c(27, 35)
               
  )
  plot(border, add = T)
  points(structure_big$JLONG, structure_big$JLAT, pch = 4)
  for (i in 1:nrow(structure_big)) {
    
    if(!(is.na(structure_big$K4_A[i]))) {
      
      print(i)
      floating.pie(
        xpos = structure_big$JLONG[i],
        ypos = structure_big$JLAT[i],
        x = c(structure_big$K4_A[i], structure_big$K4_B[i], structure_big$K4_C[i], structure_big$K4_D[i]),
        radius = 0.3,
        col = c(blue,green, red,yellow)
      )
    }
  }
  #dev.off()
  
  
  
}
dev.off()

pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_K5.pdf")
for (s in 1:length(unique(structure_all$SP))) {
  spp = as.character(unique(structure_all$SP)[s])
  print(spp)
  structure_big = structure_all[structure_all$SP == spp, ]
  
  #png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/",
  #          spp,"_Metadata_K5.png",
  #          sep=""),width=760,
  #    height=460)
  
  
  
  raster::plot(bg,
               col = "grey",
               colNA = "white",
               legend = F,
               main = spp#,xlim=c(-117, -98),ylim=c(27, 35)
               
  )
  plot(border, add = T)
  points(structure_big$JLONG, structure_big$JLAT, pch = 4)
  for (i in 1:nrow(structure_big)) {
    
    if(!(is.na(structure_big$K5_A[i]))) {
      
      print(i)
      floating.pie(
        xpos = structure_big$JLONG[i],
        ypos = structure_big$JLAT[i],
        x = c(structure_big$K5_A[i], structure_big$K5_B[i], structure_big$K5_C[i], structure_big$K5_D[i], structure_big$K5_E[i]),
        radius = 0.3,
        col = c(blue,green, red,yellow,grey)
      )
    }
  }
  #dev.off()
  
  
  
}
dev.off()

# combine both NGS and PCA

#pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/AllSpeciesMetadata_BothPCandK.pdf")
png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/",
          "All_Metadata_PCandK.png",
          sep=""),width=2280,
    height=1520)
par(mfrow=c(4,3))
for (s in 1:length(unique(structure_all$SP))) {
  spp = as.character(unique(structure_all$SP)[s])
  print(spp)
  structure_big = structure_all[structure_all$SP == spp, ]

  #png(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT/",
  #         spp,"_Metadata_PCandK.png",
  #         sep=""),width=760,
  #   height=460)

  raster::plot(bg,
               col = "grey",
               colNA = "white",
               legend = F,
               main = spp#,xlim=c(-117, -98),ylim=c(27, 35)
               
  )
  plot(border, add = T)
  points(structure_big$JLONG, structure_big$JLAT, pch = 4)
  for (i in 1:nrow(structure_big)) {
    
    if(!(is.na(structure_big$K2_A[i]))) {
      if(!(is.na(structure_big$P2_A[i]))) {  
      print(i)
      floating.pie(
        xpos = structure_big$JLONG[i],
        ypos = structure_big$JLAT[i],
        x = c(structure_big$K2_A[i], structure_big$K2_B[i],structure_big$P2_A[i], structure_big$P2_B[i]),
        radius = 0.3,startpos = pi/2,
        col = c(blue,green,red,yellow)
      )
      }
    }
  }
  #dev.off()
  
  
  
}
dev.off()




#####


#states = rgdal::readOGR("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/cb_2016_us_state_500k/cb_2016_us_state_500k_SUBSET.shp")
#plot(states)
#mex = rgdal::readOGR("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/mexstates/mexstates_SUBSET.shp")
#plot(mex)

individuals  = list.files(path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/bamlists",
                          pattern = "\\.indlist$",
                          full.names = T)
qopts = list.files(path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/QOPT",
                   pattern = "\\.qopt$",
                   full.names = T)
species_to_output = c(
  "bilineata",
  "flaviceps",
  "brunneicapillus",
  "sinuatus",
  "fusca",
  "nitens",
  "melanura",
  "crissale",
  "curvirostre",
  "bellii"
)
