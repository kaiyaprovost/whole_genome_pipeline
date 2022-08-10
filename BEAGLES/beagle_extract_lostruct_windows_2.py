## python script to create windows from a vcf 

import numpy as np
import sys, os, glob
import gzip
import pandas as pd
import io
import copy
from datetime import datetime

try:
    beaglepath = str(sys.argv[1])
    print("\tBeagle path is: ",beaglepath,str(datetime.now()))
except:
    print("Beagle path not given, defaulting",str(datetime.now()))
    beaglepath="/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/BEAGLE/Vireo_bellii/WINDOWS/"
    print("\tBeagle path is: ",beaglepath,str(datetime.now()))

try:
    windowfile = str(sys.argv[2])
    print("\tWindow file is: ",windowfile,str(datetime.now()))
except:
    print("Window file not given, quitting",str(datetime.now()))
    exit()
    #windowfile="~/nas3/bel_mds_coords_COLORLOSTRUCT.csv"

windows = pd.read_csv(windowfile,usecols=["file","ccols"])

unique_colors = list(windows["ccols"].unique())

## get all the files in the directory that are beagles 
print("Finding beagle files in filepath",str(datetime.now()))
os.chdir(beaglepath)
beagle_files = glob.glob("**/*beagle.gz",recursive=True)
beagle_basenames = [os.path.basename(i) for i in beagle_files]

## loop through colors

for color in unique_colors:
    print(color,str(datetime.now()))
    outfile = beaglepath+os.path.basename(windowfile)+"."+color+".SUBSET.beagle.gz"
    thiscol = windows[windows["ccols"] == color]
    thisfiles = list(thiscol["file"].unique())
    shortfiles = [i.replace("Vireo_bellii_PCAngsd.1.cov","beagle.gz") for i in thisfiles]
    shortfiles = [i.replace("Amphispiza_bilineata_PCAngsd.1.cov","beagle.gz") for i in shortfiles]
    shortfiles = [i.replace("Campylorhynchus_brunneicapillus_PCAngsd.1.cov","beagle.gz") for i in shortfiles]
    shortfiles = [i.replace("Toxostoma_crissale_PCAngsd.1.cov","beagle.gz") for i in shortfiles]
    shortfiles = [i.replace("Toxostoma_curvirostre_PCAngsd.1.cov","beagle.gz") for i in shortfiles]
    shortfiles = [i.replace("Auriparus_flaviceps_PCAngsd.1.cov","beagle.gz") for i in shortfiles]
    shortfiles = [i.replace("Melozone_fusca_PCAngsd.1.cov","beagle.gz") for i in shortfiles]
    shortfiles = [i.replace("Polioptila_melanura_PCAngsd.1.cov","beagle.gz") for i in shortfiles]
    shortfiles = [i.replace("Phainopepla_nitens_PCAngsd.1.cov","beagle.gz") for i in shortfiles]
    shortfiles = [i.replace("Cardinalis_sinuatus_PCAngsd.1.cov","beagle.gz") for i in shortfiles]
    these_beagles_index = [beagle_basenames.index(i) for i in shortfiles]
    these_beagles = [beagle_files[i] for i in these_beagles_index]
    print("WRITING",str(datetime.now()))
    for beagle in these_beagles:
        _ = os.system("cat "+str(beaglepath)+str(beagle)+" >> "+outfile)
    
    