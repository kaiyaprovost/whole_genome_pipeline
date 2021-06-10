## python script to create windows from a vcf 

import numpy as np
import sys, os, glob
import gzip
import pandas as pd
import io
import copy
from datetime import datetime

try:
    beaglefile = sys.argv[1]
    print("\tFile is: ",beaglefile,str(datetime.now()))
except:
    print("Filename not given, quitting",str(datetime.now()))
    #exit()
    beaglefile="/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/BEAGLE/Auriparus-flaviceps-NOWEIRD.beagle.gz"

try:
    vcffolder = str(sys.argv[2])
    print("\VCF folder is: ",vcffolder,str(datetime.now()))
except:
    print("VCF folder not given, quitting",str(datetime.now()))
    #exit()
    vcffolder="/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/BEAGLE/"

## beagle files only need the first column to work, as well as the header
## the first column will be <info>_chrom_position
## just need to extract the chroms and positions 

print("Subsetting species",str(datetime.now()))
beagspp = os.path.basename(beaglefile).split(".")[0].replace("-NOWEIRD","").replace("called","")

## this is lengthy
print("Reading Beagle",str(datetime.now()))
beagle = pd.read_csv(beaglefile,sep="\t",compression="gzip",dtype="str")

print(beagle)

## need to fix the beagle columns that have . in them 
print("Fixing beagle names",str(datetime.now()))
beagle.columns = [j.split(".")[0] for j in beagle.columns]
beagle["chrom"] = [i.split("_")[-2] for i in beagle["marker"]]
beagle["pos"] = [i.split("_")[-1] for i in beagle["marker"]]
beagle["chrompos"] = beagle["chrom"]+"_"+beagle["pos"]

print("Getting VCF list",str(datetime.now()))
vcffilelist = glob.glob(vcffolder+"*"+beagspp+"*vcf.gz") # +glob.glob(vcffolder+"*"+beagspp+"*vcf")

for vcffile in vcffilelist:
    print(vcffile,str(datetime.now()))
    vcf = pd.read_csv(vcffile,sep="\t",header=2,usecols=[0,1],compression="gzip",dtype="str")
    subset = vcffile.split("_")[-1].split(".")[0]
    vcf=vcf.rename(columns={"#CHROM":"chrom","POS":"pos"})
    vcf["chrompos"] = vcf["chrom"]+"_"+vcf["pos"]
    print("Subsetting Beagle",str(datetime.now()))
    beagle_subset=beagle[beagle['chrompos'].isin(list(set(vcf["chrompos"])))]
    outfilename =  os.getcwd()+"/"+subset+"_"+os.path.basename(beaglefile)
    print(outfilename,str(datetime.now()))
    beagle_subset=beagle_subset.drop(["chrom","pos","chrompos"],axis=1)
    beagle_subset.to_csv(outfilename, sep='\t',compression="gzip")
    