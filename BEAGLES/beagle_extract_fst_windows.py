## python script to create windows from a vcf 

import numpy as np
import sys, os, glob
import gzip
import pandas as pd
import io
import copy
from datetime import datetime

# cd /vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/BEAGLE/
# cd BRUNNEICAPILLUS
# ## might need to index gzipped files. 
# beaglefile=Campylorhynchus-brunneicapillus-called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.vcf.gz
# windowfile=/home/kprovost/nas3/LOSTRUCT/Campylorhynchus.brunneicapillus.called.geno.PseudoNC.all.fixedchroms.converted.sorted.sorted.sorted.nospace.coords_SMALL.csv
# python3 /home/kprovost/nas3/LOSTRUCT/vcf_extract_lostruct_windows.py $beaglefile $windowfile

try:
    beaglefile = sys.argv[1]
    print("\tFile is: ",beaglefile,str(datetime.now()))
except:
    print("Filename not given, quitting",str(datetime.now()))
    #exit()
    beaglefile="/Users/kprovost/Downloads/test.beagle"

try:
    windowfile = str(sys.argv[2])
    print("\tWindow file is: ",windowfile,str(datetime.now()))
except:
    print("Window file not given, quitting",str(datetime.now()))
    #exit()
    windowfile="/Users/kprovost/Downloads/test.lostruct.csv"

## beagle files only need the first column to work, as well as the header
## the first column will be <info>_chrom_position
## just need to extract the chroms and positions 


## to do this
## generate list of chrom_position from the windows file
## then do a bunch of greps on the beagle file?

## chrom,start,end,ccols
## PseudoNC_1,406,900405,black
## PseudoNC_1,900406,1000405,#1B9E77

## first grep the lines that contain the chromosome
## then go through. while statement perhaps? while value equal to or between values? 
## then 

windows = pd.read_csv(windowfile,sep="\t")
windows["chr"]=windows["chr"].str.replace("PseudoNC_","")
windows["chr"]=windows["chr"].str.replace("\d+.1_Tgut_","")

unique_colors = list(windows["zscoresppFst"].unique())
unique_chroms = list(windows["chr"].unique())

if beaglefile[-3:] == ".gz":
    with gzip.open(beaglefile,"rb") as infile:
        lines = infile.readlines()
    lines = [i.decode('utf-8') for i in lines]
else:
    with open(beaglefile,"r") as infile:
        lines = infile.readlines()

header = lines[0]
data = lines[1:]

locations=[line.strip().split("\t")[0] for line in data]
chroms = [loc.split("_")[-2] for loc in locations]
pos = [int(loc.split("_")[-1]) for loc in locations]
linenum=list(range(1,len(lines)))

df = pd.DataFrame(list(zip(chroms, pos, linenum)), 
               columns =['chroms', 'pos', "linenum"]) 

color_dict = {}
empty_lines = []

for window_chr in unique_chroms:
    print("CHROMOSOME:",window_chr)
    window_subset = windows.loc[windows["chr"]==window_chr]
    df_subset = df.loc[df["chroms"]==window_chr]
    unused_positions = list(df_subset["linenum"])
    for index, row in window_subset.iterrows():
        start = row["start"]
        try:
            end = row["stop"]
        except:
            end = row["end"]
        color = row["zscoresppFst"]
        this_df = df_subset.loc[(df_subset["pos"] >= start) & (df_subset["pos"] <= end)]
        if this_df.size > 0:
            line_value = color_dict.get(color,[])
            line_value += list(this_df["linenum"])
            color_dict[color] = line_value
            unused_positions = [x for x in unused_positions if x not in list(this_df["linenum"])]
    if len(unused_positions) > 0:
        empty_lines += unused_positions

for color in color_dict.keys():
    lines_numbers = color_dict[color]
    lines_to_write = [lines[i] for i in lines_numbers]
    lines_to_write = [header] + lines_to_write
    outfile_name = str(beaglefile)+"_"+str(color)+".vcf"
    with open(outfile_name,"a") as outfile:
        _ = outfile.writelines(lines_to_write)

lines_to_write = [lines[i] for i in empty_lines]
lines_to_write = [header] + lines_to_write
outfile_name = str(beaglefile)+"_empty.vcf"
with open(outfile_name,"a") as outfile:
    _ = outfile.writelines(lines_to_write)