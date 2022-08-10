## python script to create windows from a vcf 
import numpy as np
import sys, os, glob
import gzip
import pandas as pd
import io
import copy
from datetime import datetime

def subset_beagle(beaglefile,window_subset,outfilebasename):
    '''extract stuff from a beagle file based on chrom and position'''
    ## initialize some stuff
    color_dict = {}
    empty_lines = []
    print(beaglefile)
    ## read beagle
    if beaglefile[-3:] == ".gz":
        with gzip.open(beaglefile,"rb") as infile:
            lines = infile.readlines()
            lines = [i.decode('utf-8') for i in lines]
    else:
        with open(beaglefile,"r") as infile:
            lines = infile.readlines()
    ## separate metadata from data
    header = lines[0]
    data = lines[1:]
    locations=[line.strip().split("\t")[0] for line in data]
    chroms = [loc.split("_")[-2] for loc in locations]
    pos = [int(loc.split("_")[-1]) for loc in locations]
    linenum=list(range(1,len(lines)))
    ## turn beagle into a dataframe
    df_subset = pd.DataFrame(list(zip(chroms, pos, linenum)), 
        columns =['chroms', 'pos', "linenum"]) 
    unused_positions = list(df_subset["linenum"])
    
    ## iterate over the windows
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
    ## put extra lines aside
    if len(unused_positions) > 0:
        empty_lines += unused_positions
    ## iterate through the dictionary to get files to write out
    for color in color_dict.keys():
        lines_numbers = color_dict[color]
        lines_to_write = [lines[i] for i in lines_numbers]
        lines_to_write = [header] + lines_to_write
        outfile_name = str(outfilebasename)+"_"+str(color)+".beagle"
        with open(outfile_name,"a") as outfile:
            _ = outfile.writelines(lines_to_write)
    lines_to_write = [lines[i] for i in empty_lines]
    lines_to_write = [header] + lines_to_write
    outfile_name = str(outfilebasename)+"_NORM.beagle"
    with open(outfile_name,"a") as outfile:
        _ = outfile.writelines(lines_to_write)

try:
    beaglefolder = sys.argv[1]
    print("\tFolder is: ",beaglefolder,str(datetime.now()))
except:
    print("Folder not given, quitting",str(datetime.now()))
    exit()
    #beaglefolder="/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/BEAGLE/Amphispiza_bilineata/FULL/"
try:
    windowfile = str(sys.argv[2])
    print("\tWindow file is: ",windowfile,str(datetime.now()))
except:
    print("Window file not given, quitting",str(datetime.now()))
    exit()
    #windowfile="/home/kprovost/nas3/ANGSD_pipeline/bil.fst-100.outlier.SMALL.27July2022.txt"
    ## CURRENTLY these files contain only the outliers highs and lows
windows = pd.read_csv(windowfile,sep="\t")
windows["chr"]=windows["chr"].str.replace("PseudoNC_","")
windows["chr"]=windows["chr"].str.replace("\d+.1_Tgut_","")
unique_colors = list(windows["zscoresppFst"].unique())
unique_chroms = list(windows["chr"].unique())

for window_chr in unique_chroms:
    print("CHROMOSOME:",window_chr)
    this_folder = beaglefolder+str(window_chr)+"/"
    window_subset = windows.loc[windows["chr"]==str(window_chr)]
    print(this_folder)
    beaglefiles = glob.glob(this_folder+"*beagle*")
    for beaglefile in beaglefiles:
        #outfilebasename=beaglefolder+"subset_for_fst"
        if "subset_for_fst" not in beaglefile:
            outfilebasename = beaglefile+"subset_for_fst"
            subset_beagle(beaglefile=beaglefile,window_subset=window_subset,outfilebasename=outfilebasename)
