## python script to create windows from a vcf 

import numpy as np
import sys, os, glob
import gzip
import pandas as pd
import io
import copy
from datetime import datetime

## THIS SCRIPT IS NOT WORKING BECAUSE OF NAME MISMATCH
## AS OF 14 JAN 2021
## the chromosome names don't match -- PseudoNC is still there, need to subset them?
## line by line not working -- the other one should be 

# cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/"
# pythonfile="/Users/kprovost/Documents/Github/whole_genome_pipeline/VCF File Conversions/vcf_extract_fstoutlier_windows.py"
# vcffile="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/lostruct_2021_done/Melozone-fusca_fixed.sorted.nospace.sorted.nospace.vcf.gz"
# coordfile="/Users/kprovost/Dropbox (AMNH)/Dissertation/fus.fstoutlier.small.csv"
# python3 "$pythonfile" "$vcffile" "$coordfile"

linebyline=False

try:
	vcffile = sys.argv[1]
	print("\tFile is: ",vcffile,str(datetime.now()))
except:
	print("Filename not given, quitting",str(datetime.now()))
	exit()
	vcffile="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/lostruct_2021_done/Melozone-fusca_fixed.sorted.nospace.sorted.nospace.vcf.gz"

try:
	windowfile = str(sys.argv[2])
	print("\tWindow file is: ",windowfile,str(datetime.now()))
except:
	print("Window file not given, quitting",str(datetime.now()))
	exit()
	windowfile="/Users/kprovost/Dropbox (AMNH)/Dissertation/fus.fstoutlier.small.csv"

def windows2vcf(vcflines,windows,header,vcffile,verbose=False):
	total_rows=windows.shape[0]
	print(total_rows,str(datetime.now()))
	for index, row in windows.iterrows():
		if index % 100 == 0 and verbose==True:
			print(str(index)+"/"+str(total_rows),str(datetime.now()))
			print(vcflines,str(datetime.now()))
		subset = subset_by_window(vcflines,row["chr"],row["start"],row["stop"])
		if index % 100 == 0 and verbose==True:
			print(subset,str(datetime.now()))
			print("CHROM START END",row["chr"],row["start"],row["stop"],str(datetime.now()))
		this_outlier = row["outlier"]
		subset.to_csv(str(vcffile)+"_"+str(this_outlier)+".vcf",mode="a",index=False,sep="\t",header=False)

## this will be for a single window
def subset_by_window(vcflines,chr,start,stop):
	subset = vcflines.loc[(vcflines["#CHROM"] == chr) & (vcflines["POS"] >= start) & (vcflines["POS"] <= stop)]
	return(subset)

#subset_by_window(vcflines,chrom="Pseudo_NC_1",start,end)

windows = pd.read_csv(windowfile,usecols=["chr","start","stop","outlier"])
unique_outliers = list(windows["outlier"].unique())


if linebyline==True:
	print("RUNNING LINE BY LINE")
	next_index=0
	header = ""
	with gzip.open(vcffile,"rb") as vcf:
		for line in vcf: 
			vcfstring = line.decode("utf-8")
			#print(vcfstring)
			if vcfstring[0] == "#": ## its a comment
				header += vcfstring
				if vcfstring[0:6] == "#CHROM":
					#column_names=line.strip().split("\t")
					for color in unique_colors:
						print(color)
						with open(vcffile+"_"+str(color)+".vcf","a") as of:
							of.write("".join(header))
					with open(vcffile+"_skip.vcf","a") as of2:
						of2.write("".join(header))
					with open(vcffile+"_early.vcf","a") as of3:
						of3.write("".join(header))
			else:
				try:
					vcfline = pd.read_csv(io.StringIO(''.join(vcfstring)), delim_whitespace=True,header=None)
				except:
					print("WARNING: Something wrong with VCF file input. Removing bad lines.")
					vcfline = pd.read_csv(io.StringIO(''.join(vcfstring)), delim_whitespace=True,header=None,error_bad_lines=False)
				#vcfline.columns = column_names
				found_window = False
				skip_mtDNA = False
				this_chrom = vcfline.iloc[0,0]
				if this_chrom[0:8] == "PseudoNC":
					chrom_split = this_chrom.split("_")
					this_chrom = chrom_split[-1]
				this_position = vcfline.iloc[0,1]
				#print(this_chrom,this_position)
				
				#while (found_window == False) and (next_index < total_rows): ## total_rows is not defined here
				while (found_window == False):
					chrom = windows.iloc[next_index,0] 
					start = windows.iloc[next_index,1] 
					end = windows.iloc[next_index,2] 
					#print(chrom,start,end)
					if(chrom=="1" and this_chrom == "mtDNA"):
						vcfline.to_csv(vcffile+"_skip.vcf",mode="a",index=False,sep="\t",header=False)
						found_window = True
					else:
						if this_position < start:
							#print("TOO EARLY")
							vcfline.to_csv(vcffile+"_early.vcf",mode="a",index=False,sep="\t",header=False)
							break
						elif (chrom == this_chrom) and (this_position >= start) and (this_position <= end):
							## get the color of the found window 
							this_color = windows.iloc[next_index,3]
							vcfline.to_csv(vcffile+"_"+this_color+".vcf",mode="a",index=False,sep="\t",header=False)
							## set found window to true 
							#print(windows.iloc[next_index])
							#print("FOUND")
							found_window = True 
						else:
							next_index += 1
							#print("NOT FOUND")
							#if next_index % 100 == 0:
							#print(str(next_index)+"/"+str(total_rows)) 
				
				## total_rows is not defined here
				#if next_index >= total_rows:
				#	next
else:
	print("STARTING TO READ VCF:",str(datetime.now()))
	with gzip.open(vcffile,"rb") as vcf:
		lines=vcf.readlines()
	
	print("DONE READING VCF:",str(datetime.now()))
	vcfstrings = [x.decode("utf-8") for x in lines]
	print("DONE CONVERTING VCF:",str(datetime.now()))
	header=vcfstrings[0:3]
	for outlier in unique_outliers:
		print(outlier)
		with open(vcffile+"_"+str(outlier)+".vcf","w") as of:
			_=of.write("".join(header))

	with open(vcffile+"_skipoutlier.vcf","w") as of2:
		_=of2.write("".join(header))
	
	with open(vcffile+"_earlyoutlier.vcf","w") as of3:
		_=of3.write("".join(header))
	
	print("read vcflines in as pandas",str(datetime.now()))
	try:
		vcflines = pd.read_csv(io.StringIO(''.join(vcfstrings[2:])), delim_whitespace=True)
	except:
		print("WARNING: Something wrong with VCF file input. Removing bad lines.")
		vcflines = pd.read_csv(io.StringIO(''.join(vcfstrings[2:])), delim_whitespace=True,error_bad_lines=False)
	print("done reading vcflines in as pandas",str(datetime.now()))

	## make sure the vcflines match -- remove any PseudoNC stuff, and remove any lines that are PseudoNW or NW or SS or NC. 
	vcflines["#CHROM"]=vcflines["#CHROM"].str.replace("PseudoNC_","")
	vcflines["#CHROM"]=vcflines["#CHROM"].str.replace("\d+.1_Tgut_","")
	windows["chr"]=windows["chr"].str.replace("PseudoNC_","")
	windows["chr"]=windows["chr"].str.replace("\d+.1_Tgut_","")
	print("done string replace of chrom fields",str(datetime.now()))

	for outlier in unique_outliers:
		print(outlier,str(datetime.now()))
		this_outlier_window = windows.loc[(windows["outlier"] == outlier)]
		unique_chroms = list(this_outlier_window["chr"].unique())
		for this_chrom in unique_chroms:
			print(this_chrom,str(datetime.now()))
			this_chrom_window = this_outlier_window.loc[(this_outlier_window["chr"] == this_chrom)]
			## is this just not printing? 
			windows2vcf(vcflines,this_chrom_window,header,vcffile)


## something is wrong with this vcffile

# windows2vcf(vcflines,windows,header,vcffile)