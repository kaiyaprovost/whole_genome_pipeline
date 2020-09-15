## python script to create windows from a vcf 

import numpy as np
import sys, os, glob
import gzip
import pandas as pd
import io
import copy

## THIS SCRIPT IS NOT WORKING 


# cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/"
# python3 /Users/kprovost/Documents/Github/whole_genome_pipeline/Unsorted/vcf_extract_lostruct_windows.py ./called_geno/Amphispiza-bilineata-called.geno.fixedchroms.converted.sorted.nospace.vcf.gz ./LOSTRUCT/lostruct_results/finished/Amphispiza.bilineata.called.geno.fixedchroms.converted.sorted.nospace.coords.csv
# python3 /Users/kprovost/Documents/Github/whole_genome_pipeline/Unsorted/vcf_extract_lostruct_windows.py ./called_geno/amphispiza-test.vcf.gz ./LOSTRUCT/lostruct_results/finished/Amphispiza.bilineata.called.geno.fixedchroms.converted.sorted.nospace.coords.csv
# python3 /Users/kprovost/Documents/Github/whole_genome_pipeline/Unsorted/vcf_extract_lostruct_windows.py ./called_geno/Auriparus-flaviceps-called.geno.PseudoNC.all.sorted.sorted.sorted.nospace.vcf.gz ./LOSTRUCT/lostruct_results/finished/Auriparus.flaviceps.called.geno.PseudoNC.all.sorted.sorted.sorted.nospace.coords.csv
# python3 /Users/kprovost/Documents/Github/whole_genome_pipeline/Unsorted/vcf_extract_lostruct_windows.py ./called_geno/Vireo-bellii-called.geno.PseudoNC.all.sorted.sorted.sorted.nospace.vcf.gz ./LOSTRUCT/lostruct_results/finished/Vireo.bellii.called.geno.PseudoNC.all.sorted.sorted.sorted.nospace.coords.csv

try:
	vcffile = sys.argv[1]
	print("\tFile is: ",vcffile)
except:
	print("Filename not given, quitting")
	exit()

try:
	windowfile = str(sys.argv[2])
	print("\tWindow file is: ",windowfile)
except:
	print("Window file not given, quitting")
	exit()

def windows2vcf(vcflines,windows,header,vcffile):
	for index, row in windows.iterrows():
		if index % 100 == 0:
			print(str(index)+"/"+str(total_rows))
		subset = subset_by_window(vcflines,chrom=row["chrom"],start=row["start"],end=row["end"])
		this_color = row["ccols"]
		subset.to_csv(vcffile+"_"+this_color+".vcf",mode="a",index=False,sep="\t",header=False)

## this will be for a single window
def subset_by_window(vcflines,chrom,start,end):
	subset = vcflines.loc[(vcflines["#CHROM"] == chrom) & (vcflines["POS"] >= start) & (vcflines["POS"] <= end)]
	return(subset)

#subset_by_window(vcflines,chrom="Pseudo_NC_1",start,end)

windows = pd.read_csv(windowfile,usecols=["chrom","start","end","ccols"])

header = ""

unique_colors = list(windows["ccols"].unique())
total_rows=windows.shape[0]

next_index=0

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
					with open(vcffile+"_"+color+".vcf","w") as of:
						of.write("".join(header))
		else:
			try:
				vcfline = pd.read_csv(io.StringIO(''.join(vcfstring)), delim_whitespace=True,header=None)
			except:
				print("WARNING: Something wrong with VCF file input. Removing bad lines.")
				vcfline = pd.read_csv(io.StringIO(''.join(vcfstring)), delim_whitespace=True,header=None,error_bad_lines=False)
			#vcfline.columns = column_names
			found_window = False
			this_chrom = vcfline.iloc[0,0]
			this_position = vcfline.iloc[0,1]
			while (found_window == False) and (next_index < total_rows):
				chrom = windows.iloc[next_index,0] 
				start = windows.iloc[next_index,1] 
				end = windows.iloc[next_index,2] 
				#print(chrom,start,end)
				if this_position < start:
					#print("TOO EARLY")
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
					print(str(next_index)+"/"+str(total_rows)) 
			if next_index >= total_rows:
				next



#vcfstrings = [x.decode("utf-8") for x in lines]


# try:
# 	vcflines = pd.read_csv(io.StringIO(''.join(vcfstrings[2:])), delim_whitespace=True)
# except:
# 	print("WARNING: Something wrong with VCF file input. Removing bad lines.")
# 	vcflines = pd.read_csv(io.StringIO(''.join(vcfstrings[2:])), delim_whitespace=True,error_bad_lines=False)

## something is wrong with this vcffile

# windows2vcf(vcflines,windows,header,vcffile)