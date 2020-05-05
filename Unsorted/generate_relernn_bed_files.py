#!/bin/python3
import glob
import sys
import os

def vcf2bed2mask(vcffile):
	with open(vcffile,"r") as infile:
		lines = infile.readlines()
	bed_dict = {}
	## iterate 
	for line in lines:
		if line[0] != "#":
			chrom,position = line.split("\t")[0:2]
			x = bed_dict.get(chrom,None)
			if x == None:
				bed_dict[chrom] = [int(position),0]
			else:
				bed_dict[chrom][1] = int(position)
	## once dict is built iterate through it and round positions to nearest 1000, round(x,-3)
	for key,value in bed_dict.items():
		print(key, value)
		## key is the chrom
		## value is [first position, last position]
		value[0] = max(-1000+round(int(value[0]),-3),1)
		value[1] = max(1000+round(int(value[1]),-3),1)
	## write mask and bed files
	with open(vcffile+".mask","w") as maskfile:
		for key,value in bed_dict.items():
			maskfile.write(str(key)+"\t0\t"+str(value[0])+"\n")
	with open(vcffile+".bed","w") as bedfile:
		for key,value in bed_dict.items():
			bedfile.write(str(key)+"\t0\t"+str(value[1])+"\n")

try:
	vcffile = sys.argv[1]
	print("Read file/folder: "+vcffile)
except:
	print("No file/folder given, exiting")
	exit()
	#vcffile = "/Users/kprovost/Documents/Github/whole_genome_pipeline/example.vcf"

## read in the file line by line and generate a dictionary
## one dictionary is starts and one dictionary is ends? 

if os.path.isdir(vcffile):  
	print("\nInput is a directory")
	folderpath=vcffile
	vcflist=glob.glob(pathname=folderpath+"/*.vcf")
	numfiles=len(vcflist)
	print("Read in "+str(numfiles)+" files")
	for file_index in range(len(vcflist)):
		thisfile = vcflist[file_index]
		#print(thisfile)
		vcf2bed2mask(thisfile)
		print("Done "+str(file_index+1)+"/"+str(numfiles))
elif os.path.isfile(vcffile):  
	print("\nInput is a normal file")
	vcf2bed2mask(vcffile)
else:  
	print("\nCannot process this kind of input! Exiting")
	exit()




