#!/usr/bin/env python3

import sys
import glob
import os 

try:
	infile = str(sys.argv[1])
	outfile = infile+".vcf"
	print("\tFile is: ",infile)
except:
	print("Filename not given, quitting")
	exit()

def geno2vcf1line(singleline,start=False):
	## read the line, separate the data types 
	first = singleline.strip()
	split = first.split("\t",4)
	chrom=split[0]
	position=split[1]
	major=split[2]
	minor=split[3]
	allinds=split[4]
	#chrom,position,major,minor,allinds=split
	major = major.upper()
	minor = minor.upper()
	missing = "N"
	
	## separate the individual genotypes and assign names
	inds = allinds.upper().split("\t")
	numinds = len(inds)
	namesinds = ["indiv" + str(s) for s in range(numinds)]
	#print(namesinds)
	
	
	## generate the first half of the data line
	vcflinefirst = chrom+"\t"+position+"\t.\t"+major+"\t"+minor+"\t.\tPASS\t.\tGT\t"
	
	## create the different kinds of genotypes possible 
	maj_maj = major+major
	maj_min = major+minor
	min_maj = minor+major
	min_min = minor+minor
	maj_mis = major+missing
	min_mis = minor+missing
	mis_maj = missing+major
	mis_min = missing+minor
	mis_mis = missing+missing
	
	## convert the actual genotypes to their corresponding vcf genotype, 0 = major, 1 = minor, . = missing
	genotypes = [maj_maj,maj_min,min_maj,min_min,maj_mis,min_mis,mis_maj,mis_min,mis_mis]
	vcftypes = ["0/0","0/1","1/0","1/1","0/.","1/.","./0","./1","./."]
	vcfinds = [ vcftypes[genotypes.index(i)] for i in inds ]
	
	## generate the final data line and output
	vcflinefinal = vcflinefirst+"\t".join(vcfinds)+"\n"
	
	if start == True:
		## generate the vcf header. 	## There are 8 fixed fields per record. All data lines are tab-delimited. In all cases, missing values are specified with a dot (”.”)
		vcfheader = "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"+"\t".join(namesinds)+"\n"
		return([vcflinefinal,vcfheader])
	else:
		return(vcflinefinal)
	
with open(infile,"r") as input:
	with open(outfile,"w") as output:
		count=1
		line = input.readline()
		convertedline,vcfheader = geno2vcf1line(line,start=True)
		output.write(vcfheader)
		output.write(convertedline)
		while line:
			line = input.readline()
			count+=1
			if len(line) != 0:
				convertedline = geno2vcf1line(line,start=False)
				output.write(convertedline)
				if (count % 1000000) == 0:
					print(count)
