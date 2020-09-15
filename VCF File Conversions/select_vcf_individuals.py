# python 

## usage: python3 select_ms_individuals.py $fullvcffile $numtosample

import random 
import numpy as np
import sys

## input the file

try:
	fullvcffile = sys.argv[1]
	print("\nRead the following filename:")
	print(fullvcffile)
except:
	print("\nNo filename given, quitting")
	sys.exit()
	fullvcffile = "/Users/kprovost/Documents/folder_for_popgenome/model1_panmixia_120k-1558968497-1-recap_1.11e-07-5e-07-4000.vcf"

try:
	numtosample = int(sys.argv[2])
	print("\nRead the following number to sample:")
	print(numtosample)
except:
	print("\nNo/improper sample number given, defaulting to 20")
	numtosample=20
	
try:
	locsfile = (sys.argv[3])
	print("\nRead the following locs file:")
	print(locsfile)
except:
	print("\nNo locs file given, quitting")
	sys.exit()
	locsfile="/Users/kprovost/Documents/folder_for_popgenome/modelnuc.locs"

## read in the file 
with open(fullvcffile,"r") as filename:
	lines = filename.readlines()

print("Numlines of fullvcffile: ",len(lines))

## get the data and the metadata
## four lines of metadata, one header
metadata = "".join(lines[0:5])
header = lines[5]
## actual data
data = lines[6:]

## input the number to sample total and then per pop
## this is num individuals, num genomes is double this
numtosampleperpop = numtosample//2

## get the number of genomes and individuals

splithead = header.split("\t")
## there are 9 slots (0-8) of pure header, rest are inds 
inds = splithead[9:]
numinds = len(inds)

## set the number of populations -- here always two
numpops = 2 

## extract individuals per populations
numperpop = numinds//2 

## extract individuals per population
firstpop = inds[:numperpop]
secondpop = inds[numperpop:]

## randomly select the sample of individuals
topickpop1 = np.random.choice(a=numperpop,size=numtosampleperpop,replace=False)
topickpop2 = numperpop+np.random.choice(a=numperpop,size=numtosampleperpop,replace=False)
## this is the number of individual, not the name

"SUBSETTING LOCS"
## get the locations and extract metadata
with open(locsfile,"r") as locname:
	locs = locname.readlines()

locheader=locs[1]
locdata=locs[2:]

## subset the individual locations
topickpopall = np.ndarray.tolist(topickpop1)+ np.ndarray.tolist(topickpop2)
topickpopall.sort()
subsetlocs = [locdata[j] for j in topickpopall]

"SUBSETTING GENOMES"
## process each line and only subset inds that are selected

tokeep = [0,1,2,3,4,5,6,7,8]
tokeep += [i+8 for i in topickpopall]


subsetdata = data
for i in range(len(subsetdata)):
	todo = subsetdata[i]
	splittodo = todo.strip().split("\t")
	## first 0-8 must be kept
	## the last data is at N+8
	## the first data is at 9
	keeptodo = [splittodo[j] for j in tokeep]
	subsetdata[i] = "\t".join(keeptodo)+"\n"

newheader = "\t".join([splithead[j] for j in tokeep])

filenameprefix = fullvcffile[:-4]
subsetfilename = filenameprefix+".subset.vcf"

print("WRITING GENOMES")
## write genomes then write locs
with open(subsetfilename,"w") as outfile:
	outfile.write(metadata)
	outfile.write(newheader)
	for i in subsetdata: 
		outfile.write(i)

print("WRITING LOCS")
subsetlocfilename = filenameprefix+".subsetlocs"
with open(subsetlocfilename,"w") as outfile2:
	outfile2.write(locheader)
	for i in subsetlocs: 
		outfile2.write(i)


