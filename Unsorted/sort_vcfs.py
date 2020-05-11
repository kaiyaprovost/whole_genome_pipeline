#!/bin/python3
import glob
import sys
import os

def vcfsort(vcffile):
	print("reading lines")
	with open(vcffile,"r") as infile:
		lines=infile.readlines()
	foundchrom=False
	linestart=0
	endheader=0
	print("finding header")
	while(foundchrom == False):
		#print(linestart)
		thisline=lines[linestart]
		if thisline[0:6] == "#CHROM":
			endheader=linestart
			print("Header lines: 0-"+str(endheader))
			foundchrom=True
		else:
			linestart+=1
	print("sorting lines and removing duplicates")
	headerlines=lines[0:endheader+1]
	linestosort = lines[endheader+1:]
	linestosort.sort()
	linestosort = [item for item in linestosort if item[0]!="#"]
	linestosort=list(set(linestosort))
	linestosort.sort()
	#y=np.asarray(x)
	for i in range(len(linestosort)):
		line = linestosort[i]
		try:
			chrom,pos,data = line.split("\t",2)
			linestosort[i]=(chrom,pos,data)
		except:
			linestosort[i]=(line)
	linestosort = sorted(linestosort, key=lambda i: (i[0], int(i[1])))
	sortedlines = ["\t".join(list(item)) for item in linestosort]
	print("generating final lines and writing")
	finallines = headerlines+sortedlines
	
	sortedname=os.path.basename(vcffile)
	sortedname=sortedname.replace(".vcf","")
	with open(sortedname+".sorted.vcf","w") as outfile:
		outfile.writelines(finallines)
	print("done")

try:
	vcffile = sys.argv[1]
	print("Read file/folder: "+vcffile)
except:
	print("No file/folder given, exiting")
	exit()
	#vcffile = "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/BRUNNEICAPILLUS/test.vcf"

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
		vcfsort(thisfile)
		print("Done "+str(file_index+1)+"/"+str(numfiles))
elif os.path.isfile(vcffile):  
	print("\nInput is a normal file")
	vcfsort(vcffile)
else:  
	print("\nCannot process this kind of input! Exiting")
	exit()




