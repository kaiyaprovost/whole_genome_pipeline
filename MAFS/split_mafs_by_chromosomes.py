#!/usr/bin/python3 -u

import sys
import os
import glob
import copy

# cd "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MAFS/"; 
# for i in Phainopepla-nitens*mafs; do 
# python3 ~/nas3/ANGSD_pipeline/split_dxy_by_chromosomes.py $i 1; 
# done;

def checkfile(path,stufftowrite,header):
	if os.path.exists(path) == False:
		with open(path,"w") as outfile:
			outfile.write(header)
	with open(path,"a") as outfile:
		outfile.writelines(stufftowrite)

try:
	toconvert = sys.argv[1]
	print("\nRead the following file to convert:")
	print(toconvert)
except:
	print("\nNo file to convert given, quitting")
	sys.exit()
	#toconvert = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/raw/cri_SON_Dxy_persite_chrfix.txt"
try:
	chromcol = int(sys.argv[2])
	print("\nRead the following chromosome column:")
	print(chromcol)
except:
	print("\nNo chromosome column given, defaulting to first column")
	chromcol = 1

print("STARTING")
with open(toconvert,"r") as infile:
	lines = infile.readlines()

nlines = len(lines)

blocksize = 100000
originalblock = copy.deepcopy(blocksize)
print("Original Blocksize:"+str(originalblock))
starts = list(range(0,nlines,blocksize))
starts[0] = 1
stops = [k -1 for k in starts]
stops=stops[1:]+[nlines]

currstart=1
currstop=blocksize-1
nextstart=min(currstop+1,nlines)
nextstop=min(nextstart+blocksize-1,nlines)

header=lines[0]
fivepercent = round(nlines/10)

#print("|    |    |    |    |    |")
#print("|",end="")

while currstop < nlines:
#for j in  range(len(starts)):
	#start = starts[j]
	#stop = stops[j]
	start = currstart
	stop = currstop
	print(str(start)+"-"+str(stop))
	startline = lines[start]
	stopline = lines[stop]
	startchr = startline.split("\t")[chromcol-1].strip()
	stopchr = stopline.split("\t")[chromcol-1].strip()
	if(startchr == stopchr):
		blocksize = copy.deepcopy(originalblock)
		path = toconvert[:-4]+"-chr"+str(startchr)+toconvert[-4:]
		checkfile(path,lines[start:stop+1],header)
		nextstart=min(currstop+1,nlines)
		nextstop=min(nextstart+blocksize-1,nlines)
		currstart=nextstart
		currstop=nextstop
	else:
		print("MISMATCH")
		print(startchr+"-"+stopchr)
		blocksize=int(blocksize*0.5)
		print("new blocksize: "+str(blocksize))
		#for i in range(start,stop,1):
		#	#if (i % fivepercent == 0):
		#	#	print("-",end="")
		#	#if (i % 100000) == 0:
		#	#	print(str(i)+"/"+str(nlines)+"\t"+str(round(float(i/nlines*100),2))+"%",end="\n")
		#	#elif (i % 1000) == 0:
		#	#	print(str(i),end=" ")
		#	line = lines[i]
		#	split=line.split("\t")
		#	chrom = split[chromcol-1]
		#	path = toconvert[:-4]+"-chr"+str(chrom)+toconvert[-4:]
		#	checkfile(path,line,header)
		currstart=currstart
		currstop=currstart+blocksize-1



