#!/usr/bin/env python3

## Written by Kaiya L. Provost, 18 September 2019
## version 1.1
## version 1.2 updated 23 April 2020 -- written with python 3.7.0

import sys
import os
import glob
import numpy as np

## usage
## python3 "/Users/kprovost/Dropbox (AMNH)/BrianParrots/positionToN.py" "/Users/kprovost/Dropbox (AMNH)/BrianParrots/mafft-nexus-clean-75p-unique-taxa.phylip.fasta" 6680 7407 78336 109962 111773 131597 132346 138455 138677 163418 175140 186913 186914 197589 197607 212745 212801 217236 220790 267545 271043 318033 319907 319924 328857 336690 355175 374574 387361 409228 435474 445516 492174 521944 525434 525876 546470 563972 590289 616279 631050 645515 653024 660129 662049 667254 708672 730923 778037 789776 856840 856849 937503 942413 955547 979224 1004491 1007727 1016923 1018567 1032195 1053733 1107345 1128088 1186450 1210570 1226231 1244506 1283454 1315189 1322437 1338838 1348272 1360033 1378039 1380174 1405727 1409473 1417442 1417479 1431862 1451368 1473085 1545473 1551272 1583604 1596039 1626125 1637256 1649084 1650320 1652088 1660426 1664885 1666181 1666182 1675466 1682843 1685674 1685678 1685705 1685720 1685891 1685935 1686040 1686101 1714091 1715862 1742380 1800191 1875777 1887749 1888337 1888338 1890170 1921845 1924611 1934298 


def readFasta(fastafile,verbose=False):
	with open(fastafile) as fasta:
		text = fasta.read()
	## split the fasta file by ">", remove first because blank
	if verbose==True:
		print("\tFinished Reading Fasta")
	split = text.split(">")[1:]
	if verbose==True:
		print("\tSplit Text (1-5):")
		print(split[0:5])	
	#fastadict = {}
	#seqarray = ""
	indlist = []
	first=True
	
	## split each value of the list by the first "\n"
	for line in split:
		linesplit = line.split("\n",1)
		individual = linesplit[0]
		if verbose==True:
			print("\t Individual:")
			print(individual)
		indlist.append(individual)
		
		seq = np.array(list(linesplit[1].replace("\n","")))
		
		if (first == True):
			seqarray = seq
			first = False
		else:
			seqarray = np.vstack((seqarray,seq))
		
	return([indlist,seqarray])

def baseToNull(seqarray,positionlist,basefile,indlist,verbose=False):
	columnNum = [int(i)-1 for i in positionlist]
	postext = [str(i) for i in positionlist]
	bases = seqarray[:,columnNum]
	numinds = len(indlist)
	
	if verbose==True:
		print("\tColumn numbers (1-5):")
		print(columnNum[0:5])
		print("\tPosition Text (1-5):")
		print(postext[0:5])
		print("\tBases (1-5):")
		print(bases[0:5])
		print("\tNuminds: "+str(numinds))
		print("\tWriting")
	with open(basefile,"a") as out:
		out.write("##POSITIONS: "+" ".join(postext)+"\n") 
		for indnum in range(numinds):
			individual = indlist[indnum]
			seq = bases[indnum,:]
			seqlist = seq.tolist()
			joined = "".join(seqlist)
			towrite = ">"+individual+"\n"+joined+"\n"
			out.write(towrite)
			
	if verbose==True:
		print("\tDone Writing")
			
	bases[bases != "N"] = "N"
	seqarray[:,columnNum] = "N"	
	if verbose==True:
		print("\tBases again (1-5):")
		print(bases[0:5])
	
		

def seqarrayToFasta(indlist,seqarray,outfile,verbose=False):
	numinds = len(indlist)
	with open(outfile,"a") as out:
		for indnum in range(numinds):
			individual = indlist[indnum]
			seq = seqarray[indnum,:]
			seqlist = seq.tolist()
			joined = "".join(seqlist)
			towrite = ">"+individual+"\n"+joined+"\n"
			out.write(towrite)


def main():
	try:
		fastafile = sys.argv[1]
		print("\nRead the following fasta:")
		print(fastafile)
	except:
		print("\nNo fasta given, quitting")
		sys.exit()
		#fastafile = "/Users/kprovost/Dropbox (AMNH)/BrianParrots/mafft-nexus-clean-75p-unique-taxa.phylip.fasta"
	
	try:
		positionlist = sys.argv[2:]
		print("\nRead the following positions:")
		print(positionlist)
	except:
		print("\nNo positions given, quitting")
		sys.exit()
	
	verbose=False 
	
	outfile = fastafile+".converted"
	print("\nOutfile will be: "+outfile)

	basefile = fastafile+".bases"
	print("\nBases file will be: "+basefile)
	
	print("\nReading in fasta file")
	indlist,seqarray = readFasta(fastafile,verbose)
	
	if verbose==True:
		print("\tIndlist (1-5):")
		print(indlist[0:5])
		print("\tSeqarray (1-5):")
		print(seqarray[0:5])
	
	print("\nConverting positions to N and extracting changed bases")
	baseToNull(seqarray,positionlist,basefile,indlist,verbose)
	
	if verbose==True:
		print("\tSeqarray again (1-5):")
		print(seqarray[0:5])
	
	print("\nWriting converted file")
	seqarrayToFasta(indlist,seqarray,outfile,verbose)
	
if __name__ == "__main__":
	main()