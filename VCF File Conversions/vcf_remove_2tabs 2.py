#!/usr/bin/env python3
## Made by Kaiya L Provost
## Last updated 4 November 2019

import sys
import gzip

try:
	infile = str(sys.argv[1])
	print("\tFile is: ",infile)
except:
	print("Filename not given, quitting")
	exit()

ending = infile[-3:]
outfile=infile+".header.vcf"

print("\treading")
if(ending == ".gz"):
	with gzip.open(infile,"rb") as input:
		text = input.read()
else:
	with open(infile,"rb") as input:
		text = input.read()
		
print("\treplacing")
if(ending == ".gz"):
	text = text.replace(b"\t\t",b"\t")
	text = text.replace(b"FORMATindiv0",b"FORMAT\tindiv0")
	text = text.replace(b" ",b"")
else:
	text = text.replace("\t\t","\t")
	text = text.replace("FORMATindiv0","FORMAT\tindiv0")
	text = text.replace(" ","")

print("\twriting")
if(ending == ".gz"):
	with gzip.open(outfile+".gz","wb") as output:
		output.write(text)
else:
	with open(outfile,"w") as output:
		output.write(text)

