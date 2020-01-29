#!/usr/bin/env python3

import sys

try:
	infile = str(sys.argv[1])
	print("\tFile is: ",infile)
except:
	print("Filename not given, quitting")
	exit()
	
with open(infile,"r") as input:
	text = input.read()

text = text.replace("\t\t","\t")
text = text.replace("FORMATindiv0","FORMAT\tindiv0")

with open(infile,"w") as output:
	output.write(text)

