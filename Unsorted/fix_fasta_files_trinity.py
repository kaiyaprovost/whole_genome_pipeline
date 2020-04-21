#!/usr/bin/env python3

import sys
import os
import glob
import copy

# for fasta in *fa; do python3 fix_fasta_files_trinity.py $fasta; done;

try:
	fasta = sys.argv[1]
	print("\nRead the following file to convert:")
	print(fasta)
	#fasta="/Users/kprovost/Dropbox (AMNH)/messed_up_assemblies/Assemblies_working/Prrli_308975AMNH-6.fa"
except:
	print("\nNo file to convert given, quitting")
	sys.exit()

with open(fasta,"r") as infile:
	lines=infile.readlines()

newlines = copy.deepcopy(lines)

count=0
for i in range(len(lines)) :
	line = lines[i]
	if line[0]==">":
		## it is header, must change
		count += 1
		newlines[i] = ">TR"+str(count)+"|c0_g1_i1\n"

with open(fasta+".trinity","w") as outfile:
	outfile.writelines(newlines)

