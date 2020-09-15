import sys
import os
import glob 

try:
	filepath = sys.argv[1]
	print("\tFile is: ",filepath)
except:
	print("Filename not given, quitting")
	#filepath = "/Users/kprovost/Documents/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/genomeresequencingFromLucas/for_AMN_245109/sedtest.txt"
	exit()

split = filepath.split("/")
filename = split[-1]
print(split)

splitfile = filename.split(".")
prefix = splitfile[0]
print(splitfile)
print(prefix)

to_replace = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t*"
replacement = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+prefix

# Read in the file
with open(filepath, 'r') as file :
  filedata = file.read()
  print("read")

# Replace the target string
filedata = filedata.replace(to_replace, replacement)
print("replaced 1")

to_replace = '##fileformat=VCFv4.0\n#CHROM'
replacement = '##fileformat=VCFv4.0\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM'

# Replace the target string
filedata = filedata.replace(to_replace, replacement)
print("replaced 2")

# Write the file out again
with open(filepath, 'w') as file:
  file.write(filedata)
  print("wrote")