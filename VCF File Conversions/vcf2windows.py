## python script to create windows from a vcf 

import numpy as np
import sys, os, glob
import copy

# cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/VCFS/"
# for i in */*converted; do
# python3 "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/scripts/vcf2windows.py" "${i}" 100000 20000
# done

try:
    vcffile = sys.argv[1]
    print("\tFile is: ",vcffile)
except:
    print("Filename not given, quitting")
    exit()
    
try:
    windowsize = int(sys.argv[2])
    print("\tWindow size is: ",windowsize)
except:
    print("Window size not given, quitting")
    exit()

try:
    overlap = int(sys.argv[3])
    print("\tOverlap is: ",overlap)
except:
    print("Overlap not given, quitting")
    exit()

suffix = "_w"+str(windowsize)+"_o"+str(overlap)

with open(vcffile,"r") as vcf:
	lines = vcf.readlines()

lastline = lines[-1]
lastpos = int(lastline.split("\t")[1])

starts = list(range(1,lastpos,overlap))
ends = [i + (overlap-1) for i in starts] 

window_index = 0

if lastpos <= windowsize:
	print("Don't need to convert! Already a window.")
	outfile_string = vcffile+suffix+"_"+str(window_index)+".window.vcf"
	with open(outfile_string,"w") as outfile:
		_ = outfile.write("".join(lines))


toprint = ""
header = ""

for line in lines:
	if line[0]=="#":
		toprint += line
		header += line
	else:
		this_start = starts[window_index]
		this_end = ends[window_index]
		this_pos = int(line.split("\t")[1])
		if this_pos >= this_start and this_pos <= this_end:
			toprint += line
		else:
			outfile_string = vcffile+suffix+"_"+str(window_index)+".window.vcf"
			with open(outfile_string,"w") as outfile:
				_ = outfile.write("".join(toprint))
			toprint = copy.deepcopy(header)
			window_index += 1
			if window_index % 100 == 0:
				print(str(window_index)+" of "+str(len(starts)))



# cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/VCFS/BELLII/"
# 
# for i in *converted; do
# echo $i
# temp=$i.temp
# outsamp1=`head -2 "${i}" |tail -1 |tr '\t' '\n' |wc -l`
# fixcols=9
# outsamp=`echo "$(($outsamp1-$fixcols))"`
# echo $outsamp
# perl "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/vcf2MS.pl" "${i}" "${temp}" ${outsamp}
# done
# 
# ## 