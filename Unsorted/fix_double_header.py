#double header fixer python script

import sys
import glob
import os
import gzip

path="/home/kprovost/nas5/slim_osg/subdirectory/DOUBLEHEADER/"

os.chdir(path)
listfiles=glob.glob("*.vcf*")

for file in listfiles:
	#file="migrate-0.001-popsize-100000-mut-2.21e-10-gen-1000000-recom-1e-7-ibd-0-seccon-0-pop-2-1583563789-1-recap_2.21e-10-1e-07-100000-0.02.vcf_bases.vcf"
	print(file)
	
	suffix=file[-3:]
	if os.path.exists(path+file+"FIX.vcf"):
		print("skipping")
	elif os.path.exists(path+file+"FIX.vcf.gz"):
		print("skipping")
	else:
		if suffix==".gz":
			with gzip.open(path+file,"rb") as infile:
				lines=infile.readlines()
		else:
			with open(path+file,"r") as infile:
				lines=infile.readlines()	
		if len(lines) <= 12:
			print("PROBLEM -- skipping, too few lines")
		else:
			linelendicti = {}
			prev=None
			breakpoints=[]
			lenbreaks=[]
			differences=[0]
			for i in range(len(lines)):
				line = lines[i]
				if suffix==".gz":
					thislength=len(line.split(b"\t"))
				else:
					thislength=len(line.split("\t"))
				entry = linelendicti.get(thislength,None)
				if entry==None:
					linelendicti[thislength] = 1
				else:
					linelendicti[thislength] += 1
				if prev != thislength or i == len(lines)-1:
					if i != 0:
						differences.append(i-breakpoints[-1])
					breakpoints.append(i)
					lenbreaks.append(thislength)
					prev = thislength
			end=differences.index(max(differences))
			headerstarts = [i for i, value in enumerate(lenbreaks) if value == 1] 
			endline = breakpoints[end]+1
			for j in reversed(headerstarts):
				if j < end:
					start=j
					break
			startline=breakpoints[start]
			newlines=lines[startline:endline]
			if suffix==".gz":
				with gzip.open(path+file+"FIX.vcf.gz","wb") as outfile:
					outfile.writelines(newlines)
			else:
				with open(path+file+"FIX.vcf","w") as outfile:
					outfile.writelines(newlines)

print("done")












