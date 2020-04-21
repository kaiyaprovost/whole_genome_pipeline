import sys
import os
import glob
import pandas as pd

os.chdir("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/MISSING")

outfile="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/MISSING/missing_data_summarized.temp"

listspps=[]
listchrs=[]
listmiss=[]
listfile=[]

for imiss in glob.iglob("*/*imiss",recursive=False):
	filen=imiss.split("/")[-1]
	species=filen.split("-")[1]
	chr=filen.split("_")[-1].split(".")[0]
	#print(imiss)
	
	df = pd.read_table(imiss)
	totsites=sum(df.N_DATA)
	missings=sum(df.N_MISS)
	
	if totsites > 0:
		percentmiss=(missings/totsites)
		listmiss.append(percentmiss)
		listspps.append(species)
		listchrs.append(chr)
		listfile.append(imiss)
	
	## indiv, ndata, nfiltered, nmissing, fmissing
	

outdf = pd.DataFrame()
outdf['filename'] = listfile
outdf['species'] = listspps
outdf['chromosome'] = listchrs
outdf['percent_missing'] = listmiss

outdf.to_csv(outfile,sep="\t",index=False)
		


