import numpy as np
import os
import glob
import sys
import pickle

os.chdir("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/qopts_1/")
for npyfile in glob.glob("*/*.npy"):
	print(npyfile)
	try:
		df = np.load(npyfile,allow_pickle=True,fix_imports=True,encoding="latin1")
		np.savetxt(npyfile+".csv", df, delimiter="\t")
	except:
		print("PROBLEM WITH THIS FILE, SKIPPING")
		
	

