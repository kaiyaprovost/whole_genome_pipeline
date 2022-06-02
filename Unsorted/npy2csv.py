import numpy as np
import os
import glob
import sys
import pickle
import gzip

os.chdir("/Users/kprovost/Downloads/")
for npyfile in glob.glob("*/*/*/*/*/*.npy",recursive=True):
	print(npyfile)
	try:
		df = np.load(npyfile,allow_pickle=True,fix_imports=True,encoding="latin1")
		np.savetxt(npyfile+".csv", df, delimiter="\t")
	except:
		print("PROBLEM WITH THIS FILE, SKIPPING")
		
	

