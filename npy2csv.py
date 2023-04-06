import numpy as np
import os
import glob
import sys
import pickle
import gzip

## this will work with npy files and npz files 

os.chdir("/Users/kprovost/Documents/Postdoc_Working/Finished_Models/9SppBalanced/")
for npyfile in glob.glob("*.npy",recursive=True):
	print(npyfile)
	try:
		df = np.load(npyfile,allow_pickle=True,fix_imports=True,encoding="latin1")
		np.savetxt(npyfile+".csv", df, delimiter="\t")
	except:
		print("PROBLEM WITH THIS FILE, SKIPPING")
for npzfile in glob.glob("*.npz",recursive=True):
	print(npzfile)
	try:
		dfz = np.load(npzfile,allow_pickle=True,fix_imports=True,encoding="latin1")
		df=dfz['arr_0']
		np.savetxt(npzfile+".csv", df, delimiter="\t")
	except:
		print("PROBLEM WITH THIS FILE, SKIPPING")
	

