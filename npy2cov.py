import numpy as np
import os
import glob
import sys
import pickle
import gzip

## this will work with npy files and npz files 

try:
    npyfile = sys.argv[1]
    print("\tNpy file is: ",npyfile)
except:
    print("tNpy file not given, quitting")
    exit()

#for npyfile in glob.glob("*cov.npy",recursive=True):
print(npyfile)
try:
	df = np.load(npyfile,allow_pickle=True,fix_imports=True,encoding="latin1")
	np.savetxt(npyfile.replace(".npy",""), df, delimiter="\t")
except:
	print("PROBLEM WITH THIS FILE, SKIPPING")
