#!/bin/python3

import os
import sys
import gzip
import glob

def check_duplicate_zip_unzip(notzippedfile,folderpath=None):
	zippedfile = notzippedfile+".gz"
	
	if (os.path.exists(notzippedfile) and os.path.exists(zippedfile)) == False:
		print("Both files do not exist, exiting.")
		#raise FileExistsError('Both files do not exist.')
	else:	
		if folderpath == None:
			folderpath = os.path.dirname(notzippedfile)
		
		if folderpath == "":
			folderpath=os.getcwd()
		
		os.chdir(folderpath)
	
		oldfilesfolderpath = folderpath+"/ZIPCHECK_OLD/"
		if os.path.exists(oldfilesfolderpath) == False:
			os.mkdir(oldfilesfolderpath)
	
		notzip_basename= os.path.basename(notzippedfile)
		zip_basename= os.path.basename(zippedfile)
	
		notzip_time = os.path.getmtime(notzippedfile)
		zip_time = os.path.getmtime(zippedfile)
	
		if notzip_time > zip_time:
			## the unzipped file is newer and should be used
			## delete or move the zipped file
			os.rename(zippedfile, oldfilesfolderpath+zip_basename)
		else:
			## the unzipped file is older or has been modified at the same time
			## delete or move the unzipped file
			os.rename(notzippedfile, oldfilesfolderpath+notzip_basename)
	

try:
	notzippedfile = sys.argv[1]
	print("Read file: "+notzippedfile)
except:
	print("No file given, exiting")
	exit()
	#notzippedfile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/VCFS/Vireo-bellii-called.geno.PseudoNC_011484.1_Tgut_20.vcf.fixedchroms.converted_w100000_o20000_289.window.vcf"

if os.path.isdir(notzippedfile):  
	print("\nInput is a directory")
	folderpath=notzippedfile
	notzipped_filelist=glob.glob(pathname=folderpath+"/*.vcf")
	numfiles=len(notzipped_filelist)
	print("Read in "+str(numfiles)+" files")
	for file_index in range(len(notzipped_filelist)):
		thisfile = notzipped_filelist[file_index]
		#print(thisfile)
		check_duplicate_zip_unzip(thisfile,folderpath)
		print("Done "+str(file_index)+"/"+str(numfiles))
elif os.path.isfile(notzippedfile):  
	print("\nInput is a normal file")
	check_duplicate_zip_unzip(notzippedfile)
else:  
	print("\nCannot process this kind of input! Exiting")
	exit()

