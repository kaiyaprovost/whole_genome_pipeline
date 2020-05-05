import glob
import sys
import os


# python3 "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/angsd_ind0_to_cat.py" \
# "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/BILINEATA/BILINEATA_distancematrix_FULLGENOME.csv" \
# "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/Amphispiza-bilineata.bamlist"

def main():
	try:
		toconvert = sys.argv[1]
		print("\nRead the following file to convert:")
		print(toconvert)
	except:
		print("\nNo file to convert given, defaulting to a BELLII file")
		toconvert = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/BELLII/BELLII_distancematrix_FIRSTNGS.csv"
		print(toconvert)
	fileend = toconvert[-3:]
	try:
		bamlist = sys.argv[2]
		print("\nRead the following bamlist:")
		print(bamlist)
	except:
		print("\nNo bamlist given, defaulting")
		bamlist = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ASSEMBLY/ANGSD/A5.bamlists/bamlist/Vireo-bellii.bamlist"
		print(bamlist)
	try:
		rg = sys.argv[3]
		print("\nRead the following RG list:")
		print(rg)
	except:
		print("\nNo RG list given, defaulting")
		rg = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/rapid_genomics_by_catalog.txt"
		print(rg)
	if os.path.exists(toconvert+".converted"):
		print("ALREADY DONE")
	else:
		with open(bamlist,"r") as bamfile:
			bamlines = bamfile.readlines()
		iddict = {}
		for i in range(len(bamlines)-1,-1,-1):
			bam = bamlines[i]
			file = bam.strip().split("/")[-1]
			id = file.split(".")[0][4:]
			if(fileend=="csv"):
				lookup = "Ind_"+str(i)
			else:
				lookup = "indiv"+str(i)
			iddict[lookup] = id
		with open(toconvert,"r") as convertfile:
			converttext = convertfile.read()
		for lookup,id in iddict.items():
			#print(lookup+" "+id)
			converttext = converttext.replace(lookup,id)
		with open(rg,"r") as rgfile:
			rglines = rgfile.readlines()
		rgdict = {}
		for j in range(len(rglines)):
			rgline = rglines[j]
			rgid,rglook = rgline.strip().split("\t")
			rgdict[rglook] = rgid
		for rglook,rgid in rgdict.items():
			#print(rglook+" "+rgid)
			converttext = converttext.replace(rglook,rgid)
		with open(toconvert+".converted","w") as outfile:
			outfile.write(converttext)

if __name__ == "__main__":
	main()