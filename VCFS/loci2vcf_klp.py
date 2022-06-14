#!/usr/bin/env python2

def main():
	infile = "/Users/kprovost/Downloads/outfiles_allgroup/cardcard_files2/cardcard16_shorttest.loci"
	indfile = "/Users/kprovost/Documents/Dissertation/ddRAD_individual_names.txt"
	
	tempDict = {}
	locusCount = 0
	currentInds = []
	
	with open(infile,"rU") as locifile:
		for line in locifile:
			#print(line)
						
			if line[0] == ">":
				print("add to dictionary")
				
				ind,seq = line.split()
				#print(ind)
				#print(seq[0:5])
				tempDict[ind] = seq
				
			elif line[0:2] == "//":
				print("convert dictionary to vcf")
				print(tempDict)
				tempDict = {}
				count += 1
				
			else:
				print("error -- not recognized")
				print(line)
				quit()
				



if __name__ == "__main__":
	main()
	
