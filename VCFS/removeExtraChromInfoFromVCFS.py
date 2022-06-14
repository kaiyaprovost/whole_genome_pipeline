import glob
import sys
import os
# for i in /Users/kprovost/Downloads/*.vcf
# do echo "${i}"
# python3 "/Users/kprovost/Documents/GitHub/whole_genome_pipeline/VCF File Conversions/removeExtraChromInfoFromVCFS.py" "${i}";
# done
def main():
	try:
		toconvert = sys.argv[1]
		print("\nRead the following file to convert chroms:")
		print(toconvert)
	except:
		print("\nNo file to convert chroms given, quitting")
		exit()
		
	try:
		namechange = sys.argv[2]
		print("\nRead the following file to convert names to:")
		print(namechange)
	except:
		print("\nNo file to convert names given")
	if os.path.isfile(toconvert+".fixedchroms"):
		print("Already converted")
	else:
		with open(toconvert,"r") as infile:
			infilelines = infile.readlines()
		with open(toconvert+".fixedchroms","w") as outfile:
			removeExtraInfo(infilelines,outfile)
		

def removeExtraInfo(infilelines,outfile):
	for i in range(len(infilelines)):
		line = infilelines[i]
		
		## check if teh line starts with # if so don't change it
		if line[0] == "#":
			newline=line
		else:
			split = line.split("\t")
			## needto change the 0th indexed thing
			tochange = split[0]
			nochange = split[1:]
			## PseudoNC_007897.1_Tgut_mtDNA
			splitchange = tochange.split("_")
			chromtype = splitchange[0] ## this is NW, PseudoNC, etc
			## if the chromtype isn't PseudoNC we don't want to print the line
			if (chromtype != "PseudoNC"):
				newline=None
			else:
				chromname = splitchange[3]
				newsplit = [chromtype,chromname]
				newline = "_".join(newsplit)+"\t"+"\t".join(nochange)
				
		if newline != None :
			outfile.write(newline)
if __name__ == "__main__":
	main()
