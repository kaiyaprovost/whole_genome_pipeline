file="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/Vireo-bellii.vcf"
out="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/Vireo-bellii_fixed.vcf"

with open(file,"r") as infile:
	with open(out,"a") as outfile:
		line = infile.readline()
		while line:
			#print(line)
			## need to split the lines by "\t" and then by "_"
			split = line.split("\t")
			split[0] = split[0].split("_")[-1]
			newline = "\t".join(split)
			_ = outfile.write(newline)
			line = infile.readline()
			
