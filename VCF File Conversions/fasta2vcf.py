#!/usr/bin/env python3

import sys
from collections import Counter

# for fasta in /Users/kprovost/Dropbox\ \(AMNH\)/CFB_review_J_Biogeo/*fa; do
# echo $fasta;
# python3 /Users/kprovost/Documents/Github/whole_genome_pipeline/fasta2vcf.py "$fasta" 0 1;
# python3 /Users/kprovost/Documents/Github/whole_genome_pipeline/fasta2vcf.py "$fasta" 1 1;
# python3 /Users/kprovost/Documents/Github/whole_genome_pipeline/fasta2vcf.py "$fasta" 0 0;
# python3 /Users/kprovost/Documents/Github/whole_genome_pipeline/fasta2vcf.py "$fasta" 1 0;
# done 

try:
	fasta = str(sys.argv[1])
	print("\tFasta file is: ",fasta)
except:
	print("Fasta file not given, quitting")
	exit()

try:
	keepmissing_int = str(sys.argv[2])
	if keepmissing_int != "0":
		convertGapToMissing = True
		print("Converting gaps to missing data")
	else:
		convertGapToMissing = False
		print("Retaining gaps as SNPS")
except:
	print("Whether to keep gaps not specified. Retaining gaps as SNPS.")
	convertGapToMissing = False

try:
	monomorphic_int = str(sys.argv[3])
	if monomorphic_int != "0":
		monomorphic = True
		print("Retaining monomorphic sites")
	else:
		monomorphic = False
		print("Dropping monomorphic sites")
except:
	print("Whether to keep monomorphic sites not specified. Dropping monomorphic sites.")
	monomorphic = False

#fasta="/Users/kprovost/Dropbox (AMNH)/CFB_review_J_Biogeo/Crotalus_scutulatus_concat.fa"
vcfoutfile = fasta+".missing"+str(convertGapToMissing)+".monomorphic"+str(monomorphic)+".generated.vcf"
print("Outfile:",vcfoutfile)

def fasta2reads_dict(fasta,verbose=False):
	with open(fasta,"r") as infile:
		lines = infile.readlines()
	## iterate through the lines, build dictionary
	fasta_dict = {}
	reads_dict = {}
	## this loop reads each line, then it checks whether all reads are the same and makes them ambiguous 
	for line_num in range(len(lines)):
		if verbose==True:
			print(line_num)
		line = lines[line_num]
		#print(line)
		if line[0]==">":
			## is a name
			ind_name = line.strip()[1:-2] ## takes off the ">", the read number and the _ before it
			## check if the name is first (0) or second (1) read
			read_number = (line.strip()[-1])
			if read_number == "a":
				read_number = 0
			if read_number == "b":
				read_number = 1
			if(int(read_number)==0):
				reads_dict[ind_name] = [None,None]
		else:
			dna_seq = line.strip().upper()
			#fasta_dict[ind_name][read_number] = dna_seq.upper()
			if int(read_number) == 0:
				reads_dict[ind_name][0] = dna_seq
			elif int(read_number) == 1:
				reads_dict[ind_name][1] = dna_seq
			else:
				print("Warning: Something went wrong",line_num,read_number)
	return(reads_dict)

## ambiguity code text, from baseDisAmbig.py on github

def ambigSimp(dataList): 
	'''	reambiguates lists of bases into their ambiguity codes'''
	if "A" in dataList:
		hasA = True
	else:
		hasA = False
	if "C" in dataList:
		hasC = True
	else:
		hasC = False
	if "G" in dataList:
		hasG = True
	else:
		hasG = False
	if "T" in dataList:
		hasT = True
	else:
		hasT = False
	if hasA:
		if hasG:
			if hasC:
				if hasT:
					simp = "N"
				else:
					simp = "V"
			else:
				if hasT:
					simp = "D"
				else:
					simp = "R"
		else:
			if hasC:
				if hasT:
					simp = "H"
				else:
					simp = "M"
			else:
				if hasT:
					simp = "W"
				else:
					simp = "A"
	else:
		if hasG:
			if hasC:
				if hasT:
					simp = "B"
				else:
					simp = "S"
			else:
				if hasT:
					simp = "K"
				else:
					simp = "G"
		else:
			if hasC:
				if hasT:
					simp = "Y"
				else:
					simp = "C"
			else:
				if hasT:
					simp = "T"
				else:
					simp = dataList
	return(simp)

## find snps from fasta

reads_dict = fasta2reads_dict(fasta)

header = "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"+"\t"

individual_names=reads_dict.keys()

header=header+"\t".join(individual_names)+"\n"

#fasta_dict,reads_dict = fasta2fasta_dict(fasta)

def reads_dict2snps(reads_dict,convertGapToMissing=False,verbose=False,monomorphic=False):
	snp_dict = {}
	unique_dna=list(set([item for sublist in reads_dict.values() for item in sublist]))
	if len(unique_dna) >=2:
		for base_num in range(len(unique_dna[0])):
			if verbose == True:
				print(base_num)
			all_bases = [x[base_num] for x in [item for sublist in reads_dict.values() for item in sublist]]
			unique_bases = set(all_bases)
			unique_bases = [x for x in unique_bases if x != "?"] ## ignores unknowns
			unique_bases = [x for x in unique_bases if x != "N"] ## ignores unknowns
			if(convertGapToMissing==True):
				unique_bases = [x for x in unique_bases if x != "-"] ## this is if you want to ignore missing data when calling snps
			if monomorphic==False:
				if 	len(unique_bases) >= 2: 
					## there is a snp here
					snp_dict[base_num] = all_bases
			else:
				snp_dict[base_num] = all_bases
				
	return(snp_dict)

snp_dict = reads_dict2snps(reads_dict,convertGapToMissing=convertGapToMissing,monomorphic=monomorphic)

#snp_dict_nomissing = reads_dict2snps(reads_dict,keepMissing=False)

#overlap = set(snp_dict.keys()).intersection(snp_dict_nomissing.keys())

## now need to convert the snp dict, assign the "ref" as the most common and the alts as the others

## write out header

## if there are no alternative alleles then "." should be used

with open(vcfoutfile,"w") as outfile:
	_ = outfile.write(header) ## prevents number of characters from getting written to screen

for snp_position in snp_dict.keys():
	if snp_position % 10000 == 0:
		print(snp_position)
	reads = snp_dict[snp_position]
	reads = [r.replace("?",".") for r in reads]
	reads = [r.replace("N",".") for r in reads]
	if convertGapToMissing == True:
		reads = [r.replace("-",".") for r in reads]
	## check if there are non-"." reads
	if len(list(set(reads))) == 0:
		keepgoing = False
	elif len(list(set(reads))) == 1:
		if list(set(reads))[0] == ".":
			keepgoing = False
		else:
			keepgoing = True
	else:
		keepgoing = True
	if keepgoing == True:
		count_dict=dict(Counter(reads)) ## generates another dictionary with the values
		##Counter({'C': 4, '-': 2})
		maximum_count=max(count_dict.values())
		snps_with_max = [key for (key, value) in count_dict.items() if value == maximum_count]	
		alt_snps = [key for (key, value) in count_dict.items() if value != maximum_count]
		alt_snps.sort()
		if len(snps_with_max) >=2:
			## multiple maxes
			snps_with_max.sort() 
			if snps_with_max[0]=="-" or snps_with_max[0]==".":
				ref_snp = snps_with_max[1]
				alt_snps.append(snps_with_max[0])
				if len(snps_with_max) > 2:
					alt_snps.append(snps_with_max[2:])
			else:
				ref_snp = snps_with_max[0]
				if len(snps_with_max) > 2:
					alt_snps.append(snps_with_max[1:])
				else:
					alt_snps.append("".join(snps_with_max[1:]))
		else:
			## need to check if the ref snp is "-"
			if snps_with_max == "-" or snps_with_max[0]==".":
				ref_snp = alt_snps[0]
				alt_snps[0] = "".join(snps_with_max)
			else:
				ref_snp = snps_with_max
		alt_snps = list(set(alt_snps))
		alt_snps.sort(reverse=True)
		#ref_dict[snp_position] = ref_snp
		#alts_dict[snp_position] = alt_snps
		#allsnp = ref_dict[snp_position]+alts_dict[snp_position]
		allsnp = list(ref_snp)+list(alt_snps)
		readnums = [ allsnp.index(x) for x in reads ]
		readnums_nomissing = [ allsnp.index(x) for x in reads ]
		if "." in allsnp:
			gap_index = allsnp.index(".")
			readnums_nomissing = [str(r).replace(str(gap_index),".") for r in readnums_nomissing]
		if "." in alt_snps and len(alt_snps) >= 2:
			alt_snps.remove(".")
		if len(alt_snps) == 0:
			if alt_snps = ["."]
		toprint = "chr\t"+str(snp_position+1)+"\t.\t"+str("".join(ref_snp))+"\t"+",".join(alt_snps)+"\t.\tPASS\tGT"
		## phased = "|"
		for ind_num in range(0,len(readnums),2):
			if convertGapToMissing == False:
				genotype=str(readnums[ind_num])+"|"+str(readnums[ind_num+1])
			else:
				genotype=str(readnums_nomissing[ind_num])+"|"+str(readnums_nomissing[ind_num+1])
			toprint+=("\t"+genotype)
		toprint+="\n"
		#print(toprint)
		with open(vcfoutfile,"a") as outfile:
			_ = outfile.write(toprint)
	


