#!/usr/bin/env python3

import sys
from collections import Counter

# cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/"
# for fasta in *.fasta; do
# python3 "/Users/kprovost/Documents/Github/whole_genome_pipeline/VCF File Conversions/fasta2ms.py" $fasta 1 1 1 
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

try:
	diploid_int = int(sys.argv[4])
	if diploid_int != "0":
		diploid = True
		print("Lines are diploid")
	else:
		diploid = False
		print("Lines are haploid")
except:
	print("Diploid not set, defaulting to haploid")
	diploid=False 

#fasta="/Users/kprovost/Dropbox (AMNH)/CFB_review_J_Biogeo/Crotalus_scutulatus_concat.fa"
msoutfile = fasta+".missing"+str(convertGapToMissing)+".monomorphic"+str(monomorphic)+".NOVCF.ms"
print("Outfile:",msoutfile)

def ambigExp(base):
	'''
	disambiguates the heterozygous ambiguity codes
	'''
	if base == "Y":
		expand = "CT"
	elif base == "R":
		expand = "AG"
	elif base == "W":
		expand = "AT"
	elif base == "S":
		expand = "GC"
	elif base == "K":
		expand = "TG"
	elif base == "M":
		expand = "CA"
	else:
		expand = base
	return expand

def ambigExpDiploid(dna_seq):
	first = ""
	second = ""
	for base in dna_seq:
		disambig = ambigExp(base)
		if len(disambig)>=2:
			first += disambig[0]
			second += disambig[1]
		else:
			first += disambig
			second += disambig
	return(first,second)

def fasta2reads_dict(fasta,verbose=False):
	with open(fasta,"r") as infile:
		lines = infile.readlines()
	## iterate through the lines, build dictionary
	#fasta_dict = {}
	reads_dict = {}
	## this loop reads each line, then it checks whether all reads are the same and makes them ambiguous 
	for line_num in range(len(lines)):
		if verbose==True:
			print(line_num)
		line = lines[line_num]
		#print(line)
		if line[0]==">":
			## is a name
			if diploid==False:
				ind_name = line.strip()[1:-2] ## takes off the ">", the read number and the _ before it
			else:
				ind_name = line.strip()[1:] ## takes off initial ">" only
			## check if the name is first (0) or second (1) read
		
			if diploid==False:
				read_number = (line.strip()[-1])
				if read_number == "a":
					read_number = 0
				if read_number == "b":
					read_number = 1
			else:
				read_number = 0
			if(int(read_number)==0):
				#if diploid==False:
				reads_dict[ind_name] = [None,None]
				#else:
				#	reads_dict[ind_name] = [None]
				#	## need to disambiguate
				#	
		else:
			dna_seq = line.strip().upper()
			#fasta_dict[ind_name][read_number] = dna_seq.upper()
			if int(read_number) == 0:
				if diploid==False:
					reads_dict[ind_name][0] = dna_seq
				else:
					first,second = ambigExpDiploid(dna_seq)
					reads_dict[ind_name][0] = first
					reads_dict[ind_name][1] = second
			elif int(read_number) == 1:
				reads_dict[ind_name][1] = dna_seq
			else:
				print("Warning: Something went wrong",line_num,read_number)
	return(reads_dict)

## ambiguity code text, from baseDisAmbig.py on github
reads_dict = fasta2reads_dict(fasta)

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

segsites = "//\nsegsites: "+str(len(snp_dict.keys()))
positions="\npositions: "+" ".join([str(x) for x in list(snp_dict.keys())])

snp_binary_dict={}

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
		try:
			alt_snps = list(set(alt_snps))
		except:
			print("PROBLEM WITH: "+str(snp_position))
			print(alt_snps)
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
		snp_binary_dict[snp_position] = readnums_nomissing

reads_binary_dict = {}

with open(msoutfile,"w") as outfile:
	_ = outfile.write(segsites)
	_ = outfile.write(positions)

for individual_number in range(len(snp_dict.keys())):
	print(individual_number)
	sequence=""
	for snp_position in snp_binary_dict.keys():
		sequence+=str(snp_binary_dict[snp_position][individual_number])
	with open(msoutfile,"a") as outfile:
		_ = outfile.write("\n"+sequence)


