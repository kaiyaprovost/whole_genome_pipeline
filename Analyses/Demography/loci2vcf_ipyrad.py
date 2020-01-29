#!/usr/bin/env python2

import time
import numpy as np
from ipyrad.assemble.util import *
from collections import Counter, OrderedDict
import itertools
import re

def most_common(L):
	return max(itertools.groupby(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]


#def make(data, samples):
def main():
	""" build a vcf file from the supercatg array and the cat.clust.gz output"""
	
	#print("generate outfile")
	#outfile = open(os.path.join(data.dirs.outfiles, data.name+".vcf"), 'w')
	outfile_file = "/Users/kprovost/Downloads/outfiles_allgroup/cardcard_files2/cardcard16_monomorphic.vcf"
	
	#print("read locifile")
	#inloci = os.path.join(data.dirs.outfiles, data.name+".loci")
	inloci = "/Users/kprovost/Downloads/outfiles_allgroup/cardcard_files1/cardcard16.loci"
	#inloci = "/Users/kprovost/Downloads/outfiles_allgroup/cardcard_files1/cardcard16_shorttest3.loci"
	
	#print("read samples")
	samples_file = "/Users/kprovost/Documents/Dissertation/ddRAD_individual_names.txt"
	#samples_file = "/Users/kprovost/Documents/Dissertation/ddRAD_individual_names_shortloci2.txt"
	
	#print("get names")
	#names = [i.name for i in samples]
	with open(samples_file,"rU") as samples:
		names = samples.readlines()
		names = [i.strip() for i in names]
	names.sort()
	#print(names)

	## TODO: Get a real version number for the current sw stack
	#version = "0.1"
	version = "klptesting"
	## TODO: This is just reporting minimum depth per base. Would it be useful to
	## report real depth of reads per base? YEAH, that's what supercatg is for.
	#mindepth = data.paramsdict["mindepth_statistical"] 
	mindepth = 6

	with open(outfile_file,"w") as outfile:

		print >>outfile, "##fileformat=VCFv4.1"
		print >>outfile, "##fileDate="+time.strftime("%Y%m%d")
		print >>outfile, "##source=ipyRAD.v."+version
		print >>outfile, "##reference=common_allele_at_each_locus"
		print >>outfile, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"
		print >>outfile, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
		print >>outfile, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"
		print >>outfile, "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">"
		print >>outfile, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
		print >>outfile, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"
		print >>outfile, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"
		print >>outfile, "\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+list(names))

		#loci = open(inloci).read().split("|")[:-1]
		#print(len(loci))
		loc = open(inloci).read()
		loci = re.split('\|\d+\|',loc)[:-1]
		#print(len(loci))
		#print(loci)
		snps = 0
		bases = 0 
		vcflist = []
		for locusnumber in range(len(loci)):
			#print("locusnum_1: "+str(locusnumber))
			samps = [i.split()[0][1:] for i in loci[locusnumber].strip().split("\n") if ">" in i]
			#print(samps[0:5])
			loc = np.array([tuple(i.split()[-1]) for i in loci[locusnumber].strip().split("\n") if ">" in i])
			#print(loc)
			NS = str(len(loc))
			DP = str(mindepth)
			for base in range(len(loc.T)):
				#print("---")
				col = []
				site = list(loc.T[base])
				#site1 = list("".join(site).replace("-","").replace("N",""))
				#print(site)
				site = list("".join(site).replace("-","").replace("N",""))
				#print(site)
				#print("----")
				
				#if not site: 
				#if not site1: 
					#print("no data at this base")
					#print(list(loc.T[base]))
				if site:
				#elif site:
				#elif site1:
					for bb in site:
						if bb in list("RKYSWM"):
							col += unstruct(bb)[0]
							col += unstruct(bb)[1]
						else:
							col += bb
					#print(col)
					REF = most_common([i for i in col if i not in list("-RKYSWMN")])
					ALT = set([i for i in col if (i in list("ATGC-N")) and (i!=REF)])
					
					
					#print("ref: "+str(REF)+"\talt: "+str(ALT))
					#print("ref")
					#print(REF)
					#print("alt")
					#print(ALT)
					if ALT:
						#print("not monomorphic")
						snps += 1
						bases += 1 ##
						GENO = [REF]+list(ALT)
						GENOS = []
						for samp in names:
							if samp in samps:
								idx = samps.index(samp)
								f = unstruct(loc.T[base][idx])
								if ('-' in f) or ('N' in f):
									GENOS.append("./.")
								else:
									GENOS.append(str(GENO.index(f[0]))+"|"+str(GENO.index(f[1])))
								#GENOS.append(str(GENO.index(f[0]))+"|"+str(GENO.index(f[1])))
							else:
								GENOS.append("./.")
						vcflist.append("\t".join([`locusnumber+1`, `base+1`, '.', REF, ",".join(ALT), "20", "PASS",
												  ";".join(["NS="+NS, "DP="+DP]), "GT"]+GENOS))
					else:
						#print("monomorphic")
						ALT = "."
						bases += 1
						GENO = [REF]+list(ALT)
						GENOS = []
						for samp in names:
							if samp in samps:
								idx = samps.index(samp)
								f = unstruct(loc.T[base][idx])
								if ('-' in f) or ('N' in f):
									GENOS.append("./.")
								else:
									GENOS.append(str("0|0"))
							else:
								GENOS.append("./.")
						vcflist.append("\t".join([`locusnumber+1`, `base+1`, '.', REF, ",".join(ALT), "20", "PASS",
												  ";".join(["NS="+NS, "DP="+DP]), "GT"]+GENOS))
						
			if not locusnumber % 1000:
				outfile.write( "\n".join(vcflist)+"\n" )
				vcflist = []
											  
						#print >>outfile, "\t".join([`locusnumber+1`, `base+1`, '.', REF, ",".join(ALT), "20", "PASS",
						#							";".join(["NS="+NS, "DP="+DP]), "GT"]+GENOS)
			if locusnumber % 100 == 0:
				print("locusnum: "+str(locusnumber))
	

		outfile.write( "\n".join(vcflist) )
	#outfile.close()

if __name__ == "__main__":
	#make(WORK, version, outname, mindepth, names)
	main()