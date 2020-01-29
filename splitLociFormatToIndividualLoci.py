import os

with open("/Users/kprovost/Documents/cardcard16_labeled.loci","r") as locifile:
	loci = locifile.readlines()
	if not os.path.exists("/Users/kprovost/Documents/cardcard16_labeled_indloci/"):
		print("creating folder: ","/Users/kprovost/Documents/cardcard16_labeled_indloci/")
		os.makedirs("/Users/kprovost/Documents/cardcard16_labeled_indloci/")

count = 0
print(str(count).zfill(5))
for line in loci:
	with open("/Users/kprovost/Documents/cardcard16_labeled_indloci/locus_"+str(count).zfill(5)+".fasta","a") as outfile:
		if line[0] == ">":
			ind,seq = line.split(" ",1)
			ind = ind.strip()
			seq=seq.strip()
			outfile.write(ind+"\n"+seq+"\n")
		elif line[0] == "/":
			count += 1 
			print(str(count).zfill(5))
		else:
			print(line)
			


