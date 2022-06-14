import glob
import sys
import os

# for i in /Users/kprovost/Dropbox\ \(AMNH\)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/NGSDIST/*ngsdist_2022.txt; do
# python3 "/Users/kprovost/Documents/GitHub/whole_genome_pipeline/Unsorted/angsd_ngsdistind_to_catalog.py" "$i"
# gzip "$i"
# done 

## python3  "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/Amphispiza_bilineata_50.bamlist.beagle.gz_ngsdist_2022.txt"

def main():
    try:
        toconvert = sys.argv[1]
        print("\nRead the following file to convert:")
        print(toconvert)
    except:
        print("\nNo file to convert given, quitting")
        exit()
        toconvert="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/Amphispiza_bilineata_50.bamlist.beagle.gz_ngsdist_2022.txt"
    try:
        bamlistpath = sys.argv[2]
        print("\nRead the following bamlist path:")
        print(bamlistpath)
    except:
        print("\nNo bamlist path given, defaulting")
        bamlistpath = "/Users/kprovost/Downloads/A5.bamlists/"
        print(bamlistpath)
    try:
        rg = sys.argv[3]
        print("\nRead the following RG list:")
        print(rg)
    except:
        print("\nNo RG list given, defaulting")
        rg = "/Users/kprovost/Downloads/rapid_genomics_by_catalog.txt"
        print(rg)
    if os.path.exists(toconvert+".converted"):
        print("ALREADY DONE")
    else:
        ## find the correct bamfile from the toconvert file
        toconvert_basename = os.path.basename(toconvert)
        if os.path.isdir(bamlistpath):
            os.chdir(bamlistpath)
            found_bamfile = toconvert_basename.split(".beagle.gz")[0]
            found_bamfile = bamlistpath+found_bamfile
        else:
            found_bamfile = bamlistpath
        ## check if file exists
        if os.path.exists(found_bamfile):
            with open(found_bamfile,"r") as bamfile:
                bamlines = bamfile.readlines()
            iddict = {}
            for i in range(len(bamlines)-1,-1,-1):
                bam = bamlines[i]
                file = bam.strip().split("/")[-1]
                id = file.split(".")[0][4:]
                lookup = "Ind_"+str(i)
                iddict[lookup] = id
            with open(toconvert,"r") as convertfile:
                converttext = convertfile.read()
            ## get the number of individuals from the first non-blank line
            convertsplit = converttext.split("\n")
            convertsplit = [i for i in convertsplit if i] ## removes blank lines
            ## first line is the number of individuals
            numinds=int(convertsplit[0])
            indrange = range(numinds)
            indseq = ["Ind_"+str(i) for i in indrange]
            convertsplit[0]="\t"+"\t".join(indseq)
            converttext="\n".join(convertsplit)
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
                converttext = converttext.replace(rglook,rgid.replace(" ","").replace("-",""))
            with open(toconvert+".converted","w") as outfile:
                _ = outfile.write(converttext)
        else:
            print("Bamfile not found:"+found_bamfile)

if __name__ == "__main__":
    main()