# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 00:21:52 2016

edited by K Provost
"""

#Aligning sequences
#Muscle software installed required: http://www.drive5.com/muscle/downloads.htm    
def getAlignTree(filename,treepath,outpath,muscle_exe):
    from Bio.Align.Applications import MuscleCommandline
    from Bio import AlignIO
    #muscle_exe = "/Users/kprovost/Documents/Publications/Parrots/ParrotPipelineRedo/SCRIPTS/muscle3.8.31_i86darwin64"
    outname = "ALIGNED_"+filename
    print("GET TREE: "+filename)
    with open(filename,"r") as infile:
        read = infile.read()
        count = read.count(">")
        if count <= 2: 
            treeName = "None"
            print("ONLY ONE SEQ, DONE")               
        else:
            treeName = outname+"_tree.tre"
            try:
                muscle_cline = MuscleCommandline(muscle_exe,input=filename,tree1=treeName)
                stdout, stderr = muscle_cline()
                print("\nTREE DONE")
            except:
                print("??? ERROR")
                print(treeName)
    return(treeName)
    
def align(filename,treepath,outpath,muscle_exe):
    from Bio.Align.Applications import MuscleCommandline
    from Bio import AlignIO
    #muscle_exe = "/Users/kprovost/Documents/Publications/Parrots/ParrotPipelineRedo/SCRIPTS/muscle3.8.31_i86darwin64"
    outname = "ALIGNED_"+filename
    print("ALIGNING: "+filename)
    with open(filename,"r") as infile:
        read = infile.read()
        count = read.count(">")
        if count <= 1: 
            with open(outpath+outname,"w") as outfile:
                outfile.write(read)
                print("ONLY ONE SEQ, DONE")
        else:
            try:
                muscle_cline = MuscleCommandline(muscle_exe, input=filename, out=outname)
                stdout, stderr = muscle_cline()
                AlignIO.read(outname, "fasta")
                print("ALIGNED")  
            except:
                print("??? ERROR")
                prrint(filename)  
    
def findRootOfTree(treeName,logFile,treepath,outpath,muscle_exe):
    import shutil
    goodtaxa = []
    with open(treeName,"r") as infile,open(treepath+logFile,"a") as outfile:
        #outfile.write("\n###NEW RUN###\n")
        lines = infile.readlines()
        root1index = int((len(lines)-2)-1)
        rootNum = float(lines[root1index][2:])
        #print(rootNum)
        root2index = int(((len(lines)-2)/2)-1)
        print(lines[root1index].strip(),lines[root2index].strip())
        if lines[root2index].strip()[2:].strip() == str(rootNum):
            #print("yay!")
            firstclade = "\n".join(lines[0:root2index])
            firsttaxa = firstclade.count("///")
            firstrev = firstclade.count("~~~REV")
            lastclade = "\n".join(lines[root2index:])
            lasttaxa = lastclade.count("///")
            lastrev = lastclade.count("~~~REV")
            if firsttaxa == lasttaxa:
                #print("yay2!")
                print("First rev:",firstrev,"Last rev:",lastrev,"Taxa:",firsttaxa)
                if firstrev == 0 or firstrev == firsttaxa:
                    outfile.write("NOREV,"+treeName+"\n")
                else:
                    print("#####")
                    outfile.write(str(firstrev)+":"+str(lastrev)+"MISMATCH,"+treeName+"\n")
                clade = firstclade.split("\n")
                for line in clade:
                    #print(line)
                    #print(line[:1])
                    if line[:1] in [",","(",")",";",":","","\n"]:
                        pass
                    else:
                        identifier = line.split(":")[0]
                        goodtaxa.append(identifier)
                ## get the 
            
            else:
                print("boo2!")
                outfile.write("MISSINGTAXA,"+treeName+"\n")
        else:
            print("boo!")
            outfile.write("WEIRDROOT,"+treeName+"\n")   
        return(goodtaxa)        

def subsetFastaFromList(filename,goodtaxa,treepath,outpath,muscle_exe):
    with open(filename,"r") as infile,open(outpath+"READY_"+filename,"w") as outfile:
        read = infile.read()
        fasta = read.split(">",1)[1].split("\n>")
        for line in fasta:
            name = line.split("\n")[0]
            #print(name)
            if name in goodtaxa:
                #print(name)
                outfile.write(">"+line+"\n")
            else:
                #print("NO")
                pass
        
def main():
    import os
    import sys
    import glob
    import shutil
   
    cwd = os.getcwd()
   
    try:
        muscle_exe = sys.argv[1]
        print("\tMuscle path exists")
    except:
        #print("Muscle defaulting to:")
        #print("/Users/kprovost/Documents/Publications/Parrots/ParrotPipelineRedo/SCRIPTS/muscle3.8.31_i86darwin64")
        print("Path to Muscle not given. Please give path to Muscle.")
        quit()

    try:
        path = sys.argv[2]
        print("\tPath is: ",path)
    except:
        print("Path to FASTA files with reverse complements not given. Please give path.")
        #path = os.getcwd()+"/withrev/"
        quit()

    treepath = cwd+"/alignmenttrees/"
    if not os.path.exists(treepath):
        print("creating folder: ",treepath)
        os.makedirs(treepath)   
        
    outpath = cwd+"/readytoalign/"
    if not os.path.exists(outpath):
        print("creating folder: ",outpath)
        os.makedirs(outpath)    

    os.chdir(path)
    logFile = "treeLog.txt"    
    for filename in glob.glob("*.fa*"):
        print("###############################")
        #print(filename)
        treeName = getAlignTree(filename,treepath,outpath,muscle_exe)
        if treeName != "None":
            goodtaxa = findRootOfTree(treeName,logFile,treepath,outpath,muscle_exe)
            if goodtaxa == "ALL":
                print("\t",filename," IS GOOD\n\n")
                shutil.copy2(filename,outpath+"READY_"+filename)
            elif goodtaxa == []:
                print("ERROR")
            else:
                print(goodtaxa)
                subsetFastaFromList(filename,goodtaxa,treepath,outpath,muscle_exe)
                
    print("\n\nDONE")
                
    
if __name__ == "__main__":
    main()