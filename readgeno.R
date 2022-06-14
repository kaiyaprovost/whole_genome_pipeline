#library(LEA)

setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/")
files = list.files(pattern=".vcf")

for (input.file in files[20]) {

lines=readLines(input.file,n=2)
writeLines(lines,headerfile)

}

for (input.file in files) {
  
  print(input.file)
  headerfile=paste(input.file,".header",sep="")
  header = readLines(headerfile)
  
  lines=readLines(input.file)
  
  lines[1:2] = header
  lines=lines[-1]
  
  writeLines(lines,input.file)
  
}
