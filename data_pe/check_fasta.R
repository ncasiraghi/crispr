library(stringr)

fasta <- readLines('/Users/ncasiraghi/Documents/unitn/darosio/reference/bwa_reference/spike.fasta')

fasta <- paste(grep(fasta,pattern = '^>',invert = TRUE,value = TRUE),collapse = "")

fasta <- str_to_upper(fasta)

nchar(fasta)

# cut site

cts <- str_locate_all(string = fasta,pattern = 'CCCACCG') 

unlist(cts)[1] - 1
