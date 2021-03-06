library(stringr)

fasta <- readLines('/BCGLAB/darosio_crispr/reference/bwa_reference/egfp.fasta')

fasta <- paste(grep(fasta,pattern = '^>',invert = TRUE,value = TRUE),collapse = "")

fasta <- str_to_upper(fasta)

nchar(fasta)

# amplicone 

ampl <- str_locate_all(string = fasta,pattern = 'GGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAAC')

ampl

# L-JON

jon <- str_locate_all(string = fasta,pattern = 'AGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCC')

jon

str_locate_all(string = fasta,pattern = 'ACCGGCAA')

# L-16

gfps16 <- str_locate_all(string = fasta,pattern = 'CTGACCTACGGCGTGCAGT')

gfps16

str_locate_all(string = fasta,pattern = 'TACGGCGT')
