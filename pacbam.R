library(tidyverse)

wd <- '/BCGLAB/darosio_crispr/pacbam'

setwd(wd)

# create toy vcf

vcf <- data.frame(CHROM = 'egfp',
                  POS = 0,
                  REF = 'A',
                  ALT = 'G',
                  QUAL = '.',
                  FILTER = '.',
                  INFO = '.',
                  stringsAsFactors = FALSE)

vcf.file <- '/BCGLAB/darosio_crispr/pacbam/toy.vcf'

colnames(vcf)[1] <- '#CHROM' 

write.table(vcf,file = vcf.file,sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)

# create bed

fasta <- '/BCGLAB/darosio_crispr/reference/egfp.fasta'

fa <- readLines(fasta)

contig <- str_replace(string = paste(grep(fa,pattern = '^>',invert = FALSE,value = TRUE),collapse = ""),pattern = '>',replacement = '')

fa <- paste(grep(fa,pattern = '^>',invert = TRUE,value = TRUE),collapse = "")
fa <- str_to_upper(fa)

nchar(fa)

bed <- data.frame(contig = contig,
                  start = 0,
                  end = nchar(fa) - 1,
                  stringsAsFactors = F)

bed.file <- gsub(basename(fasta),pattern = '\\.fasta$',replacement = '.bed')

write.table(bed,file = bed.file,sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)

bam_files <- list.files('/BCGLAB/darosio_crispr/bam',full.names = TRUE,pattern = 'sorted.bam$')

for(bam in bam_files){
  
  message(bam)
  
  cmd <- paste0('/CIBIO/sharedRL/Projects/PaCBAM/git_repo/pacbam/pacbam',' bam=',bam,' fasta=',fasta,' vcf=',vcf.file,' bed=',file.path(wd,bed.file),' mode=4',' threads=20')
  
  system(cmd)
  
}



