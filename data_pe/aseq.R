library(tidyverse)

# compute pileup

wd <- '/mnt/profile/darosio'

fasta <- file.path(wd,'reference/bwa_reference/spike.fasta')

bfs <- c("/mnt/profile/darosio/data/bowtie2",
         "/mnt/profile/darosio/data/bwa-aln-sampe",
         "/mnt/profile/darosio/data/bwa-mem")

for(bam.folder in bfs){
  
  message(bam.folder)
  
  if(!exists(file.path(bam.folder,'pileup.aseq'))){
    dir.create(file.path(bam.folder,'pileup.aseq'),showWarnings = FALSE)
  }
  
  fa <- readLines(fasta)
  
  contig <- str_replace(string = paste(grep(fa,pattern = '^>',invert = FALSE,value = TRUE),collapse = ""),pattern = '>',replacement = '')
  
  fa <- paste(grep(fa,pattern = '^>',invert = TRUE,value = TRUE),collapse = "") %>% str_to_upper()
  
  vcf <- data.frame(CHROM = 'SPIKE',
                    POS = seq(0,nchar(fa)-1,1),
                    ID = '.',
                    REF = unlist(str_extract_all(fa,pattern = "")),
                    ALT = '.',
                    QUAL = '.',
                    FILTER = '.',
                    INFO = '.',
                    stringsAsFactors = FALSE)
  
  vcf.file <- file.path(file.path(bam.folder,'pileup.aseq'),'spike.vcf')
  
  colnames(vcf)[1] <- '#CHROM' 
  
  write.table(vcf,file = vcf.file,sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)
  
  bams <- list.files(bam.folder,pattern = '\\.bam$',full.names = TRUE)
  
  for(bam in bams){
    
    message(bam)
    
    aseq.path <- '/home/casiraghi/aseq/ASEQ'
    
    out <- paste0('out=',file.path(bam.folder,'pileup.aseq'))
    
    cmd <- paste0(aseq.path,' vcf=',vcf.file,' bam=',bam,' threads=1 mbq=30 mrq=30 ',out)
    
    system(cmd)
    
  }
}
