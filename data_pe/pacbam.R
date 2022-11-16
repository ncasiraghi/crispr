library(tidyverse)

# compute pileup

wd <- '/mnt/profile/darosio'

fasta <- file.path(wd,'reference/bwa_reference/spike.fasta')

bfs <- c("/mnt/profile/darosio/data/bowtie2",
         "/mnt/profile/darosio/data/bwa-aln-sampe",
         "/mnt/profile/darosio/data/bwa-mem")

for(bam.folder in bfs){
  
  message(bam.folder)
  
  if(!exists(file.path(bam.folder,'pileup'))){
    dir.create(file.path(bam.folder,'pileup'),showWarnings = FALSE)
  }
  
  vcf <- data.frame(CHROM = 'SPIKE',
                    POS = 0,
                    REF = 'A',
                    ALT = 'G',
                    QUAL = '.',
                    FILTER = '.',
                    INFO = '.',
                    stringsAsFactors = FALSE)
  
  vcf.file <- file.path(file.path(bam.folder,'pileup'),'toy.vcf')
  
  colnames(vcf)[1] <- '#CHROM' 
  
  write.table(vcf,file = vcf.file,sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)
  
  fa <- readLines(fasta)
  
  contig <- str_replace(string = paste(grep(fa,pattern = '^>',invert = FALSE,value = TRUE),collapse = ""),pattern = '>',replacement = '')
  
  fa <- paste(grep(fa,pattern = '^>',invert = TRUE,value = TRUE),collapse = "")
  fa <- str_to_upper(fa)
  
  bed <- data.frame(contig = contig,
                    start = 0,
                    end = nchar(fa) - 1,
                    stringsAsFactors = F)
  
  bed.file <- file.path(bam.folder,'pileup',gsub(basename(fasta),pattern = '\\.fasta$',replacement = '.bed'))
  
  write.table(bed,file = bed.file,sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)
  
  bams <- list.files(bam.folder,pattern = '\\.sorted.bam$',full.names = TRUE)
  
  for(bam in bams){
    
    message(bam)
    
    pacbam.path <- '/home/casiraghi/Documents/unitn/darosio/pacbam/pacbam'
    
    out <- paste0('out=',file.path(bam.folder,'pileup'))
    
    cmd <- paste0(pacbam.path,' bam=',bam,' fasta=',fasta,' vcf=',vcf.file,' bed=',bed.file,' mode=4',' threads=1 ',out)
    
    system(cmd)
    
  }
}


