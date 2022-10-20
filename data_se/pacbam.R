library(GenomicAlignments)
library(tidyverse)

if(FALSE){
  
  # load biased reads
  load('/BCGLAB/darosio_crispr/out/rdata/reads_to_exclude.RData')
  
  # FilterSamReads
  wd <- '/BCGLAB/darosio_crispr/bam'
  
  setwd(wd)
  
  bamlist = list.files("/BCGLAB/darosio_crispr/bam",pattern = ".sorted.bam$",full.names = TRUE)
  
  for(file.bam in bamlist){
    
    message(file.bam)
    
    sn <- gsub(basename(file.bam),pattern = '.sorted.bam',replacement = '')
    
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=NA),what=c("qname","pos", "cigar","mapq"))
    bam <- scanBam(file.bam, param=param)[[1]]
    
    rte <- as.data.frame(bam) %>% 
      filter(mapq >= 30) %>% 
      mutate(sample = sn) %>% 
      separate(qname,sep = ':',into = letters[1:7],remove = FALSE) %>% 
      unite('id',e,f,g,sample,sep = ':',remove = TRUE) %>% 
      filter(id %in% bias$id) %>% 
      pull(qname)
    
    fileConn <- file('read_names.txt')
    writeLines(as.character(rte),fileConn)
    close(fileConn)
    
    picard <- 'java -jar /BCGLAB/Tools/picard.jar FilterSamReads'
    
    cmd <- paste0(picard,
                  ' I=',file.bam,
                  ' O=',gsub(file.bam,pattern = '\\.bam',replacement = '.filtered.bam'),
                  ' READ_LIST_FILE=read_names.txt',
                  ' FILTER=excludeReadList')
    system(cmd)
    
    cmd <- paste('samtools index',gsub(file.bam,pattern = '\\.bam',replacement = '.filtered.bam'))
    system(cmd)
    
  }
  
}

# compute pileup
# wd <- '/BCGLAB/darosio_crispr/bowtie2'
wd <- '/BCGLAB/darosio_crispr/bwa'

if(!file.exists(file.path(wd,'pacbam'))){
  dir.create(file.path(wd,'pacbam'))
}

setwd(file.path(wd,'pacbam'))

# create toy vcf

vcf <- data.frame(CHROM = 'egfp',
                  POS = 0,
                  REF = 'A',
                  ALT = 'G',
                  QUAL = '.',
                  FILTER = '.',
                  INFO = '.',
                  stringsAsFactors = FALSE)

vcf.file <- file.path(wd,'pacbam','toy.vcf')

colnames(vcf)[1] <- '#CHROM' 

write.table(vcf,file = vcf.file,sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)

# create bed

fasta <- '/BCGLAB/darosio_crispr/reference/bwa_reference/egfp.fasta'

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

bam_files <- list.files(file.path(wd,'bam'),full.names = TRUE,pattern = '\\.sorted\\.realigned\\.bam$')

for(bam in bam_files){
  
  message(bam)
  
  cmd <- paste0('/CIBIO/sharedRL/Projects/PaCBAM/git_repo/pacbam/pacbam',' bam=',bam,' fasta=',fasta,' vcf=',vcf.file,' bed=',file.path(wd,'pacbam',bed.file),' mode=4',' threads=20')
  
  system(cmd)
  
}
