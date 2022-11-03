library(tidyverse)

wd <- '/Users/ncasiraghi/Documents/unitn/darosio'

# bwa index spike.fasta
bwa_fasta <- "/Users/ncasiraghi/Documents/unitn/darosio/reference/bwa_reference/spike.fasta"

setwd(wd)

fastq_R1 <- list.files(file.path(wd,'data/fastq'),full.names = TRUE,pattern = '_L001_R1_001.fastq.gz$', recursive = TRUE)
fastq_R2 <- list.files(file.path(wd,'data/fastq'),full.names = TRUE,pattern = '_L001_R2_001.fastq.gz$', recursive = TRUE)

fastq_R1 <- data.frame(sample = sapply(strsplit(basename(fastq_R1),split = '_'),FUN = '[',1),
                       fastq = fastq_R1,
                       read = 1,
                       stringsAsFactors = FALSE)

fastq_R2 <- data.frame(sample = sapply(strsplit(basename(fastq_R2),split = '_'),FUN = '[',1),
                       fastq = fastq_R2,
                       read = 2,
                       stringsAsFactors = FALSE)

fastq <- rbind(fastq_R1,fastq_R2) %>% arrange(sample,read)

# trimgalore trimming

if(FALSE){
  
  message("trimming fastq files")
  
  trim_galore <- 'trim_galore --path_to_cutadapt /Users/ncasiraghi/opt/miniconda3/envs/crispr/bin/cutadapt'
  
  cmd <- paste(trim_galore,'-q 30 --paired --retain_unpaired',paste(fastq$fastq,collapse = ' '),'-o',file.path(wd,'data/fastq_trimmed'))
  
  system(cmd)
  
}

trim_R1 <- list.files(file.path(wd,'data/fastq_trimmed'),full.names = TRUE,pattern = '\\_val_1.fq\\.gz$')
trim_R2 <- list.files(file.path(wd,'data/fastq_trimmed'),full.names = TRUE,pattern = '\\_val_2.fq\\.gz$')

trim_R1 <- data.frame(sample = sapply(strsplit(basename(trim_R1),split = '_'),FUN = '[',1),
                      fastq = trim_R1,
                      stringsAsFactors = FALSE)

trim_R2 <- data.frame(sample = sapply(strsplit(basename(trim_R2),split = '_'),FUN = '[',1),
                      fastq = trim_R2,
                      stringsAsFactors = FALSE)

trim <- full_join(x = trim_R1,y = trim_R2, by = 'sample',suffix = c("_1","_2"))

for(i in seq_len(nrow(trim))){
  
  message(trim$sample[i])
  
  sam <- file.path(wd,'data/sam_bwa',paste0(trim$sample[i],'.sam'))
  
  id.sample <- trim$sample[i]
  
  RG <- paste('-R',paste0('"@RG\\tID:',id.sample,'\\tSM:',id.sample,'"'))
  
  cmd <- paste('bwa mem -t 1 -Y',RG,bwa_fasta,trim$fastq_1[i],trim$fastq_2[i],'>',sam)
  system(cmd)
  
  bam <- file.path(wd,'data/bam_bwa',paste0(trim$sample[i],'.bam'))
  
  cmd <- paste('samtools view -S -b',sam,'>',bam)
  system(cmd)
  
  sorted.bam <- file.path(wd,'data/bam_bwa',paste0(trim$sample[i],'.sorted.bam'))
  cmd <- paste('samtools sort','-o',sorted.bam,bam)
  system(cmd)
  
  cmd <- paste('samtools index',sorted.bam)
  system(cmd)
  
}  

