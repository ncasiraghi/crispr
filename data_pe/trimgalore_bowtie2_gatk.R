library(tidyverse)

wd <- '/Users/ncasiraghi/Documents/unitn/darosio'

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

# bowtie2 reference

setwd(file.path(wd,'reference/bowtie2_reference'))

bowtie_fasta <- file.path(wd,'reference/bowtie2_reference/spike.fasta')

cmd <- paste('bowtie2-build -f',bowtie_fasta,gsub(basename(bowtie_fasta),pattern = '\\.fasta',replacement = ''))
system(cmd)

# gatk reference

gatk_fasta <- file.path(wd,'reference/gatk_reference/spike.fasta')

cmd <- paste('picard CreateSequenceDictionary',paste('-R',gatk_fasta))
system(cmd)

cmd <- paste('samtools faidx',gatk_fasta)
system(cmd)

# trimgalore trimming

setwd(wd)

trim_galore <- 'trim_galore --path_to_cutadapt /Users/ncasiraghi/opt/miniconda3/envs/darosio/bin/cutadapt'

cmd <- paste(trim_galore,'-q 30 --paired --retain_unpaired',paste(fastq$fastq,collapse = ' '),'-o',file.path(wd,'data/fastq_trimmed'))
system(cmd)
  
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
  
  sam <- file.path(wd,'data/sam_bowtie2',paste0(trim$sample[i],'.sam'))
  
  id.sample <- trim$sample[i]
  
  cmd <- paste('bowtie2 -x',str_remove(bowtie_fasta,pattern = '\\.fasta$'),'-1',trim$fastq_1[i],'-2',trim$fastq_2[i],'-S',sam,'--rg-id',id.sample,'--rg',paste0('SM:',trim$sample[i]),'--rg PL:ILLUMINA')
  system(cmd)
  
  bam <- file.path(wd,'data/bam_bowtie2',paste0(trim$sample[i],'.bam'))
  cmd <- paste('samtools view -S -b',sam,'>',bam)
  system(cmd)
  
  sorted.bam <- file.path(wd,'data/bam_bowtie2',paste0(trim$sample[i],'.sorted.bam'))
  cmd <- paste('samtools sort','-o',sorted.bam,bam)
  system(cmd)
  
  cmd <- paste('samtools index',sorted.bam)
  system(cmd)
  
  # GATK
  if(TRUE){
    
    gatk <- 'gatk'
    
    intervals <- gsub(sorted.bam,pattern = '\\.sorted\\.bam$',replacement = '.bam.intervals')
    cmd <- paste(gatk,
                 '-T RealignerTargetCreator',
                 '-S SILENT',
                 '-I',sorted.bam,
                 '-R',gatk_fasta,
                 '-nt 10',
                 '-o',intervals)
    system(cmd)
    
    realigned.bam <- gsub(sorted.bam,pattern = '\\.bam$',replacement = '.realigned.bam')
    cmd <- paste(gatk,
                 '-S SILENT',
                 '-I',sorted.bam,
                 '-R',gatk_fasta,
                 '-T IndelRealigner',
                 '-maxReads 5000000',
                 '-targetIntervals',intervals,
                 '-o',realigned.bam)
    system(cmd)
    
    cmd <- paste('samtools index',realigned.bam)
    system(cmd)
    
  }
  
}
  

