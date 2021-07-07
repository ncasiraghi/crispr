
wd <- '/BCGLAB/darosio_crispr'

setwd(wd)

fasta <- file.path(wd,'reference/egfp.fasta')

fastq_files <- list.files(file.path(wd,'fastq'),full.names = TRUE,pattern = '_L001_R1_001.fastq.gz$')

if(TRUE){
  cmd <- paste('bwa index',fasta) 
  system(cmd)
}

for(fastq in fastq_files){
  sam <- file.path(wd,'sam',paste0(unlist(strsplit(basename(fastq),split = '_'))[1],'.sam'))
  
  cmd <- paste('bwa mem',fasta,fastq,'>',sam)
  system(cmd)
  
  bam <- file.path(wd,'bam',paste0(unlist(strsplit(basename(fastq),split = '_'))[1],'.bam'))
  cmd <- paste('samtools view -S -b',sam,'>',bam)
  system(cmd)
  
  sorted.bam <- file.path(wd,'bam',paste0(unlist(strsplit(basename(fastq),split = '_'))[1],'.sorted.bam'))
  cmd <- paste('samtools sort','-o',sorted.bam,bam)
  system(cmd)

  cmd <- paste('samtools index',sorted.bam)
  system(cmd)
}
