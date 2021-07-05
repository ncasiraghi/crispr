
wd <- '/BCGLAB/darosio_crispr'

fasta <- file.path(wd,'reference/pires-puro3-0-egfp-nostro.fasta')

fastq_files <- list.files(file.path(wd,'fastq'),full.names = TRUE,pattern = '_R1_001.fastq.gz$')

cmd <- paste('bwa index',fasta) 
system(cmd)

for(fastq in fastq_files){
  sam <- file.path(wd,'sam',paste0(unlist(strsplit(basename(fastq),split = '_'))[1],'.sam'))
  
  cmd <- paste('bwa mem',fasta,fastq,'>',sam)
  system(cmd)
  
  bam <- file.path(wd,'bam',paste0(unlist(strsplit(basename(fastq),split = '_'))[1],'.bam'))
  cmd <- paste('samtools view -S -b',sam,'>',bam)
  system(cmd)
  
  sorted.bam <- file.path(wd,'bam',paste0(unlist(strsplit(basename(fastq),split = '_'))[1],'.sorted.bam'))
  cmd <- paste('samtools sort',bam,'>',sorted.bam)
  system(cmd)

  cmd <- paste('samtools index',sorted.bam)
  system(cmd)
}
