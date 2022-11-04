if(build.reference){

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

# bwa index spike.fasta
