
wd <- '/BCGLAB/darosio_crispr'

setwd(wd)

bwa_fasta  <- file.path(wd,'reference/bwa_reference/egfp.fasta')
gatk_fasta <- file.path(wd,'reference/gatk_reference/egfp.fasta')

fastq_files <- list.files(file.path(wd,'fastq'),full.names = TRUE,pattern = '_L001_R1_001.fastq.gz$')

# trimming fastqc with trim_galore

trim_galore <- '/BCGLAB/Tools/TrimGalore-0.6.6/trim_galore --path_to_cutadapt ~/.local/bin/cutadapt'

cmd <- paste(trim_galore,'-q 30',paste(fastq_files,collapse = ' '),'-o',unique(dirname(fastq_files)))
system(cmd)

# create index file

cmd <- paste('bwa index',bwa_fasta)
system(cmd)

# Creating the FASTA sequence dictionary file

picard <- 'java -jar /BCGLAB/Tools/picard.jar'
cmd <- paste(picard,'CreateSequenceDictionary',paste0('R=',gatk_fasta))
system(cmd)

cmd <- paste('samtools faidx',gatk_fasta)
system(cmd)

trimmed_fastq_files <- list.files(file.path(wd,'fastq'),full.names = TRUE,pattern = '_trimmed\\.fq\\.gz$')

for(fastq in trimmed_fastq_files){
  
  message(fastq)
  
  sam <- file.path(wd,'sam',paste0(unlist(strsplit(basename(fastq),split = '_'))[1],'.sam'))
  
  id.sample <- gsub(basename(sam),pattern = '\\.sam',replacement = '')
  
  RG <- paste('-R',paste0('"@RG\\tID:',id.sample,'\\tSM:',id.sample,'"'))
  
  cmd <- paste('bwa mem -t 5',RG,bwa_fasta,fastq,'>',sam)
  system(cmd)
  
  bam <- file.path(wd,'bam',paste0(unlist(strsplit(basename(fastq),split = '_'))[1],'.bam'))
  cmd <- paste('samtools view -S -b',sam,'>',bam)
  system(cmd)
  
  sorted.bam <- file.path(wd,'bam',paste0(unlist(strsplit(basename(fastq),split = '_'))[1],'.sorted.bam'))
  cmd <- paste('samtools sort','-o',sorted.bam,bam)
  system(cmd)
  
  cmd <- paste('samtools index',sorted.bam)
  system(cmd)
  
  # GATK
  if(TRUE){
    
    gatk <- 'java -Xmx16G -jar /BCGLAB/Tools/gatk.jar'
    
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
