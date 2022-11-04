library(tidyverse)

wd <- '/Users/ncasiraghi/Documents/unitn/darosio/'

bowtie_fasta  <- "/Users/ncasiraghi/Documents/unitn/darosio/reference/bowtie2_reference/spike.fasta"
gatk_fasta    <- "/Users/ncasiraghi/Documents/unitn/darosio/reference/gatk_reference/spike.fasta"

setwd(wd)

fastq_R1 <- list.files(file.path(wd,'data/fastq'),full.names = TRUE,pattern = '_L001_R1_001.fastq.gz$', recursive = TRUE)
fastq_R2 <- list.files(file.path(wd,'data/fastq'),full.names = TRUE,pattern = '_L001_R2_001.fastq.gz$', recursive = TRUE)

fastq_R1 <- data.frame(sample = sapply(strsplit(basename(fastq_R1),split = '_'),FUN = '[',1),
                       fastq = fastq_R1,
                       stringsAsFactors = FALSE)

fastq_R2 <- data.frame(sample = sapply(strsplit(basename(fastq_R2),split = '_'),FUN = '[',1),
                       fastq = fastq_R2,
                       stringsAsFactors = FALSE)

fastq <- full_join(x = fastq_R1,y = fastq_R2, by = 'sample',suffix = c("_1","_2"))

for(i in seq_len(nrow(fastq))){
  
  message(fastq$sample[i])
  
  message("Step 1	Alignment - Map to Reference")
  
  sam <- file.path(wd,'data/sam_bowtie2',paste0(fastq$sample[i],'.sam'))
  
  id.sample <- fastq$sample[i]
  
  cmd <- paste('bowtie2 -x',str_remove(bowtie_fasta,pattern = '\\.fasta$'),
               '-1',fastq$fastq_1[i],
               '-2',fastq$fastq_2[i],
               '-S',sam,'--rg-id',
               id.sample,'--rg',
               paste0('SM:',fastq$sample[i]),
               '--rg PL:ILLUMINA')
  system(cmd)
  
  message("Step 2	from SAM to sorted BAM")
  
  bam <- file.path(wd,'data/bam_bowtie2',paste0(fastq$sample[i],'.bam'))
  cmd <- paste('samtools view -S -b',sam,'>',bam)
  system(cmd)
  
  sorted.bam <- file.path(wd,'data/bam_bowtie2',paste0(fastq$sample[i],'.sorted.bam'))
  cmd <- paste('samtools sort','-o',sorted.bam,bam)
  system(cmd)
  
  message("Step 3	Mark Duplicates")
  
  mets <- str_replace(string = bam,pattern = '\\.bam$',replacement = ".dedup_metrics.txt")
  
  dedup.bam <- file.path(wd,'data/bam_bowtie2',paste0(fastq$sample[i],'.dedup.sorted.bam'))
  
  cmd <- paste('gatk MarkDuplicates -I',sorted.bam,'-M',mets,'-O',dedup.bam)
  system(cmd)
  
  message("Step 4	index BAM")
  
  cmd <- paste('samtools index',sorted.bam)
  system(cmd)
  
  message("Step 5	HaplotypeCaller")
  
  vcf <- str_replace(string = basename(sorted.bam),pattern = '\\.bam$',replacement = "_raw_variants.vcf")
  
  cmd <- paste('gatk HaplotypeCaller -R',gatk_fasta,'-I',sorted.bam,'-O',file.path(wd,'data/haplotypecaller/bowtie',vcf))
  system(cmd)
  
  message("Step 6	BaseRecalibrator")
  
  recal <- str_replace(sorted.bam,pattern = "\\.bam$",replacement = "_recal_data.table")
  
  cmd <- paste('gatk BaseRecalibrator -R',gatk_fasta,
               '-I',sorted.bam,
               "--known-sites",file.path(wd,'data/haplotypecaller/bowtie',vcf),
               '-O',recal)
  system(cmd)
  
  message("Step 7	ApplyBQSR")  
  
  recal.bam <- str_replace(sorted.bam,pattern = "\\.bam$",replacement = ".recal.bam")
  
  cmd <- paste('gatk ApplyBQSR -R',gatk_fasta,
               '-I',sorted.bam,
               "-bqsr",recal,
               '-O',recal.bam)
  system(cmd)
  
  cmd <- paste('samtools index',recal.bam)
  system(cmd)
  
  message("Step 8	HaplotypeCaller on recal BAM")
  
  vcf <- str_replace(string = basename(recal.bam),pattern = '\\.bam$',replacement = "_raw_variants.vcf")
  
  cmd <- paste('gatk HaplotypeCaller -R',gatk_fasta,'-I',recal.bam,'-O',file.path(wd,'data/haplotypecaller/bowtie/',vcf))
  system(cmd)
  
  message("done")
  
}
