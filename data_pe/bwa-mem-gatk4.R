library(tidyverse)

# https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/

wd <- '/Users/ncasiraghi/Documents/unitn/darosio'

# bwa index spike.fasta
bwa_fasta  <- "/Users/ncasiraghi/Documents/unitn/darosio/reference/bwa_reference/spike.fasta"
gatk_fasta <- "/Users/ncasiraghi/Documents/unitn/darosio/reference/gatk_reference/spike.fasta"

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

for(i in seq_len(nrow(fastq))[1]){
  
  message(fastq$sample[i])
  
  message("Step 1	Alignment - Map to Reference")
  
  sam <- file.path(wd,'data/sam_bwa',paste0(fastq$sample[i],'.sam'))
  
  id.sample <- fastq$sample[i]
  
  RG <- paste('-R',paste0('"@RG\\tID:',id.sample,'\\tSM:',id.sample,'"'))
  
  cmd <- paste('bwa mem -Y',RG,bwa_fasta,fastq$fastq_1[i],fastq$fastq_2[i],'>',sam)
  system(cmd)
  
  message("Step 2	Mark Duplicates + Sort")
  
  mets <- str_replace(string = sam,pattern = '\\.sam$',replacement = ".dedup_metrics.txt")
  
  bam <- file.path(wd,'data/bam_bwa',paste0(fastq$sample[i],'_sorted_dedup_reads.bam'))
  
  cmd <- paste('gatk MarkDuplicatesSpark -I',sam,'-M',mets,'-O',bam)
  
  message("Step 3	Call Variants")
  
  vcf <- str_replace(string = basename(bam),pattern = '\\.bam$',replacement = "_raw_variants.vcf")
  
  cmd <- paste('gatk HaplotypeCaller -R',gatk_fasta,'-I',bam,'-o',file.path(wd,'data/haplotypecaller',vcf))
  
  message("done")
  
}  

