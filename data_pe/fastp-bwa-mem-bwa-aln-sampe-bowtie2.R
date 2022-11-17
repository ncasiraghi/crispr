
wd <- '/mnt/profile/darosio'

# bwa index spike.fasta
bwa_fasta    <- file.path(wd,"reference/bwa_reference/spike.fasta")
gatk_fasta   <- file.path(wd,"reference/gatk_reference/spike.fasta")
bowtie_fasta <- file.path(wd,"reference/bowtie2_reference/spike.fasta")

setwd(wd)

# create out folders

if(!exists(file.path(wd,'data/fastp'))){
  dir.create(file.path(wd,'data/fastp'),showWarnings = TRUE)
}

if(!exists(file.path(wd,'data/bwa-mem'))){
  dir.create(file.path(wd,'data/bwa-mem'),showWarnings = TRUE)
}

if(!exists(file.path(wd,'data/bwa-aln-sampe'))){
  dir.create(file.path(wd,'data/bwa-aln-sampe'),showWarnings = TRUE)
}

if(!exists(file.path(wd,'data/bowtie2'))){
  dir.create(file.path(wd,'data/bowtie2'),showWarnings = TRUE)
}

# trimming fastq with fastp

if(FALSE){

# system("conda activate fastp")

# fastq
  
fastq_R1 <- list.files(file.path(wd,'data/fastq'),full.names = TRUE,pattern = '_R1_001.fastq.gz$',recursive = TRUE)
fastq_R2 <- list.files(file.path(wd,'data/fastq'),full.names = TRUE,pattern = '_R2_001.fastq.gz$',recursive = TRUE)

fastq_R1 <- data.frame(sample = sapply(strsplit(basename(fastq_R1),split = '_'),FUN = '[',1),
                      fastq = fastq_R1,
                      stringsAsFactors = FALSE)

fastq_R2 <- data.frame(sample = sapply(strsplit(basename(fastq_R2),split = '_'),FUN = '[',1),
                      fastq = fastq_R2,
                      stringsAsFactors = FALSE)

fastq <- merge(x = fastq_R1,y = fastq_R2, by = 'sample',suffix = c("_1","_2"),all = FALSE)

for(i in seq_len(nrow(fastq))){
  
  id.sample <- fastq$sample[i]
  
  message(id.sample)
  
  # trimming
  
  out.R1.fq.gz <- file.path(wd,'data/fastp',paste0(id.sample,'_R1_001.fastq.gz'))
  out.R2.fq.gz <- file.path(wd,'data/fastp',paste0(id.sample,'_R2_001.fastq.gz'))
  
  cmd <- paste("fastp -i",fastq$fastq_1[i],"-I",fastq$fastq_2[i],"-o",out.R1.fq.gz,"-O",out.R2.fq.gz)
  
  system(cmd)
  
  # trimming + merge PE reads
  
  
  
}

}

if(FALSE){

# use trimmed fastq

fastq_R1 <- list.files(file.path(wd,'data/fastp'),full.names = TRUE,pattern = '_R1_001.fastq.gz$',recursive = TRUE)
fastq_R2 <- list.files(file.path(wd,'data/fastp'),full.names = TRUE,pattern = '_R2_001.fastq.gz$',recursive = TRUE)

fastq_R1 <- data.frame(sample = sapply(strsplit(basename(fastq_R1),split = '_'),FUN = '[',1),
                       fastq = fastq_R1,
                       stringsAsFactors = FALSE)

fastq_R2 <- data.frame(sample = sapply(strsplit(basename(fastq_R2),split = '_'),FUN = '[',1),
                       fastq = fastq_R2,
                       stringsAsFactors = FALSE)

fastq <- merge(x = fastq_R1,y = fastq_R2, by = 'sample',suffix = c("_1","_2"),all = FALSE)

fastq

# perform alignments
# system("conda activate crispr")

for(i in seq_len(nrow(fastq))){
  
  id.sample <- fastq$sample[i]
  
  message(id.sample)
  
  # bwa mem
  out.folder <- file.path(wd,'data/bwa-mem')
  
  sam <- file.path(out.folder,paste0(id.sample,'.sam'))

  RG <- paste('-R',paste0('"@RG\\tID:',id.sample,'\\tPL:ILLUMINA','\\tSM:',id.sample,'"'))
  
  cmd <- paste('bwa mem -Y',RG,bwa_fasta,fastq$fastq_1[i],fastq$fastq_2[i],'>',sam)
  system(cmd)

  bam <- file.path(out.folder,paste0(id.sample,'.bam'))
  
  cmd <- paste('samtools view -S -b',sam,'>',bam)
  system(cmd)
  
  sorted.bam <- file.path(out.folder,paste0(id.sample,'.sorted.bam'))
  cmd <- paste('samtools sort','-o',sorted.bam,bam)
  system(cmd)
  
  cmd <- paste('samtools index',sorted.bam)
  system(cmd)
  
  system(paste("rm",sam,bam))
  
  # bwa aln + sampe
  out.folder <- file.path(wd,'data/bwa-aln-sampe')
  
  aln_sa1.sai <- file.path(out.folder,paste0(id.sample,'_aln_sa1.sai'))
  aln_sa2.sai <- file.path(out.folder,paste0(id.sample,'_aln_sa2.sai'))
  
  cmd <- paste("bwa aln",bwa_fasta,fastq$fastq_1[i],">",aln_sa1.sai)
  system(cmd)
  
  cmd <- paste("bwa aln",bwa_fasta,fastq$fastq_2[i],">",aln_sa2.sai)
  system(cmd)
  
  sam <- file.path(out.folder,paste0(id.sample,'.sam'))
  
  cmd <- paste("bwa sampe",bwa_fasta,aln_sa1.sai,aln_sa2.sai,fastq$fastq_1[i],fastq$fastq_2[i],">",sam)
  system(cmd)
  
  bam <- file.path(out.folder,paste0(id.sample,'.bam'))
  
  cmd <- paste('samtools view -S -b',sam,'>',bam)
  system(cmd)
  
  sorted.bam <- file.path(out.folder,paste0(id.sample,'.sorted.bam'))
  cmd <- paste('samtools sort','-o',sorted.bam,bam)
  system(cmd)
  
  cmd <- paste('samtools index',sorted.bam)
  system(cmd)
  
  system(paste("rm",sam,bam))
  
  # bowtie2
  out.folder <- file.path(wd,'data/bowtie2')
  
  sam <- file.path(out.folder,paste0(id.sample,'.sam'))
  
  cmd <- paste('bowtie2 -x',gsub(bowtie_fasta,pattern = '\\.fasta$',replacement = ""),
               '-1',fastq$fastq_1[i],
               '-2',fastq$fastq_2[i],
               '-S',sam,'--rg-id',
               id.sample,'--rg',
               paste0('SM:',id.sample),
               '--rg PL:ILLUMINA')
  system(cmd)
  
  bam <- file.path(out.folder,paste0(id.sample,'.bam'))
  
  cmd <- paste('samtools view -S -b',sam,'>',bam)
  system(cmd)
  
  sorted.bam <- file.path(out.folder,paste0(id.sample,'.sorted.bam'))
  cmd <- paste('samtools sort','-o',sorted.bam,bam)
  system(cmd)
  
  cmd <- paste('samtools index',sorted.bam)
  system(cmd)
  
  system(paste("rm",sam,bam))
  
}  

}
message("[done.]")