
wd <- "/mnt/profile/darosio"

faidx_fasta <- file.path(wd,'reference/samtools_reference/spike.fasta')

for(bam.folder in c("bowtie2","bwa-aln-sampe","bwa-mem")){
  
  message(bam.folder)

  if(!exists(file.path(bam.folder,'mpileup'))){
    dir.create(file.path(bam.folder,'mpileup'),showWarnings = FALSE)
  }
  
  bams <- list.files(bam.folder,pattern = '\\.sorted.bam$',full.names = TRUE)
  
  for(bam in bams){
    
    message(bam)
    
    mpileup <- file.path(bam.folder,"mpileup",gsub(basename(bam),pattern = '\\.bam$',replacement = ".mpileup"))
    
    cmd <- paste("bcftools mpileup","-f",faidx_fasta,"-q",30,'-Q',30,"-Ov -o",mpileup,bam)
    
    system(cmd)
    
  }
  
}