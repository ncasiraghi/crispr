library(GenomicAlignments)
library(tidyverse)

# data <- 'bowtie'
data <- 'bwa'

if(data == 'bowtie'){
  wd <- '/BCGLAB/darosio_crispr/bowtie2/out'
  bam_folder <- file.path(wd,'bam')
}

if(data == 'bwa'){
  wd <- '/BCGLAB/darosio_crispr/bwa/out'
  bam_folder <- file.path(wd,'bam')
}

rdatafolder <- file.path(wd,'rdata')

if(!file.exists(file.path(rdatafolder))){
  dir.create(file.path(rdatafolder),showWarnings = FALSE)
} else{
  unlink(file.path(rdatafolder),recursive = TRUE)
  dir.create(file.path(rdatafolder), showWarnings = FALSE)
}

setwd(rdatafolder)

# generate data

## Function to get the INDELs events
getIndels <- function(i,di.reads){
  
  cigar <- di.reads$cigar[i]
  
  pos <- di.reads$pos[i]
  
  tab <- as.data.frame(cigarRangesAlongReferenceSpace(cigar, with.ops=TRUE)[[1]])
  
  x = as.data.frame(cigarRangesAlongQuerySpace(cigar, with.ops=TRUE)[[1]])
  tab$tot.width = x$width+tab$width
  
  tab$start <- tab$start + pos
  tab$end <- tab$end + pos
  
  out <- tab %>%
    filter(names %in% c('D','I')) %>%
    mutate(qname = di.reads$qname[i])
  
  return(out)
  
}

## Code to extract INDELs statistics from BAM files

pattern <- "\\.sorted\\.realigned\\.bam$"

bamlist <- list.files(bam_folder,pattern = pattern,full.names = TRUE) 

bamlist <- grep(bamlist,pattern = 'Undetermined',value = TRUE,invert = TRUE)

for(file.bam in bamlist){
  
  message(file.bam)
  
  id <- gsub(basename(file.bam),pattern = pattern,replacement = '')
  
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=NA),what=c("qname","pos", "cigar","mapq"))
  bam <- scanBam(file.bam, param=param)[[1]]
  
  tot.reads <- as.data.frame(bam) %>% 
    filter(mapq >= 30) %>% 
    mutate(sample = id) %>% 
    separate(qname,sep = ':',into = letters[1:7],remove = TRUE) %>% 
    select(e:sample) %>% 
    unite('qname',e:g,sep = ":",remove = TRUE)
  
  rm(bam)
  
  save(tot.reads,file = paste0('tot.reads_',id,'.RData'),compress = TRUE)
  
  di.reads <- tot.reads %>% 
    filter(grepl(cigar,pattern = "D|I" )) %>% 
    filter(!grepl(cigar,pattern = "S|H"))
  
  rm(tot.reads)
  
  save(di.reads,file = paste0('di.reads_',id,'.RData'),compress = TRUE)
  
  indels <- mclapply(seq_len(nrow(di.reads)),FUN = getIndels,di.reads,mc.cores = 30)
  
  rm(di.reads)
  
  indels <- do.call(rbind,indels) %>% mutate(sample = id)
  
  save(indels,file = paste0('indels_',id,'.RData'),compress = TRUE)
  
  rm(indels)
  
}

df.totals <- list()
for( rdata in list.files(rdatafolder,pattern = 'tot.reads',full.names = TRUE)){
  message(rdata)
  load(rdata)
  df.totals[[basename(rdata)]] <- tot.reads
}
df.totals <- do.call(rbind, df.totals)

df.reads <- list()
for( rdata in list.files(rdatafolder,pattern = 'di.reads',full.names = TRUE)){
  message(rdata)
  load(rdata)
  df.reads[[basename(rdata)]] <- di.reads
}
df.reads <- do.call(rbind, df.reads)

df.indels <- list()
for( rdata in list.files(rdatafolder,pattern = 'indels',full.names = TRUE)){
  message(rdata)
  load(rdata)
  df.indels[[basename(rdata)]] <- indels
}
df.indels <- do.call(rbind, df.indels)

save(df.totals,df.reads,df.indels,file = 'fulldata.RData',compress = TRUE)
