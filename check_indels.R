library(GenomicAlignments)
library(tidyverse)

setwd('/BCGLAB/darosio_crispr/out/')

## This code was designed for BAM files aligned against a specific reference genome "amplicons-deepseq.fasta"
## BAM files should be first sorted by read position and indexed with SAMTools

## Function to get the INDELs events
getIndels <- function(i,reads){
  
  cigar <- reads$cigar[i]
  
  pos <- reads$pos[i]
  
  tab <- as.data.frame(cigarRangesAlongReferenceSpace(cigar, with.ops=TRUE)[[1]])
  
  x = as.data.frame(cigarRangesAlongQuerySpace(cigar, with.ops=TRUE)[[1]])
  tab$tot.width = x$width+tab$width
  
  tab$start <- tab$start + pos
  tab$end <- tab$end + pos

  out <- tab %>% filter(names %in% c('D','I'))
  
  return(out)

}

## Code to extract INDELs statistics from BAM files

bamlist = list.files("/BCGLAB/darosio_crispr/bam",pattern = ".sorted.bam$",full.names = TRUE)

df.totals <- list()
df.reads <- list()
df.indels <- list()

for(file.bam in bamlist){
  
  message(file.bam)
  
  id <- gsub(basename(file.bam),pattern = '.sorted.bam',replacement = '')
  
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=NA),what=c("rname", "pos", "cigar","mapq"))
  bam <- scanBam(file.bam, param=param)[[1]]
  
  tot.reads <- as.data.frame(bam) %>% 
    filter(mapq >= 30) %>% 
    mutate(sample = id)
  
  df.totals[[id]] <- tot.reads
  
  reads <- tot.reads %>% 
    filter(grepl(cigar,pattern = "D|I" )) %>% 
    filter(!grepl(cigar,pattern = "S|H")) %>% 
    mutate(sample = id)
  
  df.reads[[id]] <- reads
  
  indels <- mclapply(seq_len(nrow(reads)),FUN = getIndels,reads,mc.cores = 30)
  
  indels <- do.call(rbind,indels) %>% 
    mutate(sample = id)
  
  df.indels[[id]] <- indels
  
}

df.totals <- do.call(rbind,df.totals)
df.reads <- do.call(rbind,df.reads)
df.indels <- do.call(rbind,df.indels)

save(df.totals,df.indels,df.reads,file = file.path('/BCGLAB/darosio_crispr/out/data_indels.RData'),compress = TRUE)

# load('/BCGLAB/darosio_crispr/out/data_indels.RData')

all <- df.totals %>%
  group_by(sample) %>%
  summarise(n = n())

p <- ggplot(data=all, aes(x=sample, y=n)) +
  geom_bar(stat="identity") + ylab('Number of reads') + ggtitle('Selected reads with mapq > 30')

ggsave(filename = 'pdf/all_reads.pdf', plot = p, width = 150,height = 150,dpi = 300,units = 'mm',device = 'pdf')

cnt <- df.reads %>%
  group_by(sample) %>%
  summarise(n = n())

p <- ggplot(data=cnt, aes(x=sample, y=n)) +
  geom_bar(stat="identity") + ylab('Number of reads') + ggtitle('Selected reads with D or I in CIGAR [ mapq > 30 ]')

ggsave(filename = 'pdf/DI_reads_cnt.pdf', plot = p, width = 150,height = 150,dpi = 300,units = 'mm',device = 'pdf')

cnt$n <- cnt$n / all$n

p <- ggplot(data=cnt, aes(x=sample, y=n)) +
  geom_bar(stat="identity") + ylab('Fraction of reads') + ggtitle('Selected reads with D or I in CIGAR [ mapq > 30 ]')

ggsave(filename = 'pdf/DI_reads_frc.pdf', plot = p, width = 150,height = 150,dpi = 300,units = 'mm',device = 'pdf')

# distribution of indels length

m <- df.indels %>%
  group_by(sample,names) %>% 
  summarise(n = n(), mean_width = mean(tot.width), median_width = median(tot.width))

p <- ggplot(m, aes(x=sample, y=n)) +
  geom_bar(stat="identity") + ylab('number of events') + 
  facet_wrap(~names)

ggsave(filename = 'pdf/DI_count_events.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- ggplot(df.indels, aes(x=sample, y=tot.width)) + 
  geom_boxplot(varwidth = TRUE) +
  facet_wrap(~names) + ylab('width')

ggsave(filename = 'pdf/DI_widths.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

# 

fasta <- readLines('/BCGLAB/darosio_crispr/reference/egfp.fasta')
fasta <- paste(grep(fasta,pattern = '^>',invert = TRUE,value = TRUE),collapse = "")
fasta <- str_to_upper(fasta)

nchar(fasta)

df.indels <- df.indels %>% mutate(start = start + 1)

p <- ggplot(df.indels, aes(x=start)) + 
  geom_histogram(binwidth=1) +
  facet_wrap(~sample,nrow = 5) + xlim(100,300) +
  geom_vline(xintercept = 154,color = "#b2182b", size=0.4) +
  geom_vline(xintercept = 202,color = "#006837", size=0.4) + ggtitle('Starting position of D,I events')

ggsave(filename = 'pdf/DI_starts.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- ggplot(df.indels %>% filter(tot.width > 1), aes(x=start)) + 
  geom_histogram(binwidth=1) +
  facet_wrap(~sample,nrow = 5) + xlim(100,300) +
  geom_vline(xintercept = 154,color = "#b2182b", size=0.4) +
  geom_vline(xintercept = 202,color = "#006837", size=0.4) + ggtitle('Starting position of D,I events with width > 1')

ggsave(filename = 'pdf/DI_starts_no_single.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

# heatmap

mat <- matrix(data = 0,nrow = 5,ncol = nchar(fasta))
row.names(mat) <- unique(df.indels$sample)

mat_nos <- matrix(data = 0,nrow = 5,ncol = nchar(fasta))
row.names(mat_nos) <- unique(df.indels$sample)

for(id in unique(df.indels$sample)){
  message(id)
  
  h <- df.indels %>% 
    filter(sample == id)
  
  for(idx in seq_len(nrow(h))){
    
    x <- sort(h$start[idx]:h$end[idx])
    
    mat[id,x] <- mat[id,x] + 1 
    
  }
  
  h <- df.indels %>% 
    filter(sample == id) %>% 
    filter(tot.width > 1)
  
  for(idx in seq_len(nrow(h))){
    
    x <- sort(h$start[idx]:h$end[idx])
    
    mat_nos[id,x] <- mat_nos[id,x] + 1 
    
  }
  
}

mydata <- c()

for(id in unique(rownames(mat))){

  df <- data.frame(sample = id,
                   value = as.numeric(mat[id,]),
                   pos = seq_len(ncol(mat)),
                   stringsAsFactors = FALSE)
  
  mydata <- rbind(mydata,df)
  
}

p <- ggplot(data=mydata, aes(x=pos, y=value)) +
  geom_bar(stat="identity",fill='grey60',color='grey60',width = 1) +
  facet_wrap(~sample,nrow = 5) + xlim(100,300) +
  geom_vline(xintercept = 154,color = "#b2182b", size=0.4) +
  geom_vline(xintercept = 202,color = "#006837", size=0.4) + ggtitle('Positions covered by a D,I events')

ggsave(filename = 'pdf/DI_pos_covered.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')


# only DI with width > 1

mydata <- c()

for(id in unique(rownames(mat_nos))){
  
  df <- data.frame(sample = id,
                   value = as.numeric(mat_nos[id,]),
                   pos = seq_len(ncol(mat_nos)),
                   stringsAsFactors = FALSE)
  
  mydata <- rbind(mydata,df)
  
}

p <- ggplot(data=mydata, aes(x=pos, y=value)) +
  geom_bar(stat="identity",fill='grey60',color='grey60',width = 1) +
  facet_wrap(~sample,nrow = 5) + xlim(100,300) +
  geom_vline(xintercept = 154,color = "#b2182b", size=0.4) +
  geom_vline(xintercept = 202,color = "#006837", size=0.4) + ggtitle('Positions covered by a D,I events with width > 1')

ggsave(filename = 'pdf/DI_pos_covered_no_singles.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')









