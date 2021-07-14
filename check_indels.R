library(GenomicAlignments)
library(tidyverse)

wd <- '/BCGLAB/darosio_crispr/out/'

setwd(wd)

nick <- data.frame(sample = c('D10A','L-16','L-JON','L-PUC'), site = c(202,202,154,154),stringsAsFactors = F)

## This code was designed for BAM files aligned against a specific reference genome "amplicons-deepseq.fasta"
## BAM files should be first sorted by read position and indexed with SAMTools

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

bamlist = list.files("/BCGLAB/darosio_crispr/bam",pattern = ".sorted.bam$",full.names = TRUE)

for(file.bam in bamlist){
  
  message(file.bam)
  
  id <- gsub(basename(file.bam),pattern = '.sorted.bam',replacement = '')
  
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=NA),what=c("qname","pos", "cigar","mapq"))
  bam <- scanBam(file.bam, param=param)[[1]]
  
  tot.reads <- as.data.frame(bam) %>% 
    filter(mapq >= 30) %>% 
    mutate(sample = id) %>% 
    separate(qname,sep = ':',into = letters[1:7],remove = TRUE) %>% 
    select(e:sample) %>% 
    unite('qname',e:g,sep = ":",remove = TRUE)
  
  rm(bam)
  
  save(tot.reads,file = file.path(wd,'rdata',paste0('tot.totals_',id,'.RData')),compress = TRUE)
  
  di.reads <- tot.reads %>% 
    filter(grepl(cigar,pattern = "D|I" )) %>% 
    filter(!grepl(cigar,pattern = "S|H"))
  
  rm(tot.reads)
  
  save(di.reads,file = file.path(wd,'rdata',paste0('di.reads_',id,'.RData')),compress = TRUE)
  
  indels <- mclapply(seq_len(nrow(di.reads)),FUN = getIndels,di.reads,mc.cores = 10)
  
  rm(di.reads)
  
  indels <- do.call(rbind,indels) %>% mutate(sample = id)
  
  save(indels,file = file.path(wd,'rdata',paste0('indels_',id,'.RData')),compress = TRUE)
  
  rm(indels)
  
}

# df.totals <- do.call(rbind,df.totals)
# df.reads <- do.call(rbind,df.reads)
# df.indels <- do.call(rbind,df.indels)

# save(df.totals,df.indels,df.reads,file = file.path('/BCGLAB/darosio_crispr/out/data_indels.RData'),compress = TRUE)

# load('/BCGLAB/darosio_crispr/out/data_indels.RData')
if(FALSE){
  

all <- df.totals %>%
  group_by(sample) %>%
  summarise(n = n())

p <- ggplot(data=all, aes(x=sample, y=n)) +
  geom_bar(stat="identity") + ylab('Number of reads') + ggtitle('Selected reads with mapq > 30')

ggsave(filename = 'pdf/all_reads.pdf', plot = p, width = 150,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- df.reads %>%
  group_by(sample) %>%
  summarise(n = n()) %>% 
  ggplot(., aes(x=sample, y=n)) +
  geom_bar(stat="identity") + ylab('Number of reads') + ggtitle('Selected reads with D or I in CIGAR [ mapq > 30 ]')

ggsave(filename = 'pdf/DI_reads_cnt.pdf', plot = p, width = 150,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- df.reads %>%
  group_by(sample) %>%
  summarise(n = n()) %>% 
  left_join(.,all,by='sample',suffix = c('_reads','_all')) %>% 
  mutate(fraction = n_reads/n_all) %>% 
  ggplot(., aes(x=sample, y=fraction)) +
  geom_bar(stat="identity") + ylab('Fraction of reads') + ggtitle('Selected reads with D or I in CIGAR [ mapq > 30 ]')

ggsave(filename = 'pdf/DI_reads_frc.pdf', plot = p, width = 150,height = 150,dpi = 300,units = 'mm',device = 'pdf')

# distribution frequency

tot <- df.indels %>% 
  group_by(sample) %>% 
  summarise(total = n())

p <- df.indels %>% 
  mutate(id = paste(start,end,names,tot.width,sep = ':')) %>% 
  group_by(id,sample) %>% 
  summarise(n = n()) %>% 
  left_join(.,tot,by = 'sample') %>% 
  mutate(freq = n/total) %>% 
  ggplot(.,aes(x=sample,y=freq)) +
  geom_boxplot(varwidth = TRUE) + ggtitle('D,I events frequency (an event is defined as start:end:D/I:width)')

ggsave(filename = 'pdf/DI_events_frequency.pdf', plot = p, width = 150,height = 150,dpi = 300,units = 'mm',device = 'pdf')

# count events

p <- df.indels %>% 
  group_by(sample,names) %>% 
  summarise(n = n()) %>% 
  ggplot(., aes(x=sample, y=n)) +
  geom_bar(stat="identity",position="dodge") + ylab('Number of events') + 
  facet_wrap(~names) 

ggsave(filename = 'pdf/DI_count_events.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

tot <- df.indels %>% 
  group_by(sample) %>% 
  summarise(total = n())

p <- df.indels %>% 
  group_by(sample,names) %>% 
  summarise(n = n()) %>%
  left_join(.,tot,by = 'sample') %>% 
  mutate(fraction = n/total) %>% 
  ggplot(., aes(x=sample, y=fraction)) +
  geom_bar(stat="identity",position="dodge") + ylab('Fraction of events') + ylim(0,1) +
  facet_wrap(~names) 

ggsave(filename = 'pdf/DI_fraction_events.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')


# widths

p <- df.indels %>% 
  ggplot(., aes(x=sample, y=tot.width)) + 
  geom_boxplot(varwidth = TRUE) +
  facet_wrap(~names) + ylab('width')

ggsave(filename = 'pdf/DI_widths.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

# locate events

fasta <- readLines('/BCGLAB/darosio_crispr/reference/bwa_reference/egfp.fasta')
fasta <- paste(grep(fasta,pattern = '^>',invert = TRUE,value = TRUE),collapse = "")
fasta <- str_to_upper(fasta)

nchar(fasta)

zb <- df.indels %>% 
  filter(sample != 'Undetermined') %>% 
  mutate(start = start + 1, end = end + 1)

p <- ggplot(zb, aes(x=start)) + 
  geom_histogram(binwidth=1) +
  facet_wrap(~sample,nrow = 5) + coord_cartesian(xlim=c(121,322)) +
  geom_vline(xintercept = c(121,143),size=0.2) +
  geom_vline(aes(xintercept = site), nick,linetype="dotted",size=0.4) + ggtitle('Starting position of D,I events')

ggsave(filename = 'pdf/DI_starts.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- ggplot(zb %>% filter(tot.width > 1), aes(x=start)) + 
  geom_histogram(binwidth=1) +
  facet_wrap(~sample,nrow = 5) + coord_cartesian(xlim=c(121,558)) +
  geom_vline(xintercept = c(121,143,538,558),linetype="dotted",size=0.4) +
  geom_vline(xintercept = 154,color = "#b2182b", size=0.4) +
  geom_vline(xintercept = 202,color = "#006837", size=0.4) + ggtitle('Starting position of D,I events with width > 1')

ggsave(filename = 'pdf/DI_starts_no_single.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

# heatmap

getIndelsCov <- function(zb,min.width=1){
  
  mat <- matrix(data = 0,nrow = 4,ncol = nchar(fasta))
  row.names(mat) <- unique(zb$sample)
  
  for(id in unique(zb$sample)){
    
    message(id)
    
    h <- zb %>% 
      filter(sample == id) %>% 
      filter(tot.width >= min.width)
    
    for(idx in seq_len(nrow(h))){
      
      x <- sort(h$start[idx]:h$end[idx])
      
      mat[id,x] <- mat[id,x] + 1 
      
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
  
  return(mydata)
  
}


covdata <- getIndelsCov(zb = zb,min.width = 1)

p <- ggplot(covdata, aes(x=pos, y=value)) +
  geom_bar(stat="identity",fill='grey50',color='grey50',width = 1) +
  facet_wrap(~sample,nrow = 5) + coord_cartesian(xlim=c(121,558)) +
  geom_vline(xintercept = c(121,143,538,558),linetype="dotted",size=0.4) +
  geom_vline(xintercept = 154,color = "#b2182b", size=0.4) +
  geom_vline(xintercept = 202,color = "#006837", size=0.4) + ggtitle('Positions covered by a D,I events')

ggsave(filename = 'pdf/DI_pos_covered.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

covdata <- getIndelsCov(zb = zb,min.width = 2)

p <- ggplot(data=covdata, aes(x=pos, y=value)) +
  geom_bar(stat="identity",fill='grey50',color='grey50',width = 1) +
  facet_wrap(~sample,nrow = 5) + coord_cartesian(xlim=c(121,558)) +
  geom_vline(xintercept = c(121,143,538,558),linetype="dotted",size=0.4) +
  geom_vline(xintercept = 154,color = "#b2182b", size=0.4) +
  geom_vline(xintercept = 202,color = "#006837", size=0.4) + ggtitle('Positions covered by a D,I events with width > 1')

ggsave(filename = 'pdf/DI_pos_covered_no_singles.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')


}





