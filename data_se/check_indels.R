library(GenomicAlignments)
library(tidyverse)

wd <- '/BCGLAB/darosio_crispr/bwa/out'
# wd <- '/BCGLAB/darosio_crispr/bowtie2/out'

if(!file.exists(file.path(wd,'pdf'))){
  dir.create(file.path(wd,'pdf'))
}

setwd(wd)

nick <- data.frame(sample = c('D10A','L-16','L-JON'), site = c(202,202,154),stringsAsFactors = F)

## This code was designed for BAM files aligned against a specific reference genome "amplicons-deepseq.fasta"
## BAM files should be first sorted by read position and indexed with SAMTools

load(file = file.path(wd,'rdata','fulldata.RData'))

# summary stats

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

ggsave(filename = 'pdf/DI_widths_boxplot.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

mu <- df.indels %>% 
  group_by(sample,names) %>% 
  summarise(grp.mean = mean(tot.width))

p <- df.indels %>% 
  ggplot(., aes(x=tot.width, color = names)) + 
  geom_histogram(fill="white", position="dodge",binwidth = 1) +
  facet_wrap(~sample,nrow = 5) + ylab('width') +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=names),linetype="dashed")

ggsave(filename = 'pdf/DI_widths_histogram.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

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
  geom_vline(xintercept = c(121,143),linetype="dotted",size=0.4) +
  geom_vline(aes(xintercept=site), nick, linetype="dashed",size=0.5,color = 'orangered') +
  ggtitle('Starting position of D,I events')

ggsave(filename = 'pdf/DI_starts.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

# heatmap

getIndelsCov <- function(zb,min.width=1){
  
  mat <- matrix(data = 0,nrow = 4,ncol = nchar(fasta))
  row.names(mat) <- unique(zb$sample)
  
  for(sn in unique(zb$sample)){
    
    message(sn)
    
    h <- zb %>% 
      filter(sample == sn) %>% 
      filter(tot.width >= min.width)
    
    for(idx in seq_len(nrow(h))){
      
      x <- sort(h$start[idx]:h$end[idx])
      
      mat[sn,x] <- mat[sn,x] + 1 
      
    }
    
  }
  
  mydata <- c()
  
  for(sn in unique(rownames(mat))){
    
    df <- data.frame(sample = sn,
                     value = as.numeric(mat[sn,]),
                     pos = seq_len(ncol(mat)),
                     stringsAsFactors = FALSE)
    
    mydata <- rbind(mydata,df)
    
  }
  
  return(mydata)
  
}

covdata <- getIndelsCov(zb = zb,min.width = 1)

p <- ggplot(covdata, aes(x=pos, y=value)) +
  geom_bar(stat="identity",fill='grey50',color='grey50',width = 1) +
  facet_wrap(~sample,nrow = 5) + coord_cartesian(xlim=c(121,322)) +
  geom_vline(xintercept = c(121,143),linetype="dotted",size=0.4) +
  geom_vline(aes(xintercept=site), nick, linetype="dashed",size=0.5,color = 'orangered') +
  ggtitle('Positions covered by D,I events')

ggsave(filename = 'pdf/DI_pos_covered.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

