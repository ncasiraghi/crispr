library(GenomicAlignments)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

wd <- '/Users/ncasiraghi/Documents/unitn/darosio/results/PE/bwa'

setwd(wd)

load(file = "/Users/ncasiraghi/Documents/unitn/darosio/data/indels/bwa/fulldata.RData")

nick <- data.frame(sample = unique(df.indels$sample), site = 2954,stringsAsFactors = F)

# summary stats

all <- df.totals %>%
  group_by(sample) %>%
  summarise(n = n())

p <- ggplot(data=all, aes(x=sample, y=n)) +
  geom_bar(stat="identity") + ylab('Number of reads') + ggtitle('Selected reads with mapq > 30')

ggsave(filename = 'all_reads.pdf', plot = p, width = 150,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- df.reads %>%
  group_by(sample) %>%
  summarise(n = n()) %>% 
  ggplot(., aes(x=sample, y=n)) +
  geom_bar(stat="identity") + ylab('Number of reads') + ggtitle('Selected reads with D or I in CIGAR [ mapq > 30 ]')

ggsave(filename = 'DI_reads_cnt.pdf', plot = p, width = 150,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- df.reads %>%
  group_by(sample) %>%
  summarise(n = n()) %>% 
  left_join(.,all,by='sample',suffix = c('_reads','_all')) %>% 
  mutate(fraction = n_reads/n_all) %>% 
  ggplot(., aes(x=sample, y=fraction)) +
  geom_bar(stat="identity") + ylab('Fraction of reads') + ggtitle('Selected reads with D or I in CIGAR [ mapq > 30 ]')

ggsave(filename = 'DI_reads_frc.pdf', plot = p, width = 150,height = 150,dpi = 300,units = 'mm',device = 'pdf')

# count events

p <- df.indels %>% 
  group_by(sample,names) %>% 
  summarise(n = n()) %>% 
  ggplot(., aes(x=sample, y=n)) +
  geom_bar(stat="identity",position="dodge") + ylab('Number of events') + 
  facet_wrap(~names) 

ggsave(filename = 'DI_count_events.pdf', plot = p, width = 250,height = 150,dpi = 300,units = 'mm',device = 'pdf')

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

ggsave(filename = 'DI_fraction_events.pdf', plot = p, width = 250,height = 150,dpi = 300,units = 'mm',device = 'pdf')

# widths

p <- df.indels %>% 
  ggplot(., aes(x=sample, y=tot.width)) + 
  geom_boxplot(varwidth = TRUE) +
  facet_wrap(~names) + ylab('width')

ggsave(filename = 'DI_widths_boxplot.pdf', plot = p, width = 250,height = 150,dpi = 300,units = 'mm',device = 'pdf')

# locate events

fasta <- readLines('/Users/ncasiraghi/Documents/unitn/darosio/reference/bwa_reference/spike.fasta')
fasta <- paste(grep(fasta,pattern = '^>',invert = TRUE,value = TRUE),collapse = "")
fasta <- str_to_upper(fasta)

nchar(fasta)

zb <- df.indels %>% mutate(start = start + 1, end = end + 1)

p <- ggplot(zb, aes(x=start)) + 
  geom_histogram(binwidth=1) +
  facet_wrap(~sample,nrow = 5) +
  geom_vline(aes(xintercept=site), nick, linetype="dashed",size=0.5,color = 'orangered') +
  ggtitle('Starting position of D,I events')

ggsave(filename = 'DI_starts.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

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
  
  return(list(mat=mat,mydata=mydata))
  
}

covdata <- getIndelsCov(zb = zb,min.width = 2)

m <- covdata$mat[,2900:3000]

m <- m[all$sample,]

m <- m / all$n

col_fun = colorRamp2(breaks = c(min(m),max(m)), c("white", "forestgreen"))

clab <- str_split(string = fasta,pattern = "",simplify = TRUE)[2900:3000]

column_ha = HeatmapAnnotation(DI = anno_barplot(colSums(m)),border = )

png(filename = "heatmap.png", width = 260,height = 50,res = 300, units = 'mm')

Heatmap(m,name = "DI",
        rect_gp = gpar(col = "black", lwd = 0.5),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col=col_fun,
        column_names_gp = grid::gpar(fontsize = 6),
        column_labels = clab,column_names_centered = TRUE,column_names_rot = 0,top_annotation = column_ha)

dev.off()
