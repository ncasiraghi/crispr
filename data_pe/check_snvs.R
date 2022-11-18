library(tidyverse)
library(ggpubr)
library(gridExtra)

wd <- '/mnt/profile/darosio/'

setwd(wd)

ps <- list.files(path = file.path(wd,'data'),pattern = '\\.pileup$',full.names = TRUE,recursive = TRUE)

df <- c()

for(p in ps){
  message(p)
  
  m <- read.delim(p,stringsAsFactors = FALSE) %>% 
    mutate(sample = gsub(basename(p),pattern = '\\.pileup',replacement = '')) %>% 
    mutate(method = dirname(dirname(p)) %>% basename()) %>% 
    filter(cov > 0)
  
  df <- rbind(df,m)
}

df <- df %>% 
  filter(grepl(sample,pattern = '.merged.sorted.clean$')) %>% 
  mutate(sample = gsub(sample,pattern = '.merged.sorted.clean$',replacement = ""))

# coverage
p <- ggplot(df, aes(x=pos, y=cov)) +
  theme_bw() +
  geom_bar(stat="identity",fill='grey70',color='grey70',width = 1) +
  facet_grid(method~sample,scales = 'free_y') +
  ylab("Coverage")

ggsave(filename = 'samples_coverage.pdf', plot = p, width = 400,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- p + scale_y_continuous(trans='log10') + ylab("log10(Coverage)")

ggsave(filename = 'samples_coverage_log10.pdf', plot = p, width = 400,height = 150,dpi = 300,units = 'mm',device = 'pdf')

# variant allelic fractions
dat <- df %>%
  gather(base,base_cov,4:7) %>% 
  arrange(pos,sample,base) %>% 
  mutate(base_af = base_cov/cov) %>% 
  filter(pos > 2954 - 100 & pos < 2954 + 100)

p <- ggplot(dat, aes(x=pos, y=cov)) +
  theme_bw() +
  geom_bar(stat="identity",fill='grey70',color='grey70',width = 1) +
  geom_vline(xintercept = 2954, linetype="dashed",size=0.2,color = 'orangered') +
  facet_grid(method~sample,scales = 'free_y') 

ggsave(filename = 'samples_coverage_zoom.pdf', plot = p, width = 400,height = 150,dpi = 300,units = 'mm',device = 'pdf')

dat <- dat %>% filter(ref != base)

p <- ggplot(dat, aes(x=cov, y=base_af)) +
  geom_point(size=2, shape=23) +
  theme_bw() +
  facet_grid(method~sample) 

ggsave(filename = 'afs_cov.pdf', plot = p, width = 400,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- ggplot(dat %>% filter(cov > 250000), aes(x=pos, y=base_af, fill=base)) +
  theme_bw() +
  geom_bar(stat="identity",width = 1) +
  geom_vline(xintercept = 2954, linetype="dashed",size=0.2,color = 'orangered') +
  facet_grid(method~sample) 

ggsave(filename = 'samples_afs_base.pdf', plot = p, width = 400,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- p + coord_cartesian(ylim = c(0,0.01))

ggsave(filename = 'samples_afs_base_zoom.pdf', plot = p, width = 400,height = 150,dpi = 300,units = 'mm',device = 'pdf')
