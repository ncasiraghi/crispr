library(tidyverse)
library(ggpubr)

setwd('/BCGLAB/darosio_crispr/out/')

nick <- data.frame(sample = c('D10A','L-16','L-JON'), site = c(202,202,154),stringsAsFactors = F)

ps <- list.files('/BCGLAB/darosio_crispr/pacbam',pattern = '\\.filtered.pileup$',full.names = TRUE)

df <- c()

for(p in ps){
  message(p)
  m <- read.delim(p,stringsAsFactors = FALSE) %>% mutate(sample = gsub(basename(p),pattern = '\\.sorted.filtered.pileup',replacement = '') )
  
  df <- rbind(df,m)
}

df <- df %>% filter(sample != 'Undetermined')

# coverage
p <- ggplot(df, aes(x=pos, y=cov)) +
  geom_bar(stat="identity",fill='grey70',color='grey70',width = 1) +
  coord_cartesian(xlim=c(121,558)) +
  geom_vline(xintercept = c(121,143),linetype="dotted",size=0.4) +
  geom_vline(aes(xintercept=site), nick, linetype="dashed",size=0.5,color = 'orangered') +
  facet_wrap(~sample,nrow = 4)

ggsave(filename = 'pdf/samples_coverage.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

dat <- df %>%
  filter(cov >= 75e4) %>% 
  gather(base,base_cov,4:7) %>% 
  arrange(pos,sample,base) %>% 
  mutate(base_af = base_cov/cov) %>% 
  filter(ref != base) %>% 
  mutate(site = case_when(sample %in% c('L-JON','L-PUC') ~ 154, sample %in% c('D10A','L-16') ~ 202)) %>% 
  mutate(location = case_when(sample %in% c('L-JON','L-PUC') & pos <= site ~ 'R1',
                              sample %in% c('L-JON','L-PUC') & pos > site ~ 'R2',
                              sample %in% c('L-16','D10A') & pos <= site ~ 'R1',
                              sample %in% c('L-16','D10A') & pos > site ~ 'R2')) %>% 
  unite('region',sample,location,remove = FALSE)

p <- ggplot(dat, aes(x=pos, y=base_af, fill=base)) +
  geom_bar(stat="identity",width = 1) +
  coord_cartesian(ylim=c(0,0.003), xlim=c(121,322)) +
  geom_vline(xintercept = c(121,143),linetype="dotted",size=0.4) +
  geom_vline(aes(xintercept=site), nick, linetype="dashed",size=0.5,color = 'orangered') +
  facet_wrap(~sample,nrow = 4) 

ggsave(filename = 'pdf/samples_afs_base.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

ctrl <- df %>% 
  filter(sample == 'L-PUC') %>% 
  select(chr,pos,ref,af,cov)

delta <- left_join(df,ctrl,by = c('chr','pos','ref'),suffix = c('_case','_ctrl')) %>% 
  mutate(delta_af = af_case - af_ctrl) %>% 
  filter(cov_case >= 75e4)
  
p <- ggplot(delta, aes(x=pos, y=af_case))+
  geom_line() +
  coord_cartesian(ylim=c(0,0.003), xlim=c(121,322)) +
  geom_vline(xintercept = c(121,143),linetype="dotted",size=0.4) +
  geom_vline(aes(xintercept=site), nick, linetype="dashed",size=0.5,color = 'orangered') +
  facet_wrap(~sample,nrow = 4) 

ggsave(filename = 'pdf/samples_afs.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- ggplot(delta %>% filter(sample != 'L-PUC'), aes(x=pos, y=delta_af))+
  geom_line() +
  coord_cartesian(xlim=c(121,322)) +
  geom_vline(xintercept = c(121,143),linetype="dotted",size=0.4) +
  geom_vline(aes(xintercept=site), nick, linetype="dashed",size=0.5,color = 'orangered') +
  facet_wrap(~sample,nrow = 4) 

ggsave(filename = 'pdf/samples_delta_afs.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- ggplot(delta %>% filter(sample != 'L-PUC'), aes(x=pos, y=delta_af, group=sample))+
  geom_line(aes(color = sample),size = 0.6) +
  coord_cartesian(xlim=c(121,322)) +
  geom_vline(xintercept = c(121,143),linetype="dotted",size=0.4) +
  geom_vline(aes(xintercept=site), nick, linetype="dashed",size=0.5,color = 'orangered')
  
ggsave(filename = 'pdf/samples_delta_afs_alltogether.pdf', plot = p, width = 300,height = 150,dpi = 300,units = 'mm',device = 'pdf')

my_comparisons <- list(c('D10A','L-PUC'),
                       c('L-16','L-PUC'),
                       c('L-JON','L-PUC'))

p <- ggboxplot(dat, x = "sample", y = "base_af", color = "sample", palette = "jco") + 
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 0.003)

ggsave(filename = 'pdf/samples_kruskal_wallis.pdf', plot = p, width = 150,height = 150,dpi = 300,units = 'mm',device = 'pdf')

my_comparisons <- list(c('D10A_R1','D10A_R2'),
                       c('L-16_R1','L-16_R2'),
                       c('L-JON_R1','L-JON_R2'),
                       c('L-PUC_R1','L-PUC_R2'))

p <- ggboxplot(dat %>% filter(!is.na(location)), x = "region", y = "base_af", color = "sample", palette = "jco",order = sort(unique(dat$region))) + 
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 0.003)

ggsave(filename = 'pdf/samples_kruskal_wallis_up_down_sites.pdf', plot = p, width = 250,height = 150,dpi = 300,units = 'mm',device = 'pdf')

# define regions

wnd <- 10

reg <- df %>%
  filter(cov >= 75e4) %>% 
  gather(base,base_cov,4:7) %>% 
  arrange(pos,sample,base) %>% 
  mutate(base_af = base_cov/cov) %>% 
  filter(ref != base) %>% 
  mutate(site = case_when(sample %in% c('L-JON','L-PUC') ~ 154, sample %in% c('D10A','L-16') ~ 202)) %>% 
  mutate(location = case_when(sample %in% c('L-JON','L-PUC') & pos <= site & pos > site - wnd ~ 'R1',
                              sample %in% c('L-JON','L-PUC') & pos <= (site + wnd) & pos > site ~ 'R2',
                              sample %in% c('L-16','D10A') & pos <= site & pos > (site - wnd) ~ 'R1',
                              sample %in% c('L-16','D10A') & pos <= (site + wnd) & pos > site ~ 'R2')) %>% 
  unite('region',sample,location,remove = FALSE)

p <- ggboxplot(reg %>% filter(!is.na(location)), x = "region", y = "base_af", color = "sample", palette = "jco",title = paste('Regions length =',wnd),order = sort(unique(dat$region))) + 
  stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(label.y = 0.003) 

ggsave(filename = 'pdf/samples_kruskal_wallis_up_down_sites_regions.pdf', plot = p, width = 250,height = 150,dpi = 300,units = 'mm',device = 'pdf')


