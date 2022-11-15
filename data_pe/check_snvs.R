library(tidyverse)
library(ggpubr)
library(gridExtra)

wd <- '/Users/ncasiraghi/Documents/unitn/darosio/results/PE/bwa'

setwd(wd)

nick <- data.frame(sample = unique(df.indels$sample), site = 2954,stringsAsFactors = F)

ps <- list.files(path = "/Users/ncasiraghi/Documents/unitn/darosio/data/pileup",pattern = '\\.pileup',full.names = TRUE)

df <- c()

for(p in ps){
  message(p)
  m <- read.delim(p,stringsAsFactors = FALSE) %>% mutate(sample = gsub(basename(p),pattern = '\\.sorted\\.recal\\.pileup',replacement = '') )
  
  df <- rbind(df,m)
}

# coverage
p <- ggplot(df, aes(x=pos, y=cov)) +
  theme_bw() +
  coord_cartesian(xlim=c(2500,3500)) +
  geom_bar(stat="identity",fill='grey70',color='grey70',width = 1) +
  facet_wrap(~sample,nrow = 4) +
  scale_y_continuous(trans='log10') 

ggsave(filename = 'samples_coverage.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

dat <- df %>%
  filter(cov > 2e5) %>% 
  gather(base,base_cov,4:7) %>% 
  arrange(pos,sample,base) %>% 
  mutate(base_af = base_cov/cov) %>% 
  filter(ref != base) 

p <- ggplot(dat, aes(x=pos, y=base_af, fill=base)) +
  theme_bw() +
  geom_bar(stat="identity",width = 1) +
  coord_cartesian(ylim=c(0,0.01), xlim=c(2900,3000)) +
  geom_vline(xintercept = 2954, linetype="dashed",size=0.5,color = 'orangered') +
  facet_wrap(~sample,nrow = 4) 

ggsave(filename = 'samples_afs_base.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')


# # differences in distribution
# 
# my_comparisons <- list(c('D10A','L-PUC'),
#                        c('L-16','L-PUC'),
#                        c('L-JON','L-PUC'))
# 
# p <- ggboxplot(dat, x = "sample", y = "base_af", color = "sample", palette = "jco") + 
#   stat_compare_means(comparisons = my_comparisons) +
#   theme(axis.text=element_text(size=8))
# 
# ggsave(filename = file.path(wd,'pdf','overall_vafs.pdf'), plot = p, width = 150,height = 150,dpi = 300,units = 'mm',device = 'pdf')
# 
# # check differences around cut sites
# 
# for(wnd in c(10,15,20,25,30)){
#   
#   plots <- list()
#   for(idx in seq_len(nrow(nick))){
#     
#     site <- nick$site[idx]
#     id <- nick$sample[idx]
#     
#     ma <- dat %>% 
#       filter(sample %in% c(id,'L-PUC')) %>% 
#       filter(pos >= (site - wnd)) %>% 
#       filter(pos <= (site + wnd)) %>% 
#       mutate(location = case_when(pos < site ~ 'R1', pos >= site ~ 'R2')) %>% 
#       unite('region',sample,location,remove = FALSE)
#     
#     ma <- ma %>% mutate(region=factor(region, levels=c(paste0(id,'_R1'),paste0(id,'_R2'),'L-PUC_R1','L-PUC_R2')))
#     
#     my_comparisons <- list(c(paste0(id,'_R1'),paste0(id,'_R2')),
#                            c('L-PUC_R1','L-PUC_R2'),
#                            c(paste0(id,'_R1'),'L-PUC_R1'),
#                            c(paste0(id,'_R2'),'L-PUC_R2'))
#     
#     p <- ggboxplot(ma, x = "region", y = "base_af", color = "sample", palette = "jco",title = paste(id,'| window =',wnd)) + 
#       stat_compare_means(comparisons = my_comparisons) +
#       theme(axis.text=element_text(size=8))
#     
#     plots[[idx]] <- p
#     
#   }
#   
#   ggsave(filename = file.path(wd,'pdf',paste0('around_site_wnd',wnd,'.pdf')), plot = do.call(grid.arrange,c(plots,nrow=1)), width = 250,height = 150,dpi = 300,units = 'mm',device = 'pdf')
#   
# }

ctrl <- df %>% 
  filter(sample == 'Negative-control') %>% 
  select(chr,pos,ref,af,cov)

delta <- left_join(df,ctrl,by = c('chr','pos','ref'),suffix = c('_case','_ctrl')) %>% 
  mutate(delta_af = af_case - af_ctrl) %>% 
  filter(cov_case >= 75e4)
  
p <- ggplot(delta, aes(x=pos, y=af_case))+
  geom_line() +
  theme_bw() +
  coord_cartesian(ylim=c(0,0.01), xlim=c(2900,3000)) +
  geom_vline(xintercept=2954, linetype="dashed",size=0.5,color = 'orangered') +
  facet_wrap(~sample,nrow = 4) 

ggsave(filename = 'samples_afs.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- ggplot(delta %>% filter(sample != 'Negative-control'), aes(x=pos, y=delta_af))+
  geom_line() +
  theme_bw() +
  coord_cartesian(xlim=c(2900,3000)) +
  geom_vline(xintercept=2954, linetype="dashed",size=0.5,color = 'orangered') +
  facet_wrap(~sample,nrow = 4) 

ggsave(filename = 'samples_delta_afs.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')
