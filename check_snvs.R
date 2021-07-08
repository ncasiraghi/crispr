library(tidyverse)

setwd('/BCGLAB/darosio_crispr/out/')

ps <- list.files('/BCGLAB/darosio_crispr/pacbam',pattern = '\\.pileup$',full.names = TRUE)

df <- c()

for(p in ps){
  message(p)
  m <- read.delim(p,stringsAsFactors = FALSE) %>% mutate(sample = gsub(basename(p),pattern = '\\.sorted.pileup',replacement = '') )
  
  df <- rbind(df,m)
}

head(df)

m <- df %>% 
  filter(sample != 'Undetermined') %>% 
  filter(cov > 5000) %>% 
  gather(base,base_cov,4:7) %>% 
  arrange(pos,sample,base) %>% 
  mutate(base_af = base_cov/cov) %>% 
  filter(base != ref)

ctrl <- m %>% filter(sample == 'L-PUC')
  
case <- m %>% filter(sample != 'L-PUC')
  
case$pt.pvalue <- NA

for(i in seq_len(nrow(case))){
  
  a <- ctrl %>% 
    filter(pos == case$pos[i]) %>% 
    filter(base == case$base[i])
  
  # permutation test
  
  x <- c(a$base_cov, case$base_cov[i])
  
  n <- c(a$cov, case$cov[i])
  
  pt <- prop.test(x=x,n=n)
  
  case$pt.pvalue[i] <- pt$p.value
  
}

case$p.adjust <- p.adjust(case$pt.pvalue,method = 'fdr')

m <- case %>% filter(p.adjust <= 0.01)

m %>% 
  group_by(sample) %>% 
  summarise(n = n_distinct(pos))

ggplot(data=m, aes(x=pos, y=cov)) +
  geom_bar(stat="identity",fill='grey50',color='grey50',width = 1) +
  facet_wrap(~sample,nrow = 4) +
  geom_vline(xintercept = 154,color = "#b2182b", size=0.4) +
  geom_vline(xintercept = 202,color = "#006837", size=0.4)

ggplot(data=m, aes(x=pos, y=base_af, fill=base)) +
  geom_bar(stat="identity",width = 1) +
  coord_cartesian(ylim=c(0,0.003)) +
  facet_wrap(~sample,nrow = 4) +
  geom_vline(xintercept = 154,color = "#b2182b", size=0.4) +
  geom_vline(xintercept = 202,color = "#006837", size=0.4)

ggplot(data=m, aes(x=pos, y=base_af, fill=base)) +
  geom_bar(stat="identity",width = 1) +
  coord_cartesian(ylim=c(0,0.003)) +
  facet_wrap(~sample,nrow = 4) +
  geom_vline(xintercept = 154,color = "#b2182b", size=0.4) +
  geom_vline(xintercept = 202,color = "#006837", size=0.4)

