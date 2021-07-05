library(GenomicAlignments)

## This code was designed for BAM files aligned against a specific reference genome "amplicons-deepseq.fasta"
## BAM files should be first sorted by read position and indexed with SAMTools

targets = read.delim("PATH/amplicons-deepseq.fasta.fai",header=F)
targets$cut.pos = c(58,63,81,76,59,60,100,72,48,66,47,57,67,80,75,85,64,63,68,68,54,61,78,59,60,58,61,57,63,59,48)
window=5

## Function to get the INDELs events
getIndels <- function(cigar,cut.pos)
{
  tab = as.data.frame(cigarRangesAlongReferenceSpace(cigar, with.ops=TRUE)[[1]])
  x = as.data.frame(cigarRangesAlongQuerySpace(cigar, with.ops=TRUE)[[1]])
  tab$tot = x$width+tab$width
  res = 0
  tmp = tab[which(tab$names%in%c("D")&tab$tot>=2),]
  if(nrow(tmp)>0)
  {
    tmp = tmp[which((tmp$start>=cut.pos-window&tmp$end<=cut.pos+window)|
                      (tmp$start<=cut.pos-window&tmp$end>=cut.pos+window)|
                      (tmp$start<=cut.pos-window&tmp$end>=cut.pos-window&tmp$end<=cut.pos+window)|
                      (tmp$end>=cut.pos+window&tmp$start>=cut.pos-window&tmp$start<=cut.pos+window)),]
    res = res + nrow(tmp)
  }
  tmp = tab[which(tab$names%in%c("I")&tab$tot>=2),]
  if(nrow(tmp)>0)
  {
    tmp = tmp[which(tmp$start>=cut.pos-window&tmp$start<=cut.pos+window),]
    res = res + nrow(tmp)
  }
  tmp = tab[which(tab$names%in%c("I","D")&tab$tot==1),]
  if(nrow(tmp)>0)
  {
    tmp = tmp[which(tmp$start>=cut.pos&tmp$start<=cut.pos+1),]
    res = res + nrow(tmp)
  }
  return(res)
}

### Code to get the coverage and mapping specificity
l = list.files("BAM_files_path/",pattern = ".sorted.bam$")
tab = targets[,1,drop=F]
index=2
for(i in l)
{
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=NA),what=c("rname", "pos", "cigar","mapq"))
  bam <- scanBam(paste("BAM_files_path/",i,sep=""), param=param)[[1]]
  tot.reads = as.data.frame(bam)
  tot.reads = tot.reads[which(tot.reads$mapq>=30),]
  
  vals = c()
  for(k in 1:nrow(tab))
  {
    reads = tot.reads[which(tot.reads$rname==tab[k,1]),]
    vals = c(vals,nrow(reads))
  }
  
  tab = cbind(tab,vals)
  colnames(tab)[index] = gsub(".sorted.bam","",i)
  index=index+1
  cat(i,"\n")
}


## Code to extract INDELs statistics from BAM files
cores=10
l = list.files("BAM_files_path/",pattern = ".sorted.bam$") 
tab = targets[,1,drop=F]
index=2
for(i in l)
{
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=NA),what=c("rname", "pos", "cigar","mapq"))
  bam <- scanBam(paste("BAM_files_path/",i,sep=""), param=param)[[1]]
  tot.reads = as.data.frame(bam)
  tot.reads = tot.reads[which(tot.reads$mapq>=30),]
  
  vals = c()
  for(k in 1:nrow(tab))
  {
    reads = tot.reads[which(tot.reads$rname==tab[k,1]),]
    tot = nrow(reads)
    reads = reads[grep("D|I",reads$cigar),]
    
    if(nrow(reads)>0)
    {  
      indels = mclapply(1:nrow(reads),function(j)
      {
        cut.pos = targets$cut.pos[which(targets[,1]==reads$rname[j])]
        return(getIndels(reads$cigar[j],cut.pos))
      },mc.cores=cores)
      
      vals = c(vals,length(which(unlist(indels)>0))/tot)
    } else
    {
      vals = c(vals,0)
    }
  }
  
  ### edit table
  tab = cbind(tab,vals)
  colnames(tab)[index] = gsub(".info","",i)
  index=index+1
  cat(i,"\n")
}

