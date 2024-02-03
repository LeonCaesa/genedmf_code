library(foreach)
library(doParallel)



setwd("/projectnb2/dmfgrp/GeneDMF")

# [get the data]

url1 = 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr'
url2 = '.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz'

numCores = detectCores()

registerDoParallel(numCores)  

foreach (chr_num=1:22) %dopar% {
  print(chr_num)
  chr_url = paste(url1, url2,  sep = paste(chr_num))  
  save_dir = paste('data/chr', '.vcf.gz', sep = paste(chr_num))
  download.file(chr_url, save_dir, method = 'wget', quite = TRUE)  
}
stopImplicitCluster()

#dataout = read.csv('data/chr22/chr22_filtered.prune.out')
#datain = read.csv('data/chr22/chr22_filtered.prune.in')