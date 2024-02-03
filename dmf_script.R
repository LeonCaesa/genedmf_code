library(vcfR)
library(data.table)
library(dmf)
library(MASS)
library(tidyverse)
library(foreach)
library(doParallel)
mainDir = "/projectnb2/dmfgrp/GeneDMF"
setwd(mainDir)


#nrows = 10000
# model2D_name = paste('pruned', "_2D.RData", sep = '')
# model3D_name = paste('pruned', "_3D.RData", sep = '')



#numCores = detectCores()
#registerDoParallel(numCores)  

#vcftools --vcf data/chr22.vcf.gz --indv 22_16050607 --extract-FORMAT-info GT --out data/chr22/test.vcf.gz

#phi_list = c(1, 5, 10)
#phi_list = c(10)
phi_list = c(0.1)
q_star = 20
for (phi in phi_list){
      #print(phi)
      # model2D_name = paste('pruned_phi', phi, "_2D.RData", sep = '')
      # model3D_name = paste('pruned_phi', phi, "_3D.RData", sep = '')
      model20D_name = paste('pruned_phi', phi, "_20D.RData", sep = '')

#foreach (chr_num=1:22) %dopar% {
      for (chr_num in 1:22){
        print(c(phi,chr_num))
        #data_dir = paste("data/chr", ".vcf.gz",  sep = paste(chr_num))  
        data_dir = paste("data/chr", '/chr', '_pruned.vcf.gz', sep = paste(chr_num))  
        # [read the data]
        #vcf <- read.vcfR(data_dir, nrows = nrows)
        vcf <- read.vcfR(data_dir)
        gt <- extract.gt(vcf, element = 'GT', as.numeric = TRUE, IDtoRowNames = T)
        X = matrix(gt, nrow = nrow(gt));colnames(X) = colnames(gt)
        
        # [estimate dispersion]
        #phi_hat = mean(X)^2/(sd(X)^2 - mean(X))
        #phi = min(0.3,  max(phi_hat, 0.1))
        factor_family = negative.binomial(phi)
        
        # [model fitting]
        # dmf_2D = dmf(X, factor_family, rank = 2)
        # dmf_3D = dmf(X, factor_family, rank = 3)
        dmf_20D = dmf(X, factor_family, rank = q_star)
      
        # [model saving]
        subDir <- paste("model/chr", chr_num, sep = '')
        dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
        #save.image(file = paste(subDir, model_name , sep = '/'))
        #save(dmf_2D, file = paste(subDir, model2D_name , sep = '/'))
        #save(dmf_3D, file = paste(subDir, model3D_name , sep = '/'))
        save(dmf_20D, file = paste(subDir, model20D_name , sep = '/'))
        
        #if (!file.exists(paste(subDir, 'PCA3D.RData', sep = '/'))){
            # [pca fitting]
            # PCA_3D = prcomp(X, center =FALSE, rank = 3)
            # PCA_3D$V = PCA_3D$rotation
            # PCA_2D = prcomp(X, center =FALSE, rank = 2)
            # PCA_2D$V = PCA_2D$rotation
            # save(PCA_3D, file = paste(subDir, 'PCA3D.RData', sep = '/'))
            # save(PCA_2D, file = paste(subDir, 'PCA2D.RData', sep = '/'))  
            # save(X, file  = paste(subDir, 'data.RData', sep = '/'))
            
            PCA_20D = prcomp(X, center =FALSE, rank = q_star)
            PCA_20D$V = PCA_20D$rotation
            
            save(PCA_20D, file = paste(subDir, 'PCA20D.RData', sep = '/'))
            save(X, file  = paste(subDir, 'data.RData', sep = '/'))
        #}
      }
}
#stopImplicitCluster()
