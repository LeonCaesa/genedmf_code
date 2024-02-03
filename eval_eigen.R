mainDir = "/projectnb2/dmfgrp/GeneDMF"
setwd(mainDir)
if(!exists("foo", mode="function")) source("utils.R")


singular_df <- data.frame(pc_idx = character(0), singular = numeric(0), model = character(0), chr_num = numeric(0),max_eigen = numeric(0))


for (chr_num in 1:22){
    chr_folder = paste('chr', chr_num, sep = '')
    data_name = paste( "/projectnb/dmfgrp/GeneDMF/model", chr_folder, 'data.RData', sep = '/')
    dmf_name = paste( "/projectnb/dmfgrp/GeneDMF/model", chr_folder, 'pruned_phi0.1_20D.RData', sep = '/')
    pca_name = paste( "/projectnb/dmfgrp/GeneDMF/model", chr_folder, 'PCA20D.RData', sep = '/')
    load(data_name)
    load(dmf_name)
    load(pca_name)
    
    dmf_centered = dmf_center(dmf_20D)
    pca_centered = pca_center(PCA_20D)
    
    singular_dmf = diag(crossprod(dmf_centered$L))
    singular_pca = diag(crossprod(pca_centered$L))
    
    
    # chr_df = data.frame(rbind(cbind(1:20, singular_dmf,'DMF', chr_num, sum(singular_dmf)),
    #                           cbind(1:20, singular_pca,'PCA', chr_num, sum(singular_pca))))
    
    chr_df = data.frame(rbind(cbind(1:20, cumsum(singular_dmf),'DMF', chr_num, sum(singular_dmf)),
                              cbind(1:20, cumsum(singular_pca),'PCA', chr_num, sum(singular_pca))))
    colnames(chr_df) = colnames(singular_df)
    
    singular_df = rbind(singular_df, chr_df)
}

singular_df$pc_idx =factor(singular_df$pc_idx, levels = c(1:20))
singular_df[, c('singular', 'max_eigen')] = apply(singular_df[, c('singular', 'max_eigen')], 2, as.numeric)

library(ggplot2)
png('/projectnb/dmfgrp/GeneDMF/figure/20D_EigenComp.png', units="in", width=5, height=3, res=300)
ggplot(singular_df) + geom_boxplot(aes(x = pc_idx, y = singular/max_eigen, colour = model)) +
  theme_bw() + theme(legend.position="bottom") + 
  ylab(' % of Explained Variance') + xlab('Num of PCs') + scale_colour_manual('Factr Model', values=c(1,2))
dev.off()





