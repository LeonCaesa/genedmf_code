library(vcfR)
library(data.table)
#library(dmf)
library(MASS)
library(tidyverse)
library(plotly)
library(grid)
library(latex2exp)
library(fpc)
library(gridExtra)
mainDir = "/projectnb2/dmfgrp/GeneDMF"
setwd(mainDir)
if(!exists("foo", mode="function")) source("utils.R")



# model2D_name = "pruned_2D.RData"
# model3D_name = "pruned_3D.RData"

data_name = 'data.RData'


# [visualization] 
make_plotdf = function(dmf_fit, X = X){
  
  plotdf = data.frame(dmf_fit$V)
  col_names =  c('PC1', 'PC2', 'PC3')  
  colnames(plotdf) = col_names[1:dim(plotdf)[2]]
  
  # [mapping table]
  chr_df = data.frame(colnames(X))
  colnames(chr_df) = 'id'
  samples <- fread("data/integrated_call_samples_v3.20130502.ALL.panel",sep="\t", header = F)
  colnames(samples) <- c("id", "pop", "super_pop", "sex")
  mapping = merge(chr_df, samples, by.x = "id", by.y = "id")
  
  plotdf$id = mapping$id
  plotdf= merge(plotdf, mapping, by.x = "id", by.y = "id")
  return(plotdf)
}




# [mapping]
samples <- fread("data/integrated_call_samples_v3.20130502.ALL.panel",sep="\t", header = F)
colnames(samples) <- c("id", "pop", "super_pop", "sex")
label_list = c("pop", "super_pop", "sex")
samples = as.matrix(samples)


CH_df <- data.frame( TwoD= numeric(0), ThreeD= numeric(0), chr = numeric(0), label = numeric(0), phi = character(0))
CH_PCAdf <- data.frame( TwoD= numeric(0), ThreeD= numeric(0), chr = numeric(0), label = numeric(0), phi = character(0))

# phi_list = c('01', '1', '5', '10')
phi_list = c('01')



for (phi in phi_list){
  model2D_name = paste('pruned', paste('_phi', phi, sep = ''), "_2D.RData", sep = '')
  model3D_name = paste('pruned', paste('_phi', phi, sep = ''), "_3D.RData", sep = '')  
      for (chr_num in 1:22){
          # [load files]  
          #chr_num = 1
          subDir <- paste("model/chr", chr_num, sep = '')
          load(paste(mainDir, subDir, data_name, sep = '/'))
          
          # [for dmf]
          load(paste(mainDir, subDir, model2D_name, sep = '/'))
          load(paste(mainDir, subDir, model3D_name, sep = '/'))
          plot2D = make_plotdf(dmf_center(dmf_2D), X)
          plot3D = make_plotdf(dmf_center(dmf_3D), X)
          # plot2D = make_plotdf(dmf_2D, X)
          # plot3D = make_plotdf(dmf_3D, X)
          
          # [for PCA]
          load(paste(mainDir, subDir, 'PCA2D.RData', sep = '/'))
          load(paste(mainDir, subDir, 'PCA3D.RData', sep = '/'))
          plotPCA3D = make_plotdf(pca_center(PCA_3D), X)
          plotPCA2D = make_plotdf(pca_center(PCA_2D), X)
          # plotPCA3D = make_plotdf(PCA_3D, X)
          # plotPCA2D = make_plotdf(PCA_2D, X)
          

          # [CH stats compute]
          for (label_idx in 1:length(label_list)){
              label_ = label_list[label_idx]
              true_label = as.integer(as.factor(unlist(samples[, 1 + label_idx])))
              n_class = dim(unique(samples[, 1 + label_idx]))[1]
              # [for dmf] 
              CH_2D = round(calinhara(plot2D[,2:3], true_label), digits = 2)
              CH_3D = round(calinhara(plot3D[,2:4], true_label), digits = 2) 
              CH_df[nrow(CH_df)+1, ] <- c(CH_2D, CH_3D, chr_num, label_, phi)
              # [for pca] 
              if (phi==phi_list[1]){
                CH_2DPCA = round(calinhara(plotPCA3D[,2:3], true_label), digits = 2)
                CH_3DPCA = round(calinhara(plotPCA3D[,2:4], true_label), digits = 2) 
                CH_PCAdf[nrow(CH_PCAdf)+1, ] <- c(CH_2DPCA, CH_3DPCA, chr_num, label_, phi)}
          }
      }
}


CH_df %>% group_by(label) %>% mutate(max_ = which.max(chr_num))

return_plotdf = function(CH_df){
  CH_df = filter(CH_df, label %in% c('super_pop', 'pop')) 
  CH_df$TwoD = as.numeric(CH_df$TwoD);
  CH_df$ThreeD = as.numeric(CH_df$ThreeD);
  CH_df$chr = as.factor(as.numeric(CH_df$chr))
  CH_df = filter(CH_df, phi %in% c('01') )
  return(CH_df)
}
CH_df = return_plotdf(CH_df);CH_df$model = 'DMF'
CH_PCAdf = return_plotdf(CH_PCAdf);CH_PCAdf$model = 'PCA'
CH_aggdf = rbind(CH_df, CH_PCAdf)
ggplot(CH_aggdf) + geom_point(size = 3, aes(x = chr, y =  log(TwoD), 
                                   colour = model, shape = label)) +
  ggtitle(TeX("2D CH Index Comparision")) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom") +
  ylab('log(CH)') + scale_shape_manual(values = c(0, 1))+ scale_colour_manual('Factr Model', values=c(1,2))


png('/projectnb/dmfgrp/GeneDMF/figure/3DCH.png', units="in", width=8, height=4, res=300)
ggplot(CH_aggdf) + geom_point(size = 3, aes(x = chr, y =  log(ThreeD), 
                                            colour = model, shape = label)) +
  ggtitle(TeX("3D Factorization CH Index Comparision")) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom") +
  xlab('Chr Num') +
  ylab('log(CH)') + scale_colour_manual('Factr Model', values=c(1,2)) + scale_shape_manual('Label', values = c(0, 1)) 
dev.off()



best_label <- function(label_name, dim_num, CH_df){
  rank_df = filter(CH_df, label == label_name)
  if (dim_num ==2){
    best_row = tail(rank_df[order(rank_df$TwoD),], 1)
    }else{
    best_row = tail(rank_df[order(rank_df$ThreeD),], 1)
    }
  return (best_row)
}
DMF_best2d_superpop = best_label('super_pop', 2, CH_df)
DMF_best2d_pop = best_label('pop', 2, CH_df)
DMF_best3d_superpop = best_label('super_pop', 3, CH_df)
DMF_best3d_pop = best_label('pop', 3, CH_df)

PCA_best2d_superpop = best_label('super_pop', 2, CH_PCAdf)
PCA_best2d_pop = best_label('pop', 2, CH_PCAdf)
PCA_best3d_superpop = best_label('super_pop', 3, CH_PCAdf)
PCA_best3d_pop = best_label('pop', 3, CH_PCAdf)



# [prepare the best 2d/3d df]
# [load the best model]
DMF_subDir2d_superpop <- paste("model/chr", DMF_best2d_superpop$chr, sep = '')
DMF_subDir3d_superpop <- paste("model/chr", DMF_best3d_superpop$chr, sep = '')
DMF_subDir2d_pop <- paste("model/chr", DMF_best2d_pop$chr, sep = '')
DMF_subDir3d_pop <- paste("model/chr", DMF_best3d_pop$chr, sep = '')

PCA_subDir2d_superpop <- paste("model/chr", PCA_best2d_superpop$chr, sep = '')
PCA_subDir3d_superpop <- paste("model/chr", PCA_best3d_superpop$chr, sep = '')
PCA_subDir2d_pop <- paste("model/chr", PCA_best2d_pop$chr, sep = '')
PCA_subDir3d_pop <- paste("model/chr", PCA_best3d_pop$chr, sep = '')

model2D_SuperPop_name = paste('pruned', paste('_phi', DMF_best2d_superpop$phi, sep = ''), "_2D.RData", sep = '')
model3D_SuperPop_name = paste('pruned', paste('_phi', DMF_best3d_superpop$phi, sep = ''), "_2D.RData", sep = '')
model2D_Pop_name = paste('pruned', paste('_phi', DMF_best2d_pop$phi, sep = ''), "_2D.RData", sep = '')
model3D_Pop_name = paste('pruned', paste('_phi', DMF_best3d_pop$phi, sep = ''), "_2D.RData", sep = '')

# [fixing name to chr 1]
DMF_subDir3d_superpop = 'model/chr1'
DMF_subDir3d_pop = 'model/chr1'
model3D_SuperPop_name = 'pruned_phi01_3D.RData'
PCA_subDir3d_superpop = 'model/chr1'
PCA_subDir3d_pop = 'model/chr1'
# [end]

# [ loading dmf models]
load(paste(mainDir, DMF_subDir2d_superpop, model2D_SuperPop_name, sep = '/'))
dmf2D_SuperPop = dmf_2D
load(paste(mainDir, DMF_subDir2d_pop, model2D_Pop_name, sep = '/'))
dmf2D_Pop = dmf_2D
load(paste(mainDir, DMF_subDir3d_superpop, model3D_SuperPop_name, sep = '/'))
dmf3D_SuperPop = dmf_3D
load(paste(mainDir, DMF_subDir3d_pop, model3D_Pop_name, sep = '/'))
dmf3D_Pop = dmf_3D

# [ loading pca models]

load(paste(mainDir, PCA_subDir2d_superpop, 'PCA2D.RData' , sep = '/'))
pca2D_SuperPop = PCA_2D
load(paste(mainDir, PCA_subDir3d_superpop, 'PCA3D.RData' , sep = '/'))
pca3D_SuperPop = PCA_3D
load(paste(mainDir, PCA_subDir2d_pop, 'PCA2D.RData' , sep = '/'))
pca2D_Pop = PCA_2D
load(paste(mainDir, PCA_subDir3d_pop, 'PCA3D.RData' , sep = '/'))
pca3D_Pop = PCA_3D

# [load the corresponding data to get column label]
load(paste(mainDir, DMF_subDir2d_superpop, data_name, sep = '/'))
X2d= X
load(paste(mainDir, DMF_subDir3d_superpop, data_name, sep = '/'))
X3d = X


# [load the PCA result]
plot2D_Pop = make_plotdf(dmf_center(dmf2D_Pop), X2d)
plot2D_SuperPop = make_plotdf(dmf_center(dmf2D_SuperPop), X2d)
plot3D_Pop = make_plotdf(dmf_center(dmf3D_Pop), X3d)
plot3D_SuperPop = make_plotdf(dmf_center(dmf3D_SuperPop), X3d)



plotPCA3D_Pop = make_plotdf(pca_center(pca3D_Pop), X3d)
plotPCA3D_SuperPop = make_plotdf(pca_center(pca3D_SuperPop), X3d)
plotPCA2D_Pop = make_plotdf(pca_center(pca2D_Pop), X2d)
plotPCA2D_SuperPop = make_plotdf(pca_center(pca2D_SuperPop), X2d)


# [dmf plotting]
gdmf2d_pop = ggplot(plot2D_Pop) + geom_point(aes(x = PC1, y = PC2, colour = pop)) + 
      ggtitle(TeX("DMF 2D (NegBin $phi = 0.1$)")) + theme(plot.title = element_text(hjust = 0.5))
gdmf2d_superpop = ggplot(plot2D_SuperPop) + geom_point(aes(x = PC1, y = PC2, colour = super_pop))+ 
  ggtitle(TeX("DMF 2D (NegBin $phi = 0.1$)")) + theme(plot.title = element_text(hjust = 0.5))


gdmf3d_pop = plot_ly(plot3D_Pop, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~pop,
                 marker = list(size = 3)) %>%layout(legend = list(orientation = "h",   xanchor = "center", x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))
gdmf3d_superpop = plot_ly(plot3D_SuperPop, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~super_pop,
                     marker = list(size = 3)) %>%layout(legend = list(orientation = "h",   xanchor = "center", x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))

# [pca plotting]
gpca2d_pop = ggplot(plotPCA2D_Pop) + geom_point(aes(x = PC1, y = PC2, colour = pop)) + 
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position="bottom", legend.justification = c(0.5), plot.title = element_text(hjust = 0.5)) + ggtitle('PCA 2D')

gpca2d_superpop = ggplot(plotPCA2D_SuperPop) + geom_point(aes(x = PC1, y = PC2, colour = super_pop)) + 
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position="bottom", legend.justification = c(0.5), plot.title = element_text(hjust = 0.5)) + ggtitle('PCA 2D')




gpca3d_pop = plot_ly(plotPCA3D_Pop, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~pop,
                     marker = list(size = 3)) %>%layout(legend = list(orientation = "h",   xanchor = "center", x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))
gpca3d_superpop = plot_ly(plotPCA3D_SuperPop, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~super_pop,
                          marker = list(size = 3)) %>%
            layout(legend = list(orientation = "h",
                                 xanchor = "center", x = 0.5),
                   margin = list(t = 0, l = 0, r= 0, b =0))


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend <- get_legend(gpca2d_pop)
grid.arrange(arrangeGrob(gdmf2d_pop + theme(legend.position="none"),
                         gpca2d_pop + theme(legend.position="none"),
                         nrow=1),
             legend, nrow=2,heights=c(10, 1))


legend <- get_legend(gpca2d_superpop)
grid.arrange(arrangeGrob(gdmf2d_superpop + theme(legend.position="none"),
                         gpca2d_superpop + theme(legend.position="none"),
                         nrow=1),
             legend, nrow=2,heights=c(10, 1))


#[END OF ANALYSIS]

# best_2d = tail(CH_df[order(CH_df$TwoD),],1)
# best_3d = tail(CH_df[order(CH_df$ThreeD),],1)
# 
# best_2dPCA = tail(CH_PCAdf[order(CH_PCAdf$TwoD),],1)
# best_3dPCA = tail(CH_PCAdf[order(CH_PCAdf$ThreeD),],1)


# [prepare the best 2d/3d df]
  # [load the best model]
#   subDir2d <- paste("model/chr", best_2d$chr, sep = '')
#   subDir3d <- paste("model/chr", best_3d$chr, sep = '')
#   
#   subDir2dPCA <- paste("model/chr", best_2dPCA$chr, sep = '')
#   subDir3dPCA <- paste("model/chr", best_3dPCA$chr, sep = '')
#   
#   model2D_name = paste('pruned', paste('_phi', best_2d$phi, sep = ''), "_2D.RData", sep = '')
#   model3D_name = paste('pruned', paste('_phi', best_3d$phi, sep = ''), "_3D.RData", sep = '')  
#   
#   load(paste(mainDir, subDir2d, model2D_name, sep = '/'))
#   load(paste(mainDir, subDir3d, model3D_name, sep = '/'))
#   # [load the corresponding data]
#   load(paste(mainDir, subDir2d, data_name, sep = '/'))
#   X2d= X
#   load(paste(mainDir, subDir3d, data_name, sep = '/'))
#   X3d = X
#   
#   # load(paste(mainDir, subDir2dPCA, data_name, sep = '/'))
#   # X2dPCA= X
#   # load(paste(mainDir, subDir3dPCA, data_name, sep = '/'))
#   # X3dPCA = X
#   
#   
#   # [load the PCA result]
#   load(paste(mainDir, subDir2dPCA, 'PCA2D.RData', sep = '/'))
#   load(paste(mainDir, subDir3dPCA, 'PCA3D.RData', sep = '/'))
#   
#   # [prepare the data]
#   # plot2D = make_plotdf(dmf_2D, X2d)
#   # plot3D = make_plotdf(dmf_3D, X3d)
#   # plotPCA3D = make_plotdf(PCA_3D, X3d)
#   # plotPCA2D = make_plotdf(PCA_2D, X2d)
#   
#   plot2D = make_plotdf(dmf_center(dmf_2D), X2d)
#   plot3D = make_plotdf(dmf_center(dmf_3D), X3d)
#   plotPCA3D = make_plotdf(pca_center(PCA_3D), X3d)
#   plotPCA2D = make_plotdf(pca_center(PCA_2D), X2d)
#   
# #  "pop", "super_pop", "sex"
# gdmf2d = ggplot(plot2D) + geom_point(aes(x = PC1, y = PC2, colour = super_pop))
# gdmf3d = plot_ly(plot3D, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~pop,
#                  marker = list(size = 3)) %>%layout(legend = list(orientation = "h",   xanchor = "center", x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))
# 
# 
# 
# gpca2d = ggplot(plotPCA2D) + geom_point(aes(x = PC1, y = PC2, colour = super_pop))
# gpca3d = plot_ly(plotPCA3D, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~pop,
#       marker = list(size = 3)) %>%layout(legend = list(orientation = "h",   xanchor = "center", x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))
# 
# # [estimate the phi]
# # eta_2d = tcrossprod(dmf_2D$L, dmf_2D$V) 
# # mu_2d = dmf_2D$family$linkinv(eta_2d)
# # phi2d = theta.ml(mu_2d, X2d)
# # 
# # mu_3d = dmf_3D$family$linkinv(tcrossprod(dmf_3D$L, dmf_3D$V))
# # phi3d = theta.ml(mu_3d, X3d)
# 
# 
# # [agg the plot]
# library(gridExtra)
# grid.arrange(gpca2d, gdmf2d, ncol=2)
# 
# print(best_2d)
# print(best_3d)
# 
# print(best_2dPCA)
# print(best_3dPCA)
# # subplot(gpca3d, gdmf3d, nrows =2)
# gdmf2d = ggplot(plot3D) + geom_point(aes(x = PC1, y = PC3, colour = super_pop))
# gpca2d = ggplot(plotPCA3D) + geom_point(aes(x = PC1, y = PC3, colour = super_pop))
# grid.arrange(gpca2d, gdmf2d, ncol=2)

