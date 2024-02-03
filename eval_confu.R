library(data.table)
#library(dmf)
library(MASS)
library(tidyverse)
library(latex2exp)
# library(nnet) # multinom
# library(rpart) # tree
# library(caret) #knn3 with full prob
# library(HandTill2001)
library(gridExtra)
mainDir = "/projectnb2/dmfgrp/GeneDMF"
setwd(mainDir)
if(!exists("foo", mode="function")) source("utils.R")


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

return_plotdf = function(root_dir, label_name, chr_num, model_factor, model_class){
  
  load_dir = paste(root_dir, label_name, chr_num, model_factor, model_class, sep = '/')
  confuM_load <-loadRData(load_dir)
  confuM_load$predicted = rownames(confuM_load)
  confuM_load <- confuM_load %>% gather(true, value, 1:(length(colnames(confuM_load))-1))
  return (confuM_load)
}


return_plotallmodel = function(root_dir, label_name, chr_num, model_factor){
  
  model_list = c('tree', 'multi', 'knn')

  for (i in 1:length(model_list)){
      if (i ==1){
         confuM_loadAgg = return_plotdf(root_dir, label_name, chr_num, model_factor, paste(model_list[i], '.RData', sep =''))
         confuM_loadAgg$model = model_list[i]
      }else{
        confuM_load = return_plotdf(root_dir, label_name, chr_num, model_factor, paste(model_list[i], '.RData', sep =''))
        confuM_load$model = model_list[i]
        confuM_loadAgg = rbind(confuM_loadAgg, confuM_load)
      }
  }
  return (confuM_loadAgg)
}
# [In Sample Result]
#root_dir = '/projectnb/dmfgrp/GeneDMF/model/confuM_result/random_split'
# root_dir = '/projectnb/dmfgrp/GeneDMF/model/confuM_result/seed_split'
# root_dir = '/projectnb/dmfgrp/GeneDMF/model/confuM_result/20D_seedsplit/'

# [Out of Sample Result]
# root_dir = '/projectnb/dmfgrp/GeneDMF/model_updated/confuM_result/3D_seedsplit/'
# root_dir = '/projectnb/dmfgrp/GeneDMF/model_updated/confuM_result/20D_seedsplit/'

# [Out of Sample Result ratio = 0.5]
# root_dir = '/projectnb/dmfgrp/GeneDMF/model_updated05/confuM_result/3D_seedsplit/'
# root_dir = '/projectnb/dmfgrp/GeneDMF/model_updated05/confuM_result/20D_seedsplit/'
ThreeD = TRUE
if (ThreeD){
    root_dir = '/projectnb/dmfgrp/GeneDMF/model05_correctmapping/confuM_result/3D_seedsplit/'
    chr_num = 2
    model_factor = 'dmf3D_phi0.1'
    comp_factor = 'PCA3D'
    plot_PCAtitle= paste('PCA3D with Chr ', chr_num, sep ='')
    plot_DMFtitle= paste('DMF3D with Chr ', chr_num, sep ='')
    
}else{
    root_dir = '/projectnb/dmfgrp/GeneDMF/model05_correctmapping/confuM_result/20D_seedsplit/'
    chr_num = 1
    model_factor = 'dmf20D_phi0.1'
    comp_factor = 'PCA20D'
    plot_PCAtitle= paste('PCA20D with Chr ', chr_num, sep ='')
    plot_DMFtitle= paste('DMF20D with Chr ', chr_num, sep ='')
}


# [uncomment for right population class]
# label_name = 'pop'
label_name = 'super_pop'


# [depreciated for single classification model plot]
# model_class = 'knn.RData'
# model_class = 'multi.RData'
# model_class = 'tree.RData'
# dmf_load = return_plotdf(root_dir, label_name, chr_num, model_factor, model_class)
# pca_load = return_plotdf(root_dir, label_name, chr_num, comp_factor, model_class)


dmf_load = return_plotallmodel(root_dir, label_name, chr_num, model_factor)
pca_load = return_plotallmodel(root_dir, label_name, chr_num, comp_factor)

if (label_name == 'pop'){
  plot_level =  rev(c("PJL", "BEB", "STU", "ITU", "GIH", "GBR",
                      "FIN", "IBS" ,"CEU", "TSI", "CHS", "CDX" ,
                      "KHV", "CHB" ,"JPT", "PUR", "CLM", "PEL" ,"MXL",
                      "ACB", "GWD", "ESN", "MSL", "YRI" ,"LWK", "ASW"))  
}else{
  plot_level = levels(factor(dmf_load$true))  
}

# ggplot(dmf_load, aes(x=true,
#                      y=predicted,
#                      fill=value)) +
#   geom_tile() + theme_bw() + coord_equal() +
#   scale_fill_distiller(palette="Greens", direction=1) +
#   xlab('True') + ylab('Predicted') +  theme(plot.title = element_text(hjust = 0.5)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   guides(fill=F) +
#   facet_wrap(~model) +
#   labs(title = plot_DMFtitle) + # using a title instead
#   geom_text( aes(label=value), color="black") # printing values


# [change to percentage]
dmf_load = dmf_load %>% group_by(true, model) %>% mutate(totl = sum(value))
dmf_load$value = round(dmf_load$value/dmf_load$totl,3) * 100
pca_load = pca_load %>% group_by(true, model) %>% mutate(totl = sum(value))
pca_load$value = round(pca_load$value/pca_load$totl,3) * 100

g1 = ggplot(dmf_load, aes(x=factor(true,level = plot_level),
                          y=factor(predicted,level = plot_level), 
                          fill=value)) +
  geom_tile() + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="Greens", limits = c(0, 200), direction=1) +
  xlab('True') + ylab('Predicted') +  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=F) + 
  facet_wrap(~model) +
  labs(title = plot_DMFtitle) + # using a title instead
  geom_text( aes(label=value), color="black") # printing values
  # geom_text(size=1.5, aes(label=value), color="black") # printing values

g2 = ggplot(pca_load, aes(x=factor(true,level = plot_level),
                          y=factor(predicted,level = plot_level), 
                          fill=value)) +
  geom_tile() + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="Greens",limits = c(0, 200), direction=1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab('True') + ylab('Predicted') +  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill=F) + 
  facet_wrap(~model) +
  labs(title = plot_PCAtitle) + # using a title instead
  geom_text(aes(label=value), color="black") # printing values
  # geom_text(size=1.5,aes(label=value), color="black") # printing values

# sum(colSums(data.frame(matrix(dmf_load$value, nrow = 26, ncol = 26))))



# [for single model plot]
# png('/projectnb/dmfgrp/GeneDMF/figure/Pop20D_ConfuOrdered.png', units="in", width=8, height=4, res=300)
# png('/projectnb/dmfgrp/GeneDMF/figure/Pop20D_ConfuOrdered.png', units="in", width=8, height=4, res=300)

# [for all three model plot]
  # [for Pop class]
# png('/projectnb/dmfgrp/GeneDMF/figure/Pop20D_ConfuOrdered.png', units="in", width=20, height=16, res=300)
# png('/projectnb/dmfgrp/GeneDMF/figure/Pop3D_ConfuOrdered.png', units="in", width=20, height=16, res=300)
  
  # [for Spop class]
png('/projectnb/dmfgrp/GeneDMF/figure/SPop3D_Confu.png', units="in", width=8, height=6, res=300)
#png('/projectnb/dmfgrp/GeneDMF/figure/SPop20D_Confu.png', units="in", width=8, height=6, res=300)
grid.arrange(g1, g2, nrow=2)
dev.off()




# [To check that the class distribution is similar to the original]
model_class = 'multi'
pca_load  = filter(pca_load, model == model_class)
dmf_load  = filter(dmf_load, model == model_class)
sample_classes = colSums(matrix(pca_load$value, nrow = length(plot_level)))
sample_classes = data.frame(sample_classes)
  # [to check which class DMF performs badly]
  diff_ = round(diag(matrix(dmf_load$value, nrow = length(plot_level)) - matrix(pca_load$value, nrow = length(plot_level))),4)
  names(diff_) = unique(dmf_load$true)
  diff_[diff_<0]

  # [to check the number of sampled pop/super_pop classes] 
  cbind(unique(pca_load$true), sample_classes * 10) #multiply by 10 since repeated 10 times


  # [ to check the percentage of sampled pop/super_pop classes]
  samples <- fread("data/integrated_call_samples_v3.20130502.ALL.panel",sep="\t", header = F)
  colnames(samples) = c('idx', 'pop', 'super_pop', 'sex')
  if (label_name == 'pop'){
    sample_classes$pop =  unique(pca_load$true)
    comp_table = samples %>% group_by(pop) %>% summarise(true = n())  
    comp_table = merge(sample_classes, comp_table, by = 'pop')
  }else{
    sample_classes$super_pop =  unique(pca_load$true)
    comp_table = samples %>% group_by(super_pop) %>% summarise(true = n())  
    comp_table = merge(sample_classes, comp_table, by = 'super_pop')
  }
  comp_table$sample = comp_table$sample_classes / sum(comp_table$sample_classes)
  comp_table$true = comp_table$true / sum(comp_table$true)
  
  # the classes are roughly in proportion, but not exactly due to random sample
  plot(comp_table$true, comp_table$sample)




# [for pheatpmap analysis]
kmeans_pca = data.frame(matrix(pca_load$value, nrow = length(plot_level)))
colnames(kmeans_pca) = unique(pca_load$true)
rownames(kmeans_pca) = unique(pca_load$predicted)

kmeans_dmf = data.frame(matrix(dmf_load$value, nrow = length(plot_level)))
colnames(kmeans_dmf) = unique(dmf_load$true)
rownames(kmeans_dmf) = unique(dmf_load$predicted)


pheatmap::pheatmap(kmeans_dmf)
pheatmap::pheatmap(kmeans_dmf)









