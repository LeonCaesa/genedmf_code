# [for plotting]
 # load('/projectnb/dmfgrp/GeneDMF/model/class_result/auc_chrRecorded.RData')
#load('/projectnb/dmfgrp/GeneDMF/model/class_result/auc20_chrRecorded.RData')
 
# load('/projectnb/dmfgrp/GeneDMF/model_updated/class_result/auc_chrRecorded.RData')
 # load('/projectnb/dmfgrp/GeneDMF/model_updated/class_result/auc20_chrRecorded.RData')

# load( '/projectnb/dmfgrp/GeneDMF/model_updated05/class_result/auc_chrRecorded.RData')
# load( '/projectnb/dmfgrp/GeneDMF/model_updated05/class_result/auc20_chrRecorded.RData')

load( '/projectnb/dmfgrp/GeneDMF/model05_correctmapping/class_result/auc20_chrRecorded.RData')
library(tidyverse)
library(gridExtra)
library(latex2exp)

auc_df[,1:3] = mutate_all(auc_df[,1:3], function(x) as.numeric(as.character(x)))
colnames(auc_df) = c("tree" ,"multi", "knn",   "y_label" ,  "model"  ,   "chr_num" ,  "repeats"  )
auc_df$chr_num = as.factor(as.numeric( auc_df$chr_num))
auc_long <- auc_df %>%
  gather(class_model, value, tree:knn)



auc_subdf3D = filter(auc_long, model %in% c('dmf20D_phi0.1', 'PCA20D'))
# auc_subdf3D = filter(auc_long, model %in% c('dmf3D_phi0.1', 'PCA3D'))

# png('/projectnb/dmfgrp/GeneDMF/figure/3DAUC.png', units="in", width=8, height=6, res=300)
png('/projectnb/dmfgrp/GeneDMF/figure/20DAUC.png', units="in", width=8, height=6, res=300)
ggplot(auc_subdf3D) + geom_boxplot( outlier.size = 0.1, aes(x = chr_num, y= value, colour = model)) +
  facet_wrap(~ class_model + y_label, ncol = 2,
             labeller = labeller(
               .multi_line = FALSE), scales = "free"
  ) + theme_bw() + theme(legend.position="bottom") + 
  ylab('Multi-AUC Score') + xlab('Chr Num') + 
  # scale_colour_manual('3D Factr Model', values=c(1,2),labels =c(TeX("DMF (NegBin $phi = 0.1$)"), "PCA"))
  scale_colour_manual('20D Factr Model', values=c(1,2),labels =c(TeX("DMF (NegBin $phi = 0.1$)"), "PCA"))
dev.off()

grid.arrange(g1,g2,nrow=2)


best_label <- function(label_name, model_name, class_name, auc_df){
  rank_df = filter(auc_df, y_label == label_name, model == model_name, class_model == class_name)
  rank_df = rank_df %>% group_by(chr_num) %>% summarize(mean_ = mean(value),
                                                        max_ = max(value))
  # return(rank_df)  
  best_row = tail(rank_df[order(rank_df$mean_),],1)
  return (best_row)
}

# best_label('pop', 'dmf3D_phi01', 'knn', auc_subdf2D)
# 
# auc_subdf2D %>% group_by(y_label, model, class_model, chr_num) %>% transform(mean_ = mean(value))

best_label('super_pop', 'dmf20D_phi0.1', 'multi', auc_subdf3D)
best_label('super_pop', 'PCA20D', 'multi', auc_subdf3D)
best_label('pop', 'dmf20D_phi0.1', 'multi', auc_subdf3D)
best_label('pop', 'PCA20D', 'multi', auc_subdf3D)

# super pop, chr = 11/ pop, chr = 2
best_label('super_pop', 'dmf3D_phi0.1', 'multi', auc_subdf3D)
best_label('super_pop', 'PCA3D', 'multi', auc_subdf3D)
best_label('pop', 'dmf3D_phi0.1', 'multi', auc_subdf3D)
best_label('pop', 'PCA3D', 'multi', auc_subdf3D)

