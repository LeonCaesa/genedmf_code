mainDir = "/projectnb2/dmfgrp/GeneDMF"
#subDir <- paste("model/chr", chr_num, sep = '')
setwd(mainDir)
if(!exists("foo", mode="function")) source("utils.R")
if(!exists("foo", mode="function")) source("/projectnb2/dmfgrp/dmf_revision/utils.R")
library(vcfR)
library(data.table)
library(MASS)
library(tidyverse)
library(plotly)
library(grid)
library(latex2exp)
library(fpc)
library(gridExtra)
library(nnet) # multinom
library(rpart) # tree
library(caret) #knn3 with full prob
library(HandTill2001)

load("/projectnb2/dmfgrp/GeneDMF/data/ChrAgg.RData")

#load( "/projectnb2/dmfgrp/GeneDMF/model/AggChr/Agged_phi01_3D.RData")
load( "/projectnb2/dmfgrp/GeneDMF/model/AggChr/Agged_phi01_20D.RData")
# load( "/projectnb2/dmfgrp/GeneDMF/model/AggChr/Agged_phi01_6D.RData")

#load("/projectnb2/dmfgrp/GeneDMF/model/AggChr/Agged_PCA_3D.RData")
load("/projectnb2/dmfgrp/GeneDMF/model/AggChr/Agged_PCA_20D.RData")
# load("/projectnb2/dmfgrp/GeneDMF/model/AggChr/Agged_PCA_6D.RData")


# count_x = as_tibble(sample(data_x, 10000)) %>% count()
# count_x = table(data_x) 
#plot(count_x)
# par(oma=c( 0,0,0,4))
# png('/projectnb/dmfgrp/GeneDMF/figure/SparseGene.png', units="in", width=5, height=5, res=300)
# image(data_x[sample(1:10000),], col = c("white", "black"),
#       xlab="Individual idx",
#       ylab="Variant idx")
# dev.off()



# plot(sample(data_x, 10000), log='y')


# plot_DMF = make_plotdf(dmf_center(dmf_fit), data_x)
# plot_PCA = make_plotdf(pca_center(PCA_3D), data_x)
plot_DMF = make_regdf(dmf_center(dmf_fit), data_x)
plot_PCA = make_regdf(pca_center(PCA_20D), data_x)
# plot_DMF = make_regdf(dmf_center(dmf_fit), data_x)
# plot_PCA = make_regdf(pca_center(PCA_6D), data_x)


plot_ly(plot_DMF, x=~PC1, y=~PC2, z=~PC3, type="scatter3d",
        mode="markers", color=~pop,
        marker = list(size = 3)) %>%
  layout(legend = list(orientation = "h", 
                       xanchor = "center", x = 0.5),
         margin = list(t = 0, l = 0, r= 0, b =0))%>%
  layout(scene = list(xaxis = list(range = c(-0.04, 0.03),
                                   yaxis = list(range = c(-0.02, 0.04)),
                                   zaxis = list(range = c(-0.04, 0.04)))))

plot_ly(plot_PCA, x=~PC1, y=~PC2, z=~PC3, type="scatter3d",
        mode="markers", color=~pop,
        marker = list(size = 3)) %>%
  layout(legend = list(orientation = "h", 
                       xanchor = "center", x = 0.5),
         margin = list(t = 0, l = 0, r= 0, b =0))%>%
  layout(scene = list(xaxis = list(range = c(-0.04, 0.03),
                      yaxis = list(range = c(-0.02, 0.04)),
                      zaxis = list(range = c(-0.04, 0.04)))))



# [classification] 
n_repeats = 10
list.dfs <- list(plot_DMF, plot_PCA)
#label_list = c( "super_pop","pop")
label_list = c("super_pop")
# label_list = c("pop")

# model_order = c(paste('dmf3D_phi','0.1', sep =''), 'PCA3D')
 model_order = c(paste('dmf20D_phi','0.1', sep =''), 'PCA20sD')
#model_order = c(paste('dmf6D_phi','0.1', sep =''), 'PCA6D')
auc_df <- data.frame( auc_tree= numeric(0), auc_multi = numeric(0), auc_knn = numeric(0), 
                      y_label = character(0), model = character(0), repeats = numeric(0))

# create_confu <- function(pred_prob, response, confuM){
#   pred_label = colnames(pred_prob)[apply(pred_prob, 1, which.max)]
#   conf_tree = confusionMatrix(data = as.factor(pred_label), reference = as.factor(response))$table
#   conf_tree = conf_tree[order(colnames(confuM)),order(colnames(confuM))]  
#   return(conf_tree)}

create_confu <- function(pred_prob, response, confuM, label_level = plot_level){
  pred_label = colnames(pred_prob)[apply(pred_prob, 1, which.max)]
  pred_label = factor(pred_label, levels = levels(label_level))
  response = factor(response, levels = levels(label_level))
  
  conf_return = confusionMatrix(data = pred_label, reference =response)$table
  conf_return = conf_return[order(colnames(confuM)), order(colnames(confuM))]  
  return(conf_return)}

confuM_tree_total = c()
confuM_multi_total = c()
confuM_knn_total= c()

for (label_idx in 1:length(label_list)){
  label_ = label_list[label_idx]
  unique_factor = c(unique(plot_DMF[label_]))[[1]]
  confuM = data.frame( matrix(0, nrow = length(unique_factor), ncol = length(unique_factor)))
  plot_level = factor(unique_factor)
  colnames(confuM) = levels(plot_level); rownames(confuM) = levels(plot_level);
  
  for (model_idx in 1:length(list.dfs)){
      confuM_tree_total[[model_idx]] = confuM
      confuM_multi_total[[model_idx]]  = confuM
      confuM_knn_total[[model_idx]]  = confuM
      for (j in 1:n_repeats){
        SplitFlag = TrainTest_Flag(plot_DMF, seed_ = j, train_ratio = 0.5)
        splitdata = TrainTest_Split(list.dfs[[model_idx]], label_, splitflag = SplitFlag, label_level = levels(plot_level)) # todo: check effect of label_level
        response = splitdata$test$y
        tree_fit = Tree_tuned(splitdata$train)
        multi_fit = multinom(y~., splitdata$train, trace = FALSE)
        knn_fit = knn3(y~., splitdata$train, k =9)
        
        
        tree_pred = predict(tree_fit, newdata = splitdata$test)
        multi_pred = predict(multi_fit, newdata = splitdata$test, type= 'prob')
        knn_pred = predict(knn_fit, newdata = splitdata$test)
        
        htauc_tree = auc(multcap(response = factor(response),predicted = tree_pred))
        htauc_multi = auc(multcap(response = factor(response), predicted = multi_pred))
        htauc_knn = auc(multcap(response = factor(response),predicted = knn_pred))
        ht_scores = c(htauc_tree, htauc_multi, htauc_knn)
        auc_df[nrow(auc_df)+1, ] <- c(round(ht_scores,4), label_, model_order[model_idx],  j )
    
    print(tail(auc_df,1))
    if (label_ == label_list[length(label_list)]){
      confuM_tree_total[[model_idx]]  = confuM_tree_total[[model_idx]] + create_confu(tree_pred, response, confuM_tree_total[[model_idx]], label_level = plot_level)/n_repeats
      confuM_multi_total[[model_idx]]  = confuM_multi_total[[model_idx]] + create_confu(multi_pred, response, confuM_tree_total[[model_idx]], label_level = plot_level)/n_repeats
      confuM_knn_total[[model_idx]]  = confuM_knn_total[[model_idx]] + create_confu(knn_pred, response, confuM_tree_total[[model_idx]], label_level = plot_level)/n_repeats
    }
  } # end of repeats
      
  } #end of models
  
}# end of label

# [Plot AUC]
auc_df[,1:3] = mutate_all(auc_df[,1:3], function(x) as.numeric(as.character(x)))
colnames(auc_df) = c("tree" ,"multi", "knn",   "y_label" ,  "model"  ,  "repeats"  )
auc_long <- auc_df %>%
  gather(class_model, value, tree:knn)
# save(auc_long, file = '/projectnb/dmfgrp/GeneDMF/model_updated05/class_result/auc6D_AllChr.RData')


#png('/projectnb/dmfgrp/GeneDMF/figure/AUC3DAll.png', units="in", width=5, height=3, res=300)
ggplot(auc_long) + geom_boxplot(aes(x = class_model, y= value, colour = model)) +
   theme_bw()+ theme(legend.position="bottom") + 
  ylab('Multi-AUC Score') + xlab('Classification Model') + scale_colour_manual('3D Factr Model', values=c(1,2),
                                                                  labels =c( "DMF", "PCA")) + facet_wrap(~y_label)
#dev.off()



# [plot Confusion Matrix]
model_idx = 1
confuM_lists = list(confuM_tree_total[[model_idx]], confuM_multi_total[[model_idx]], confuM_knn_total[[model_idx]])
class_name = c('tree', 'multi', 'knn')
remove(confuM_load)
for (i in c(1:length(class_name))){
  confuM_load = confuM_lists[[i]]
  confuM_load$predicted = rownames(confuM_load)
  confuM_load <- confuM_load %>% gather(true, value, 1:(length(colnames(confuM_load))-1))  
  confuM_load$class_model =  class_name[i]
  if (exists('confuM_agg')){
    confuM_agg = rbind(confuM_agg, confuM_load)
  }else{
    confuM_agg = confuM_load
  }
}
#save(confuM_agg, file = '/projectnb/dmfgrp/GeneDMF/agg_confusion/correct_label/SuperPop/DMF_AllChrConfu20D.RData')

# [for pop class only]
# mapping = unique(plot_DMF[, c('super_pop', 'pop')])
# confuM_agg = merge(confuM_agg, mapping, by.x = 'true', by.y ='pop' )
# confuM_agg = confuM_agg[order(confuM_agg$super_pop, decreasing = TRUE),]

# [filter]
# confuM_agg = filter(confuM_agg, class_model %in% c('multi'))
# confuM_agg = filter(confuM_agg, class_model %in% c('knn'))

# plot_level = unique(plot_DMF[order(plot_DMF$super_pop, decreasing = FALSE),]$pop)
# confuM_agg = filter(confuM_agg, class_model %in% c('knn'))
# plot_level =  rev(c("PJL", "BEB", "STU", "ITU", "GIH", "GBR",
#                     "FIN", "IBS" ,"CEU", "TSI", "CHS", "CDX" ,
#                     "KHV", "CHB" ,"JPT", "PUR", "CLM", "PEL" ,"MXL",
#                     "ACB", "GWD", "ESN", "MSL", "YRI" ,"LWK", "ASW"))

# plot_level = levels(factor(confuM_agg$true))



  
g1 = ggplot(confuM_agg, aes(x=factor(true,level = plot_level),
                            y=factor(predicted,level = plot_level), 
                            fill=value)) +
  geom_tile() + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="Greens", direction=1) +
  xlab('True') + ylab('Predicted') +  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill= 'none') + # removing legend for `fill`
  labs(title = model_order[1]) + # using a title instead
  # geom_text(size=1.5, aes(label=value), color="black") +facet_wrap(~class_model) +
  geom_text( aes(label=value), color="black") +facet_wrap(~class_model) +
  ggtitle("DMF3D with All Chr")



model_idx = 2
confuM_lists = list(confuM_tree_total[[model_idx]], confuM_multi_total[[model_idx]], confuM_knn_total[[model_idx]])
class_name = c('tree', 'multi', 'knn')
remove(confuM_agg)

for (i in c(1:3)){
  confuM_load = confuM_lists[[i]]
  confuM_load$predicted = rownames(confuM_load)
  confuM_load <- confuM_load %>% gather(true, value, 1:(length(colnames(confuM_load))-1))  
  confuM_load$class_model =  class_name[i]
  if (exists('confuM_agg')){
    confuM_agg = rbind(confuM_agg, confuM_load)
  }else{
    confuM_agg = confuM_load
  }
}
#save(confuM_agg, file = '/projectnb/dmfgrp/GeneDMF/agg_confusion/correct_label/SuperPop/PCA_AllChrConfu20D.RData')
# confuM_agg = filter(confuM_agg, class_model %in% c('multi'))
# confuM_agg = filter(confuM_agg, class_model %in% c('knn'))

# ggplot(confuM_agg, aes(x=true,
#                        y=predicted, 
#                        fill=value)) +
#   geom_tile() + theme_bw() + coord_equal() +
#   scale_fill_distiller(palette="Greens", direction=1) +
#   xlab('True') + ylab('Predicted') +  theme(plot.title = element_text(hjust = 0.5)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   guides(fill= 'none') + # removing legend for `fill`
#   labs(title = model_order[1]) + # using a title instead
#   # geom_text(size=1.5, aes(label=value), color="black") +facet_wrap(~class_model) +
#   geom_text( aes(label=value), color="black") +facet_wrap(~class_model) +
#   ggtitle("DMF3D with All Chr")



g2 = ggplot(confuM_agg, aes(x=factor(true,level = plot_level),
                            y=factor(predicted,level = plot_level), 
                            fill=value)) +
  geom_tile() + theme_bw() + coord_equal() +
  scale_fill_distiller(palette="Greens", direction=1) +
  xlab('True') + ylab('Predicted') +  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill= 'none') + # removing legend for `fill`
  labs(title = model_order[2]) + # using a title instead
  geom_text( aes(label=value), color="black") +facet_wrap(~class_model)+
  # geom_text(size=1.5, aes(label=value), color="black") +facet_wrap(~class_model)+
  ggtitle("PCA3D with All Chr")

# png('/projectnb/dmfgrp/GeneDMF/figure/Pop20DAll_Confu.png', units="in", width=12, height=8, res=300)
# png('/projectnb/dmfgrp/GeneDMF/figure/SPop20DAll_Confu.png', units="in", width=8, height=6, res=300)
png('/projectnb/dmfgrp/GeneDMF/figure/SPop6DAll_Confu.png', units="in", width=8, height=6, res=300)
# png('/projectnb/dmfgrp/GeneDMF/figure/Pop3DAll_ConfuOrdered.png', units="in", width=8, height=4, res=300)
# png('/projectnb/dmfgrp/GeneDMF/figure/Pop3DAll_ConfuOrdered.png', units="in", width=20, height=16, res=300)
# png('/projectnb/dmfgrp/GeneDMF/figure/SPop3DAll_Confu.png', units="in", width=8, height=6, res=300)
# margin = theme(plot.margin = unit(c(0,0.5,0,0.5), "cm"))
margin = theme(plot.margin = unit(c(0,0,0,0), "cm"))
grid.arrange(g1 + margin, g2 + margin, nrow=2)
dev.off()