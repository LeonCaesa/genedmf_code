library(data.table)
#library(dmf)
library(MASS)
library(tidyverse)
library(nnet) # multinom
library(rpart) # tree
library(caret) #knn3 with full prob
library(HandTill2001)
mainDir = "/projectnb2/dmfgrp/GeneDMF"
setwd(mainDir)
if(!exists("foo", mode="function")) source("utils.R")
if(!exists("foo", mode="function")) source("/projectnb2/dmfgrp/dmf_revision/utils.R")

# if (! require("remotes")) install.packages("remotes")
# remotes::install_gitlab("fvafrCU/HandTill2001")



# subDir <- paste("model/chr", chr_num, sep = '')
create_confu <- function(pred_prob, response, confuM, label_level = plot_level){
  pred_label = colnames(pred_prob)[apply(pred_prob, 1, which.max)]
  pred_label = factor(pred_label, levels = levels(label_level))
  response = factor(response, levels = levels(label_level))
  
  conf_return = confusionMatrix(data = pred_label, reference =response)$table
  conf_return = conf_return[order(colnames(confuM)), order(colnames(confuM))]  
  return(conf_return)}


data_name = 'data.RData'
n_repeats = 10
#phi_list = c('01')
phi_list = c('0.1')
label_list = c("pop", "super_pop")

save_name = 'model05_correctmapping'
auc_df <- data.frame( auc_tree= numeric(0), auc_multi = numeric(0), auc_knn = numeric(0), 
                      y_label = character(0), model = character(0), 
                      chr_num = numeric(0),  
                      repeats = numeric(0))

dir.create(paste(mainDir, save_name, sep = '/'), showWarnings = FALSE)
dir.create(paste(mainDir, save_name, 'confuM_result', sep = '/'), showWarnings = FALSE)

# plot_level =  rev(c("PJL", "BEB", "STU", "ITU", "GIH", "GBR",
#                     "FIN", "IBS" ,"CEU", "TSI", "CHS", "CDX" ,
#                     "KHV", "CHB" ,"JPT", "PUR", "CLM", "PEL" ,"MXL",
#                     "ACB", "GWD", "ESN", "MSL", "YRI" ,"LWK", "ASW"))

for (phi in phi_list){
  model20D_name = paste('pruned', paste('_phi', phi, sep = ''), "_20D.RData", sep = '')
  #model3D_name = paste('pruned', paste('_phi01', sep = ''), "_3D.RData", sep = '')
  for (chr_num in 1:22){
    # [load files]  
    subDir <- paste("model/chr", chr_num, sep = '')
    load(paste(mainDir, subDir, data_name, sep = '/'))
    
    # [load the factorized results]
    # load(paste(mainDir, subDir, model3D_name, sep = '/'))
    # load(paste(mainDir, subDir, 'PCA3D.RData', sep = '/'))
    load(paste(mainDir, subDir, model20D_name, sep = '/'))
    load(paste(mainDir, subDir, 'PCA20D.RData', sep = '/'))
  
      # [center and create regression df]
    # [3D Experiment]
    # plot3D = make_regdf(dmf_center(dmf_3D), X)
    # plotPCA3D = make_regdf(pca_center(PCA_3D), X)
    # list.dfs <- list( plot3D, plotPCA3D)
    # model_order = c(paste('dmf3D_phi', phi, sep =''),
    #                 'PCA3D')
    
    # [20D Experiment]
    plot20D = make_regdf(dmf_center(dmf_20D), X)
    plotPCA20D = make_regdf(pca_center(PCA_20D), X)
    list.dfs <- list(plot20D,plotPCA20D)
    model_order = c(paste('dmf20D_phi', phi, sep =''),
                    'PCA20D')
    
    
    
    for (label_idx in 1:length(label_list)){
        label_ = label_list[label_idx]
      
        dir.create(paste(mainDir, save_name, 'confuM_result', label_, sep = '/'), showWarnings = FALSE)
        dir.create(paste(mainDir, save_name, 'confuM_result', label_, chr_num, sep = '/'), showWarnings = FALSE)

        unique_factor = c(c(unique(plot20D[label_]))[[1]])
        # unique_factor = factor(c(unique(plot3D[label_]))[[1]])
        confuM = data.frame( matrix(0, nrow = length(unique_factor), ncol = length(unique_factor)))
        plot_level = factor(unique_factor)
        colnames(confuM) = levels(plot_level); rownames(confuM) = levels(plot_level);
        
      for (model_idx in 1:length(list.dfs)){
        
          confuM_tree_total = confuM
          confuM_multi_total = confuM
          confuM_knn_total = confuM
          
          dir.create(paste(mainDir, save_name, 'confuM_result', label_, chr_num, model_order[model_idx], sep = '/'), showWarnings = FALSE)        
        
          for (j in 1:n_repeats){
              # SplitFlag = TrainTest_Flag(plot3D, seed_ =j, train_ratio =0.5)
              SplitFlag = TrainTest_Flag(plot20D, seed_ =j, train_ratio = 0.5)
        
              splitdata = TrainTest_Split(list.dfs[[model_idx]], label_, splitflag = SplitFlag, label_level = levels(plot_level))
              response = splitdata$test$y
              
              tree_fit = Tree_tuned(splitdata$train)
              multi_fit = multinom(y~., splitdata$train, trace = FALSE)
              knn_fit = knn3(y~., splitdata$train, k =9)
        
              
              # create_confu <- function(pred_prob, response, confuM, label_level = plot_level){
              #   pred_label = colnames(pred_prob)[apply(pred_prob, 1, which.max)]
              #   pred_label = factor(pred_label, levels = label_level)
              #   response = factor(response, levels = label_level)
              # 
              #   conf_return = confusionMatrix(data = pred_label, reference =response)$table
              #   conf_return = conf_return[order(colnames(confuM)), order(colnames(confuM))]  
              #   return(conf_return)}
              

          tree_pred = predict(tree_fit, newdata = splitdata$test)
          multi_pred = predict(multi_fit, newdata = splitdata$test, type= 'prob')
          knn_pred = predict(knn_fit, newdata = splitdata$test)
          
          confuM_tree_total = confuM_tree_total + create_confu(tree_pred, response, confuM_tree_total, label_level = plot_level)/n_repeats
          confuM_multi_total = confuM_multi_total + create_confu(multi_pred, response, confuM_tree_total, label_level = plot_level)/n_repeats
          confuM_knn_total = confuM_knn_total + create_confu(knn_pred, response, confuM_tree_total, label_level = plot_level)/n_repeats
          
          print( c(chr_num, label_, j))
        } # end of repeats
      save_dir = paste(mainDir, save_name, 'confuM_result', label_, chr_num, sep = '/') 
            
      save_tree= paste(save_dir, paste(model_order[model_idx], '/tree.RData', sep = ''), sep = '/')
      save_multi= paste(save_dir, paste(model_order[model_idx], '/multi.RData', sep = ''), sep = '/')
      save_knn= paste(save_dir, paste(model_order[model_idx], '/knn.RData', sep = ''), sep = '/')
      save(confuM_tree_total, file = save_tree)
      save(confuM_multi_total, file = save_multi)
      save(confuM_knn_total, file = save_knn)

      } #end of four models
    }# end of label
  }# end of chr_num
}#end of phi



n_lists = rep(0, 22)
for (chr_num in 1:22){
  # [load files]  
  subDir <- paste("model/chr", chr_num, sep = '')
  load(paste(mainDir, subDir, data_name, sep = '/'))
  n_lists[chr_num] = dim(X)[1]
}

n_df = data.frame(n_lists); colnames(n_df) = 'gene_count'
n_df$chr_num = c(1:22)
png('/projectnb/dmfgrp/GeneDMF/figure/gene_count.png', units="in",  width=5, height=5, res=300)
ggplot(n_df, aes(x =chr_num ,y= gene_count)) + xlab('Chr Num') + ylab('Variant Count') + 
  theme_bw() + geom_bar(stat="identity") + theme(plot.margin = margin(2, 0, 2, 0, "cm"))
dev.off()