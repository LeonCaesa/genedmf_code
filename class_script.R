library(data.table)
# library(dmf)
library(MASS)
library(tidyverse)
library(nnet) # multinom
library(rpart) # tree
library(caret) #knn3 with full prob
mainDir = "/projectnb2/dmfgrp/GeneDMF"
setwd(mainDir)
if(!exists("foo", mode="function")) source("utils.R")
if(!exists("foo", mode="function")) source("/projectnb2/dmfgrp/dmf_revision/utils.R")
library(HandTill2001)
# if (! require("remotes")) install.packages("remotes")
# remotes::install_gitlab("fvafrCU/HandTill2001")




data_name = 'data.RData'
n_repeats = 10
#phi_list = c('01', '1', '5', '10')
phi_list = c('0.1')
label_list = c("pop", "super_pop")

auc_df <- data.frame( auc_tree= numeric(0), auc_multi = numeric(0), auc_knn = numeric(0), 
                      y_label = character(0), model = character(0), 
                      chr_num = numeric(0),  
                      repeats = numeric(0))

for (phi in phi_list){
  # model2D_name = paste('pruned', paste('_phi', phi, sep = ''), "_2D.RData", sep = '')
  # model3D_name = paste('pruned', paste('_phi01', sep = ''), "_3D.RData", sep = '')
  model20D_name = paste('pruned', paste('_phi', phi, sep = ''), "_20D.RData", sep = '')
  
  for (chr_num in 1:22){
      # [load files]  
      subDir <- paste("model/chr", chr_num, sep = '')
      load(paste(mainDir, subDir, data_name, sep = '/'))
      
      # [load the factorized results]
      
      # [3D Experiment]
      # load(paste(mainDir, subDir, model3D_name, sep = '/'))
      # load(paste(mainDir, subDir, 'PCA3D.RData', sep = '/'))
      # plot3D = make_regdf(dmf_center(dmf_3D), X)
      # plotPCA3D = make_regdf(pca_center(PCA_3D), X)
      # list.dfs <- list( plot3D,plotPCA3D )
      # model_order = c(paste('dmf3D_phi', phi, sep =''),'PCA3D')
      
      
      # [20D Experiment]
      load(paste(mainDir, subDir, model20D_name, sep = '/'))
      load(paste(mainDir, subDir, 'PCA20D.RData', sep = '/'))
      plot20D = make_regdf(dmf_center(dmf_20D), X)
      plotPCA20D = make_regdf(pca_center(PCA_20D), X)
      list.dfs <- list(plot20D,plotPCA20D)
      model_order = c(paste('dmf20D_phi', phi, sep =''),
                      'PCA20D')
      
      for (j in 1:n_repeats){
        SplitFlag = TrainTest_Flag(plot20D, seed_ = j, train_ratio =0.5)
        # SplitFlag = TrainTest_Flag(plot3D, seed_ = j, train_ratio =0.5)
          for (label_idx in 1:length(label_list)){
              label_ = label_list[label_idx]
              
              # unique_factor = factor(c(unique(plot3D[label_]))[[1]])
              # unique_factor = factor(c(unique(plot20D[label_]))[[1]])
              # plot_level = levels(unique_factor)
              
              unique_factor = c(c(unique(plot20D[label_]))[[1]])
              # unique_factor = factor(c(unique(plot3D[label_]))[[1]])
              confuM = data.frame( matrix(0, nrow = length(unique_factor), ncol = length(unique_factor)))
              plot_level = factor(unique_factor)
              colnames(confuM) = levels(plot_level); rownames(confuM) = levels(plot_level);
              
              
              for (model_idx in 1:length(list.dfs)){
                
                  splitdata = TrainTest_Split(list.dfs[[model_idx]], label_, splitflag = SplitFlag, label_level = levels(plot_level))
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
                  auc_df[nrow(auc_df)+1, ] <- c(round(ht_scores,4), label_, model_order[model_idx], chr_num,  j )
                  
                  print(tail(auc_df,1))
            } #end of four models
        }# end of label
    }# end of repeats
  }# end of chr_num
}#end of phi
#save(auc_df, file = '/projectnb/dmfgrp/GeneDMF/model/class_result/auc_chrRecorded.RData')
#save(auc_df, file = '/projectnb/dmfgrp/GeneDMF/model_updated/class_result/auc20_chrRecorded.RData')
# save(auc_df, file = '/projectnb/dmfgrp/GeneDMF/model_updated/class_result/auc_chrRecorded.RData')

#save(auc_df, file = '/projectnb/dmfgrp/GeneDMF/model_updated05/class_result/auc20_chrRecorded.RData')
# save(auc_df, file = '/projectnb/dmfgrp/GeneDMF/model_updated05/class_result/auc_chrRecorded.RData')

save(auc_df, file = '/projectnb/dmfgrp/GeneDMF/model05_correctmapping/class_result/auc20_chrRecorded.RData')
