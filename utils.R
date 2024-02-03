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
  
  plotdf$id = colnames(X)
  plotdf= merge(plotdf, mapping, by.x = "id", by.y = "id")
  return(plotdf)
}

make_regdf <- function(dmf_fit, X){
  
  plotdf = data.frame(dmf_fit$V)
  colnames(plotdf) = paste('PC', 1:dim(plotdf)[2], sep = '')
  
  # [mapping table]
  chr_df = data.frame(colnames(X))
  colnames(chr_df) = 'id'
  samples <- fread("data/integrated_call_samples_v3.20130502.ALL.panel",sep="\t", header = F)
  colnames(samples) <- c("id", "pop", "super_pop", "sex")
  mapping = merge(chr_df, samples, by.x = "id", by.y = "id")
  
  #plotdf$id = mapping$id
  plotdf$id = colnames(X)
  plotdf= merge(plotdf, mapping, by.x = "id", by.y = "id")
  
}


pca_center = function(pca_result){
  pca_result$L = pca_result$x
  pca_result$family = gaussian()
  pca_result = dmf_center(pca_result)
  return(pca_result)
}

TrainTest_Flag<-function(plot_df, seed_, train_ratio =0.7){
  train_size = floor(train_ratio * nrow(plot_df))
  set.seed(seed_)
  split1<- sample(c(rep(0, train_size), rep(1, nrow(plot_df)- train_size)))
  list(train = split1==0, test = split1==1)}

TrainTest_Split<-function(plot_df, label_name, splitflag, label_level){
  if (is.element( 'PC20', colnames(plot_df))){
    # subset_cols = c(c('PC1', 'PC2', 'PC3'), label_name)}
    subset_cols = c(paste('PC', 1:20, sep = ''), label_name)}
  else{
    subset_cols = c(c('PC1', 'PC2', 'PC3'), label_name)}
  reg_df = plot_df[subset_cols]
  p = dim(reg_df)[2]; colnames(reg_df)[p] = 'y'
  reg_df$y = factor(reg_df$y, levels= label_level)
  
  # return (list(train = reg_df[splitflag$train, ], test = reg_df[!splitflag$test, ]))}
  return (list(train = reg_df[splitflag$train, ], test = reg_df[splitflag$test, ]))}

Tree_tuned<-function(data){
  tree_fit = rpart(y ~ ., data = data, parms = list(split ="information"), 
                   control = c(cp = 0, xval =10), method= 'class') 

  best_cp = tree_fit$cptable[which.min(tree_fit$cptable[, "xerror"]),"CP"]

  tree_fit = prune(tree_fit, best_cp)
  return(tree_fit)}

dmf_center <- function (lv, x0 = rep(1, nrow(lv$L)), reorder = TRUE) {
  q0 <- qr(x0)
  Lc <- qr.resid(q0, lv$L); Vc <- lv$V
  if (reorder) {
    so <- order(apply(Lc, 2, norm2), decreasing = TRUE)
    Lc <- Lc[, so]; Vc <- Vc[, so]
  }
  sl <- svd(Lc) # re-orthogonalize
  Lc <- sweep(sl$u, 2, sl$d, `*`)
  Vc <- Vc %*% sl$v
  list(L = Lc, V = Vc, deviance = lv$deviance, family = lv$family,
       center = drop(tcrossprod(lv$V, qr.coef(q0, lv$L))))
}

zeros <- function (x) matrix(0, nrow = nrow(x), ncol = ncol(x))
symmetrize <- function (x) .5 * (x + t(x))
norm2 <- function (x) norm(as.matrix(x[!is.na(x)]), "2")
normalize <- function (x, margin = 2)
  sweep(x, margin, apply(x, margin, norm2), `/`)