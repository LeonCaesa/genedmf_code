mainDir = "/projectnb2/dmfgrp/dmf_revision/"
#subDir <- paste("model/chr", chr_num, sep = '')
setwd(mainDir)
library(MASS)
if(!exists("foo", mode="function")) source("/projectnb2/dmfgrp/dmf_revision/utils.R")


# [load data-x]
glm_family = negative.binomial(0.1)
load("/projectnb2/dmfgrp/GeneDMF/data/ChrAgg.RData")

# q_star = 3
# PCA_3D = prcomp(data_x, center =FALSE, rank = q_star)
# PCA_3D$V = PCA_3D$rotation
# save(PCA_3D, file = "/projectnb2/dmfgrp/GeneDMF/model/AggChr/Agged_PCA_3D.RData")
#remove(PCA3D)

# q_star = 20
# PCA_20D = prcomp(data_x, center =FALSE, rank = q_star)
# PCA_20D$V = PCA_20D$rotation
# save(PCA_20D, file = "/projectnb2/dmfgrp/GeneDMF/model/AggChr/Agged_PCA_20D.RData")


# q_star = 3
# dmf_fit = dmf(data_x, glm_family, q_star )
# save(dmf_fit, file = "/projectnb2/dmfgrp/GeneDMF/model/AggChr/Agged_phi01_3D.RData")
# remove(dmf_fit)

# q_star = 20
# dmf_fit = dmf(data_x, glm_family, q_star )
# save(dmf_fit, file = "/projectnb2/dmfgrp/GeneDMF/model/AggChr/Agged_phi01_20D.RData")

q_star = 6
PCA_6D = prcomp(data_x, center =FALSE, rank = q_star)
PCA_6D$V = PCA_6D$rotation
save(PCA_6D, file = "/projectnb2/dmfgrp/GeneDMF/model/AggChr/Agged_PCA_6D.RData")

q_star = 6
dmf_fit = dmf(data_x, glm_family, q_star )
save(dmf_fit, file = "/projectnb2/dmfgrp/GeneDMF/model/AggChr/Agged_phi01_6D.RData")



# [not relevant]
# data_name = 'data.RData'
# # n_repeats = 5
# phi = 0.1
# glm_family = negative.binomial(phi)
# label_list = c("pop", "super_pop")

# 
# for (chr_num in 1:22){
#     # [load files]  
#     subDir <- paste("model/chr", chr_num, sep = '')
#     load(paste(mainDir, subDir, data_name, sep = '/'))
#     if (exists('data_x')){data_x = rbind(data_x, X)
#     }else{
#       load(paste(mainDir, subDir, data_name, sep = '/'))
#       data_x = X}
# }

# save(data_x, file = "/projectnb2/dmfgrp/GeneDMF/data/ChrAgg.RData")

    