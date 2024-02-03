mainDir = "/projectnb2/dmfgrp/GeneDMF/"
saveDir = "/projectnb2/dmfgrp/GeneDMF/familytest"
setwd(mainDir)
if(!exists("foo", mode="function")) source("/projectnb/dmfgrp/dmf/R/dmf.R")
library(MASS)

load("/projectnb/dmfgrp/GeneDMF/data/ChrAgg.RData")
rank <- 3
max_data <- max(data_x)

negbin_dmf <- dmf(data_x, rank = rank, family = negative.binomial(0.1))
save(negbin_dmf, paste(saveDir, '/negbinom.RData', sep = ''))

binom_dmf <- dmf(data_x/max_data, rank = rank, family = binomial(), weights = max_data)
save(binom_dmf, paste(saveDir, '/binomial.RData', sep = ''))

gaussian_dmf <- dmf(data_x, rank = rank, family = gaussian())
save(binom_dmf, paste(saveDir, '/gaussian.RData', sep = ''))





  
