library(vcfR)
library(data.table)
library(dmf)
library(MASS)
library(tidyverse)
library(plotly)
library(foreach)
library(doParallel)

#vcf <- read.vcfR("ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", nrows = 20000)
#vcf <- read.vcfR("ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")

setwd("/projectnb2/dmfgrp/GeneDMF")

# read vcf
vcf <- read.vcfR("data/chr1.vcf.gz", nrows = 10000)

# genotype matrix 
gt <- extract.gt(vcf, element = 'GT', as.numeric = TRUE, IDtoRowNames = T)

# construct X matrix
X = matrix(gt, nrow = nrow(gt))
colnames(X) = colnames(gt)

#phi = mean(X)^2/(sd(X)^2 - mean(X))
phi = max(0.3,  mean(X)^2/(sd(X)^2 - mean(X)))
factor_family = negative.binomial(phi)


# [2D plot]
dmf_2D = dmf(X, factor_family, rank = 2)

mainDir = "/projectnb2/dmfgrp/GeneDMF"
subDir <- "model/chr1"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
save.image(file = paste(subDir,"10k_2D.RData", sep = '/'))


dmf_3D = dmf(X, factor_family, rank = 3)



# mapping table
chr_df = data.frame(colnames(X))
colnames(chr_df) = 'id'

# file with sample id and populations
samples <- fread("data/integrated_call_samples_v3.20130502.ALL.panel",sep="\t", header = F)
colnames(samples) <- c("id", "pop", "super_pop", "sex")
mapping = merge(chr_df, samples, by.x = "id", by.y = "id")



plotdf = data.frame(dmf_result$V)
colnames(plotdf) = c('PC1', 'PC2')
plotdf$id = mapping$id
plotdf= merge(plotdf, mapping, by.x = "id", by.y = "id")

ggplot(plotdf) + geom_point(aes(x = PC1, y = PC2, colour = pop))


# [3D plot]
dmf_3D = dmf(X, factor_family, rank = 3)

plot3D = data.frame(dmf_center(dmf_3D)$V)
colnames(plot3D) = c('PC1', 'PC2', 'PC3')
plot3D$id = mapping$id
plot3D= merge(plot3D, mapping, by.x = "id", by.y = "id")


# "pop", "super_pop", "sex" 
# [3D plot]
plot_ly(plot3D, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~super_pop,
        marker = list(size = 3)) %>%
  # layout(scene = list(xaxis = list(range = c(-1, 4)),
  #                     yaxis = list(range = c(-1, 2)),
  #                     zaxis = list(range = c(-0.2, 0.6)))
  # )%>%
  layout(legend = list(orientation = "h",   # show entries horizontally
                       xanchor = "center",  # use center of legend as anchor
                       x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))


# [compare to pca]
PCA_3D = prcomp(X, center =TRUE)
plotPCA3D = data.frame(PCA_3D$rotation[,1:3])
colnames(plotPCA3D) = c('PC1','PC2','PC3' )

plotPCA3D$id = mapping$id
plotPCA3D= merge(plotPCA3D, mapping, by.x = "id", by.y = "id")


# "pop", "super_pop", "sex" 
# [3D plot]
plot_ly(plotPCA3D, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~pop,
        marker = list(size = 3)) %>%
  # layout(scene = list(xaxis = list(range = c(-1, 4)),
  #                     yaxis = list(range = c(-1, 2)),
  #                     zaxis = list(range = c(-0.2, 0.6)))
  # )%>%
  layout(legend = list(orientation = "h",   # show entries horizontally
                       xanchor = "center",  # use center of legend as anchor
                       x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))

# [---------------------- fpkm------------------------]
library("DESeq2")

m <- matrix(1e6 * rep(c(.125, .25, .25, .5), each=4),
            ncol=4, dimnames=list(1:4,1:4))
mode(m) <- "integer"
se <- SummarizedExperiment(list(counts=m), colData=DataFrame(sample=1:4))
dds <- DESeqDataSet(se, ~ 1)

# create 4 GRanges with lengths: 1, 1, 2, 2.5 Kb
gr1 <- GRanges("chr1",IRanges(1,1000)) # 1kb
gr2 <- GRanges("chr1",IRanges(c(1,1001),c( 500,1500))) # 1kb
gr3 <- GRanges("chr1",IRanges(c(1,1001),c(1000,2000))) # 2kb
gr4 <- GRanges("chr1",IRanges(c(1,1001),c(200,1300))) # 500bp
rowRanges(dds) <- GRangesList(gr1,gr2,gr3,gr4)

fpkm(dds)