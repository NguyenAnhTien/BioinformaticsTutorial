BiocManager::install("affy")
library(affy)
rawdata <- affy::ReadAffy()
exprs <- pm(rawdata)
image(rawdata[,1])
length(exprs)
ncol(exprs)
nrow(exprs)

limma::plotMA(exprs[,c(1,3)])