setwd("/data/projects/BioinformaticsTutorial/src/microarrays")

rawdata <- affy::ReadAffy()
exprs <- affy::pm(rawdata)
normdata <- affy::rma(rawdata)
normexprs <- affy::exprs(normdata)

design <- cbind(c(1,1,1,1,1,1),c(0,0,0,1,1,1))

samples <- affy::sampleNames(normdata)
probesets <- affy::featureNames(normdata)

fit <- limma::lmFit(normexprs, design)
fit <- limma::eBayes(fit)

#BiocManager::install("annotate")
#BiocManager::install("hgu133plus2.db")

#library(annotate)
library(hgu133plus2.db)

fit.out <- limma::topTable(fit, coef=2, number=20)
fit.out$Symbol <- unlist(mget(rownames(fit.out), hgu133plus2.db::hgu133plus2SYMBOL))



