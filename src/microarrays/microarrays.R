rawdata <- affy::ReadAffy()
exprs <- affy::pm(rawdata)
normdata <- affy::rma(rawdata)
normexprs <- affy::exprs(normdata)

design <- cbind(c(1,1,1,1,1,1),c(0,0,0,1,1,1))

samples <- affy::sampleNames(normdata)
probesets <- affy::featureNames(normdata)

fit <- limma::lmFit(normexprs, design)
fit <- limma::eBayes(fit)