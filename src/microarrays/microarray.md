# Microarray with R
## Dependency

```
install.packages("BiocManager")
BiocManager::install("affy")
BiocManager::install("limma")
```

## Prepare Data

### Downloading data
* http://www.ncbi.nlm.nih.gov/geo/
* Enter GSE25191 into the field *keyword or GEO Accession* and press enter button
* or using the link https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25191
* Uncompress the downloaded file *GSE25191_RAW.tar* into the same folder containing the code files

## Read Data
```
rawdata <- affy::ReadAffy() #read all .CEL.gz files
```

## Extract expression data
```
exprs <- affy::pm(rawdata)
```
The number of columns of exprs is the number of .CEL.gz files.

## Data Visualization
## Plot MA plot
plot the first and the third chips
```
limma::plotMA(exprs[,c(1,3)])
```

## Scatterplots
```
plot(x=exprs[,1],y=exprs[,3])
```

## Compare 1st, 3th, and 6th chips
```
plot(x=exprs[,1],y=exprs[,6],col="blue")
points(x=exprs[,1],y=exprs[,3])
```

## Box plots
```
boxplot(exprs)
```
![Box Plot](boxplot.png)

*Box plot of microarray expression*

```
boxplot(log(exprs,base=2))
```

![Box Plot of Log Values](logarithmboxplot.png)

*Box plot of Logarithm values of the microarray expression*


## Normalizing data
```
normdata <- affy::rma(rawdata)
normexprs <- exprs(normdata)
```

## R Matrix
```
design <- cbind(c(1, 1, 1, 1, 1, 1), c(0, 0, 0, 1, 1, 1))
design
#[, 1] [, 2]
#[1, ] 1 0
#[2, ] 1 0
#[3, ] 1 0
#[4, ] 1 1
#[5, ] 1 1
#[6, ] 1 1
```

## Differential expression

### Design Matrix
```
design <- cbind(c(1,1,1,1,1,1),c(0,0,0,1,1,1))
```
* The first column is intercept
* The second column indicates whether the corresponding microarray profiled samples from side-populations (value = 1) or non-side-population (value = 0)

### Fitting linear models
```
rawdata <- affy::ReadAffy()
exprs <- affy::pm(rawdata)
normdata <- affy::rma(rawdata)
normexprs <- affy::exprs(normdata)
design <- cbind(c(1,1,1,1,1,1),c(0,0,0,1,1,1))
fit <- limma::lmFit(normexprs, design)
fit <- limma::eBayes(fit)
```

#### Annother example for Design Matrix

```
design2 <- cbind(intercept=1,pair2=c(0,1,0,0,1,0), pair3=c(0,0,1,0,0,1))
```

* The second column encodes does this sample have the variable 'Pair2'?
* The third column encodes does this sample have the variable 'Pair3'?

#### Fitting linear models with The New Design Matrix

```
fit2 <- limma::lmFit(normexprs,design2)
fit2 <- limma::eBayes(fit2)
```

#### Noting
The number of rows of the **design** matrix must equal to the number of columns of the **normexprs**.


### Interpreting Results
#### Target

* Identify lists of probesets with significant differential expression

#### Platform Annotation

* Must have annotation information for the microarray platform that has been used for the experiment
* Annotation information can be download from http://www.bioconductor.org/packages/release/data/annotation/
* The annotation files have the name ending in .db
* Ex: hgu133plus2.db provides detailed information about the hgu133plus2 platform

#### Loading Annotation

**Must install the following package**

```
BiocManager::install("annotate")
BiocManager::install("hgu133plus2.db")
```

**Using topTable command to take the annotated table of differentially-expressed genes**

```
library(hgu133plus2.db)

fit.out <- limma::topTable(fit, coef=2, number=20)
fit.out$Symbol <- unlist(mget(rownames(fit.out), hgu133plus2.db::hgu133plus2SYMBOL))
```

* fit.out is a data frame
* mget functions look up the gene symbol annotations for each rowname in the fit.out, then add this as an additional column.

**Diggin the few of rows of the fit.out**

```
fit.out[1:10, c(1, 7, 2:5)]
```

```
              logFC    Symbol  AveExpr         t      P.Value   adj.P.Val
244829_at   -3.512223 LINC00518 4.280104 -24.25504 1.294201e-07 0.007076045
201667_at   -4.290522      GJA1 7.151514 -12.56086 8.528587e-06 0.121794106
241079_at    2.113882      <NA> 4.534789  12.33570 9.551411e-06 0.121794106
1560422_at   1.938365      <NA> 4.240167  11.57664 1.419890e-05 0.121794106
204115_at    1.489710     GNG11 7.129100  11.41084 1.553265e-05 0.121794106
223875_s_at  2.110507      EPC1 6.656454  11.16331 1.780061e-05 0.121794106
207086_x_at -1.401677      <NA> 6.445216 -11.12425 1.819208e-05 0.121794106
207663_x_at -1.865800      <NA> 4.392493 -11.10017 1.843839e-05 0.121794106
241855_s_at  1.408529      <NA> 6.137816  10.67595 2.347218e-05 0.121794106
214844_s_at -2.815481      DOK5 6.797450 -10.53924 2.541796e-05 0.121794106
```

* The first column is the **ID** - this is the identifier for the probe-set from which the corresponding measurements come.
* **Symbol** - the official symbol for the gene into which the probe-set maps.
* **logFC** - the log fold-change associated with the contrast, this is equivalent to the coefficients of the linear model.

**Export fit.out to file**

```
write.table(fit.out,file="fittable.txt",sep="",quote=FALSE,row.names=FALSE)
```