# load GM12878 and HepG2 objects from ERBS package
library(ERBS)
data(GM12878)
data(HepG2)

# inspect HepG2 GRanges object
class(HepG2)
HepG2
values(HepG2)

#access row using indices and bracket
HepG2[1:10, ]

#using $ to access all values of a specific columns
mean(HepG2$signalValue)

median(HepG2$signalValue)

max(HepG2$signalValue)

#find row which has max value
HepG2[which.max(HepG2$signalValue)]

#get all chromosome names
chr = seqnames(HepG2)
as.character(chr)

#count number of regions of chromosome 16
table(chr)[16]

#histogram width
#what is width function?
median(width(HepG2))
hist(width(HepG2))
