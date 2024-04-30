if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")
library(impute)
library(dplyr)
library(caret)
library(methylKit)
library(GenomicRanges)
library(genomation)


# load filtered subsets - methylbase objects

load("/scratch/sb61937/work/ZEBRA_SHARK/R/filteredByCov_ouput/zebrashark.5x.meth.destr.53S.RData")

# COVERED IN 53 of 66 samples
zebrashark.5x.percmeth.53S<-percMethylation(zebrashark.5x.meth.destr.53S, rowids = TRUE)
zebrashark.5x.percmeth.53S.df<-data.frame(zebrashark.5x.percmeth.53S)
head(zebrashark.5x.percmeth.53S.df)
colSums(is.na(zebrashark.5x.percmeth.53S.df))

zebrashark.5x.percmeth.53S.df.t<-t(zebrashark.5x.percmeth.53S.df)

zebrashark.5x.percmeth.53S.impute<-impute.knn(zebrashark.5x.percmeth.53S.df.t, k = 5, rowmax = 0.5, colmax = 0.7, rng.seed=362436069)
save(zebrashark.5x.percmeth.53S.impute, file="/scratch/sb61937/work/ZEBRA_SHARK/R/filteredByCov_ouput/zebrashark.5x.percmeth.53S.impute.RData")
zebrashark.5x.percmeth.53S.impute.df<-data.frame(zebrashark.5x.percmeth.53S.impute$data)
head(zebrashark.5x.percmeth.53S.impute.df)
