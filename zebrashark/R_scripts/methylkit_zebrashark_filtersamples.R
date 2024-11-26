library(dplyr)
library(caret)
library(methylKit)
library(GenomicRanges)
library(genomation)

#-------------------------------------------------------------------------------#
## FILTER BY COVERAGE ##

# Filter by coverage
load("/scratch/sb61937/work/ZEBRA_SHARK/R/filtered.zebrashark.obj.5x.RData")

# NORMALIZE COVERAGE
filtered.zebrashark.obj.5x.norm <- normalizeCoverage(filtered.zebrashark.obj.5x)

# Merge all samples
# creates MethylBase object with CpGs covered by at least 10x in ALL samples
# destranding can increase the number of CpGs available for analysis in some cases by providing better coverage of individual CpGs
# this function merges reads on both strands of a CpG dinucleotide
# setting destrand=TRUE is only recommended when looking at CpG methylation
# this function will only work when operating on a base-pair resolution

# CpGs covered in ALL individuals
zebrashark.5x.meth.destr.all <- unite(filtered.zebrashark.obj.5x.norm, destrand=TRUE)
dim(zebrashark.5x.meth.destr.all)
save(zebrashark.5x.meth.destr.all, file = "zebrashark.5x.meth.destr.all.RData")

# CpGs covered in AT LEAST 30 SAMPLES
zebrashark.5x.meth.destr.30S<- unite(filtered.zebrashark.obj.5x.norm, min.per.group=30L, destrand=TRUE)
dim(zebrashark.5x.meth.destr.30S)
save(zebrashark.5x.meth.destr.30S, file = "zebrashark.5x.meth.destr.30S.RData")

# CpGs covered in AT LEAST 28 SAMPLES
zebrashark.5x.meth.destr.28S<- unite(filtered.zebrashark.obj.5x.norm, min.per.group=28L, destrand=TRUE)
dim(zebrashark.5x.meth.destr.28S)
save(zebrashark.5x.meth.destr.28S, file = "zebrashark.5x.meth.destr.28S.RData")

# CpGs covered in AT LEAST 21 SAMPLES
zebrashark.5x.meth.destr.21S<- unite(filtered.zebrashark.obj.5x.norm, min.per.group=21L, destrand=TRUE)
dim(zebrashark.5x.meth.destr.21S)
save(zebrashark.5x.meth.destr.21S, file = "zebrashark.5x.meth.destr.21S.RData")
