if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='http://cran.us.r-project.org')


BiocManager::install("methylKit")
BiocManager::install("GenomicRanges")

library(methylKit)
library(GenomicRanges)
#-------------------------------------------------------------------------------#
## READ IN FILES & PROCESS FILES ##
# First step is to read the sorted bam files into methylkit
file.GN22909 <-list("/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22909_sort.bam")

zebrashark.obj.GN22909 <-processBismarkAln(location=file.GN22909, sample.id= list("GN22909"),
                                     assembly="sSteTig4", save.folder=NULL, save.context=NULL, read.context="CpG",
                                     nolap=FALSE, mincov=1, minqual=20, phred64=FALSE,
                                     treatment=c(1))

save(zebrashark.obj.GN22909, file="/dir/work/ZEBRA_SHARK/R/methylRaw.GN22909.RData")
