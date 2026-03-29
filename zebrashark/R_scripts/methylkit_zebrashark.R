if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='http://cran.us.r-project.org')


BiocManager::install("methylKit")
BiocManager::install("GenomicRanges")

library(methylKit)
library(GenomicRanges)
#-------------------------------------------------------------------------------#
## READ IN FILES & PROCESS FILES ##
# First step is to read the sorted bam files into methylkit
file.list.zebrashark <- list("/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22646_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22647_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22648_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22649_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22650_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22651_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22652_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22867_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22868_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22869_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22870_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22871_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22872_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22873_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22874_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22875_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22876_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22877_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22878_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22879_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22880_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22881_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22882_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22883_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22884_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22885_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22886_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22887_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22888_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22889_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22890_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22891_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22892_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22893_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22894_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22895_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22896_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22897_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22898_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22899_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22900_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22901_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22902_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22903_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22904_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22905_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22906_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22907_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22908_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22909_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22910_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22911_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22912_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22913_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22914_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22915_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22916_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22917_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22918_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22920_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22921_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22922_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22923_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22926_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22927_sort.bam",
                            "/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/GN22928_sort.bam")

zebrashark.obj <-processBismarkAln(location=file.list.zebrashark, sample.id= list("GN22646", "GN22647", "GN22648", "GN22649", "GN22650", "GN22651", "GN22652",
                                                                                  "GN22867", "GN22868", "GN22869", "GN22870", "GN22871", "GN22872", "GN22873",
                                                                                  "GN22874", "GN22875", "GN22876", "GN22877", "GN22878", "GN22879", "GN22880",
                                                                                  "GN22881", "GN22882", "GN22883", "GN22884", "GN22885", "GN22886", "GN22887",
                                                                                  "GN22888", "GN22889", "GN22890", "GN22891", "GN22892", "GN22893", "GN22894",
                                                                                  "GN22895", "GN22896", "GN22897", "GN22898", "GN22899", "GN22900", "GN22901",
                                                                                  "GN22902", "GN22903", "GN22904", "GN22905", "GN22906", "GN22907", "GN22908",
                                                                                  "GN22909", "GN22910", "GN22911", "GN22912", "GN22913", "GN22914", "GN22915",
                                                                                  "GN22916", "GN22917", "GN22918", "GN22920", "GN22921", "GN22922", "GN22923",
                                                                                  "GN22926", "GN22927", "GN22928"),
                                  assembly="sSteTig4", save.folder=NULL, save.context=NULL, read.context="CpG",
                                  nolap=FALSE, mincov=1, minqual=20, phred64=FALSE,
                                  treatment=c(1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1, 1, 1, 1, 1,
                                              1, 1, 1))

#-------------------------------------------------------------------------------#
## FILTER BY COVERAGE ##

# Filter by coverage
# Remove CpGs with <1x or <5x read coverage, and CpGs in the 99.9%tile of coverage (PCR bias)

filtered.zebrashark.obj.1x <- filterByCoverage(zebrashark.obj, lo.count=1, lo.perc=NULL, hi.count=NULL, high.perc=99.9)
str(filtered.zebrashark.obj.1x)

save(filtered.zebrashark.obj.1x, file = "filtered.zebrashark.obj.1x.RData")

filtered.zebrashark.obj.5x <- filterByCoverage(zebrashark.obj, lo.count=5, lo.perc=NULL, hi.count=NULL, high.perc=99.9)
str(filtered.zebrashark.obj.5x)

save(filtered.zebrashark.obj.5x, file = "filtered.zebrashark.obj.5x.RData")

rm(list = ls(all.names = TRUE))
