(!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='http://cran.us.r-project.org')


BiocManager::install("methylKit")
BiocManager::install("GenomicRanges")

library(methylKit)
library(GenomicRanges)

# remove old objects from environment
rm(list=ls())

load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22646.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22647.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22648.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22649.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22650.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22651.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22652.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22867.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22868.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22869.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22870.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22871.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22872.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22873.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22874.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22875.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22876.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22877.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22878.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22879.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22880.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22881.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22882.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22883.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22884.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22885.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22886.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22887.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22888.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22889.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22890.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22891.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22892.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22893.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22894.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22895.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22896.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22897.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22898.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22899.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22900.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22901.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22902.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22903.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22904.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22905.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22906.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22907.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22908.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22909.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22910.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22911.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22912.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22913.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22914.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22915.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22916.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22917.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22918.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22920.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22921.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22922.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22923.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22926.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22927.RData")
load("/scratch/sb61937/work/ZEBRA_SHARK/R/methylRaw_files/methylRaw.GN22928.RData")



zebrashark.obj.collate <- methylRawList(c(zebrashark.obj.GN22646, zebrashark.obj.GN22647, zebrashark.obj.GN22648, zebrashark.obj.GN22649, zebrashark.obj.GN22650, zebrashark.obj.GN22651,
                               zebrashark.obj.GN22652, zebrashark.obj.GN22867, zebrashark.obj.GN22868, zebrashark.obj.GN22869, zebrashark.obj.GN22870, zebrashark.obj.GN22871,
                               zebrashark.obj.GN22872, zebrashark.obj.GN22873, zebrashark.obj.GN22874, zebrashark.obj.GN22875, zebrashark.obj.GN22876, zebrashark.obj.GN22877,
                               zebrashark.obj.GN22878, zebrashark.obj.GN22879, zebrashark.obj.GN22880, zebrashark.obj.GN22881, zebrashark.obj.GN22882, zebrashark.obj.GN22883,
                               zebrashark.obj.GN22884, zebrashark.obj.GN22885, zebrashark.obj.GN22886, zebrashark.obj.GN22887, zebrashark.obj.GN22888, zebrashark.obj.GN22889,
                               zebrashark.obj.GN22890, zebrashark.obj.GN22891, zebrashark.obj.GN22892, zebrashark.obj.GN22893, zebrashark.obj.GN22894, zebrashark.obj.GN22895,
                               zebrashark.obj.GN22896, zebrashark.obj.GN22897, zebrashark.obj.GN22898, zebrashark.obj.GN22899, zebrashark.obj.GN22900, zebrashark.obj.GN22901,
                               zebrashark.obj.GN22902, zebrashark.obj.GN22903, zebrashark.obj.GN22904, zebrashark.obj.GN22905, zebrashark.obj.GN22906, zebrashark.obj.GN22907,
                               zebrashark.obj.GN22908, zebrashark.obj.GN22909, zebrashark.obj.GN22910, zebrashark.obj.GN22911, zebrashark.obj.GN22912, zebrashark.obj.GN22913,
                               zebrashark.obj.GN22914, zebrashark.obj.GN22915, zebrashark.obj.GN22916, zebrashark.obj.GN22917, zebrashark.obj.GN22918, zebrashark.obj.GN22920,
                               zebrashark.obj.GN22921, zebrashark.obj.GN22922, zebrashark.obj.GN22923, zebrashark.obj.GN22926, zebrashark.obj.GN22927, zebrashark.obj.GN22928),
                               treatment=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1))


filtered.zebrashark.obj.1x <- filterByCoverage(zebrashark.obj, lo.count=10, lo.perc=NULL, hi.count=NULL, high.perc=99.9)
str(filtered.zebrashark.obj.1x)

save(filtered.zebrashark.obj.1x, file = "filtered.zebrashark.obj.1x.RData")

# remove objects from environment
rm(list=ls())
