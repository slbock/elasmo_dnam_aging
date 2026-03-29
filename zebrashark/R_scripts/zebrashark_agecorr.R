library(impute)
library(dplyr)
library(caret)
library(methylKit)
library(GenomicRanges)
library(genomation)
library(WGCNA)
install.packages("reshape2", repos='http://cran.us.r-project.org', dependencies = TRUE)
library(reshape2)

load("/dir/work/ZEBRA_SHARK/R/filteredByCov_ouput/zebrashark.5x.percmeth.53S.impute.RData")

zebrashark.5x.percmeth.53S.impute.data <- zebrashark.5x.percmeth.53S.impute$data
zebrashark.5x.percmeth.53S.impute.data.df <- data.frame(zebrashark.5x.percmeth.53S.impute.data)
rownames(zebrashark.5x.percmeth.53S.impute.data.df)



# subset dataframe to only include captive individuals
# wild individuals to be removed: GN22646, GN22648, GN22876, GN22882, GN22886, GN22890, GN22896, GN22899, GN22901, GN22903, GN22908, GN22911, GN22912, GN22914, GN22915, GN22916, GN22922, GN22923

captive <- c("GN22647", "GN22649", "GN22650", "GN22651", "GN22652", "GN22867",
             "GN22868", "GN22869", "GN22870", "GN22871", "GN22872", "GN22873",
             "GN22874", "GN22875", "GN22877", "GN22878", "GN22879", "GN22880",
             "GN22881", "GN22883", "GN22884", "GN22885", "GN22887", "GN22888",
             "GN22889", "GN22891", "GN22892", "GN22893", "GN22894", "GN22895",
             "GN22897", "GN22898", "GN22900", "GN22902", "GN22904", "GN22905",
             "GN22906", "GN22907", "GN22909", "GN22910", "GN22913", "GN22917",
             "GN22918", "GN22920", "GN22921", "GN22926", "GN22927", "GN22928")

zebrashark.subset.percmeth.captive.df <- zebrashark.5x.percmeth.53S.impute.data.df[row.names(zebrashark.5x.percmeth.53S.impute.data.df) %in% captive,]
row.names(zebrashark.subset.percmeth.captive.df)
nrow(zebrashark.subset.percmeth.captive.df)
ncol(zebrashark.subset.percmeth.captive.df)

zebrashark.subset.percmeth.captive.df.t <- t(zebrashark.subset.percmeth.captive.df)

write.csv(zebrashark.subset.percmeth.captive.df.t, file="/dir/work/ZEBRA_SHARK/R/filteredByCov_ouput/zebrashark.5x.percmeth.53S.impute.captive.csv")

# test for and add age correlation information
# based on the old example, the input for the corAndPvalue function is a matrix of values, in this case perc methylation values, where each column is a locus, and each row is a sample
# the second input is the variable to test a correlation with, if temperature, then this a dataframe with one column of 24 temperature values
zebrashark_factors_cc<-read.csv(file="/dir/work/ZEBRA_SHARK/R/zebrashark_factors_cc.csv")

age<-zebrashark_factors_cc$est.age
age_d<-as.data.frame(age)
str(age_d)


# TEST FOR SPEARMAN CORRELATIONS
age.cp.spearman<-corAndPvalue(zebrashark.subset.percmeth.captive.df, age_d, method="spearman")
age.cp.spearman.pval<-age.cp.spearman$p
age.cp.spearman.pval.adj<-data.frame(p.adjust(age.cp.spearman$p, method="fdr"))
age.cp.spearman.cor<-age.cp.spearman$cor

age.cp.spearman.df<-cbind(age.cp.spearman.pval, age.cp.spearman.pval.adj,age.cp.spearman.cor)
colnames(age.cp.spearman.df)<-c("pval", "fdr", "cor")

write.csv(age.cp.spearman.df, file="/dir/work/ZEBRA_SHARK/R/age.cp.spearman.all.csv")

# TEST FOR PEARSON CORRELATIONS
age.cp.pearson<-corAndPvalue(zebrashark.subset.percmeth.captive.df, age_d, method="pearson")
age.cp.pearson.pval<-age.cp.pearson$p
age.cp.pearson.pval.adj<-data.frame(p.adjust(age.cp.pearson$p, method="fdr"))
age.cp.pearson.cor<-age.cp.pearson$cor

age.cp.pearson.df<-cbind(age.cp.pearson.pval, age.cp.pearson.pval.adj,age.cp.pearson.cor)
colnames(age.cp.pearson.df)<-c("pval", "fdr", "cor")

write.csv(age.cp.pearson.df, file="/dir/work/ZEBRA_SHARK/R/age.cp.pearson.all.csv")
