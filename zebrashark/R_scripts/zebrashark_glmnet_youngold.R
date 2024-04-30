library(impute)
library(dplyr)
library(caret)
library(methylKit)
library(GenomicRanges)
library(genomation)
library(WGCNA)
library(reshape2)
library(glmnet)

load("/scratch/sb61937/work/ZEBRA_SHARK/R/filteredByCov_ouput/zebrashark.5x.percmeth.53S.impute.RData")

zebrashark.5x.percmeth.53S.impute.data <- zebrashark.5x.percmeth.53S.impute$data
zebrashark.5x.percmeth.53S.impute.data.df <- data.frame(zebrashark.5x.percmeth.53S.impute.data)
rownames(zebrashark.5x.percmeth.53S.impute.data.df)

## Set this option so genomic coordinates do not get printed as scientific notation
options(scipen=50)

# load percent methylation data for covered, filtered loci (n = 15,578,679 CpGs)
zebrashark.subset.percmeth.df<-data.frame(t(zebrashark.5x.percmeth.53S.impute.data.df))
head(zebrashark.subset.percmeth.df)

rm(zebrashark.5x.percmeth.53S.impute, zebrashark.5x.percmeth.53S.impute.data, zebrashark.5x.percmeth.53S.impute.data.df)

zebrashark_factors <- read.csv("/scratch/sb61937/work/ZEBRA_SHARK/R/zebrashark_factors.csv")

#### divide dataframe into pre- and post-maturity ('young' and 'old' respectively) ####
## YOUNG ##
young.list <- c("GN22651", "GN22652", "GN22872", "GN22877", "GN22878", "GN22879",
                "GN22880", "GN22881", "GN22887", "GN22888", "GN22889", "GN22891",
                "GN22892", "GN22894", "GN22895", "GN22900", "GN22902", "GN22904",
                "GN22927", "GN22928")
# young data only
zebrashark.subset.percmeth.young.df <- zebrashark.subset.percmeth.df[,names(zebrashark.subset.percmeth.df) %in% young.list]
head(zebrashark.subset.percmeth.young.df)
nrow(zebrashark.subset.percmeth.young.df)
ncol(zebrashark.subset.percmeth.young.df)

# get metadata for subset
zebrashark.young.factors <- zebrashark_factors[zebrashark_factors$gn.number %in% young.list, ]
# create dataframe only including age data and GN identifiers
zebrashark.young.age<-subset(zebrashark.young.factors, select=c(gn.number, est.age))

## OLD - CAPTIVE ##
old.c.list <- c("GN22647", "GN22649", "GN22650", "GN22867", "GN22868", "GN22869",
                "GN22870", "GN22871", "GN22873", "GN22874", "GN22875", "GN22883",
                "GN22884", "GN22885", "GN22893", "GN22897", "GN22898", "GN22905",
                "GN22906", "GN22907", "GN22909", "GN22910", "GN22913", "GN22917",
                "GN22918", "GN22920", "GN22921", "GN22926")
# old.c data only
zebrashark.subset.percmeth.old.c.df <- zebrashark.subset.percmeth.df[,names(zebrashark.subset.percmeth.df) %in% old.c.list]
head(zebrashark.subset.percmeth.old.c.df)
nrow(zebrashark.subset.percmeth.old.c.df)
ncol(zebrashark.subset.percmeth.old.c.df)

# get metadata for subset
zebrashark.old.c.factors <- zebrashark_factors[zebrashark_factors$gn.number %in% old.c.list, ]
# create dataframe only including age data and GN identifiers
zebrashark.old.c.age<-subset(zebrashark.old.c.factors, select=c(gn.number, est.age))

## OLD - WILD ##
old.w.list <- c("GN22646", "GN22648", "GN22876", "GN22882", "GN22886", "GN22890",
                "GN22896", "GN22899", "GN22901", "GN22903", "GN22908", "GN22911",
                "GN22912", "GN22914", "GN22915", "GN22916", "GN22922", "GN22923")
# old.w data only
zebrashark.subset.percmeth.old.w.df <- zebrashark.subset.percmeth.df[,names(zebrashark.subset.percmeth.df) %in% old.w.list]
head(zebrashark.subset.percmeth.old.w.df)
nrow(zebrashark.subset.percmeth.old.w.df)
ncol(zebrashark.subset.percmeth.old.w.df)

# get metadata for subset
zebrashark.old.w.factors <- zebrashark_factors[zebrashark_factors$gn.number %in% old.w.list, ]
# create dataframe only including age data and GN identifiers
zebrashark.old.w.age<-subset(zebrashark.old.w.factors, select=c(gn.number, est.age))

#### single training set - YOUNG ####
# create transposed matrix of methylation data
meth_young_t <- as.matrix(t(zebrashark.subset.percmeth.young.df))
meth_oldc_t <- as.matrix(t(zebrashark.subset.percmeth.old.c.df))
meth_oldw_t <- as.matrix(t(zebrashark.subset.percmeth.old.w.df))

young.model.cv<-cv.glmnet(meth_young_t, as.matrix(zebrashark.young.age$est.age),
                           alpha=0.5, nfolds=5, family="gaussian")
young_best_lambda<-young.model.cv$lambda.min
young.model<-glmnet(meth_young_t, as.matrix(zebrashark.young.age$est.age),
                     alpha=0.5, family="gaussian", nlambda=100)
young_en_coeff<-coef(young.model, s=young_best_lambda)

young_predict_age_young<-predict(young.model, meth_young_t, s=young_best_lambda)
young_predict_age_oldc<-predict(young.model, meth_oldc_t, s=young_best_lambda)
young_predict_age_oldw<-predict(young.model, meth_oldw_t, s=young_best_lambda)


#### single training set - OLD-C ####

oldc.model.cv<-cv.glmnet(meth_oldc_t, as.matrix(zebrashark.old.c.age$est.age),
                           alpha=0.5, nfolds=5, family="gaussian")
oldc_best_lambda<-oldc.model.cv$lambda.min
oldc.model<-glmnet(meth_oldc_t, as.matrix(zebrashark.old.c.age$est.age),
                     alpha=0.5, family="gaussian", nlambda=100)
oldc_en_coeff<-coef(oldc.model, s=oldc_best_lambda)

oldc_predict_age_young<-predict(oldc.model, meth_young_t, s=oldc_best_lambda)
oldc_predict_age_oldc<-predict(oldc.model, meth_oldc_t, s=oldc_best_lambda)
oldc_predict_age_oldw<-predict(oldc.model, meth_oldw_t, s=oldc_best_lambda)

rm(list=setdiff(ls(), c("young_en_coeff", "young_predict_age_young", "young_predict_age_oldc", "young_predict_age_oldw", "young_best_lambda",
                        "oldc_en_coeff", "oldc_predict_age_young", "oldc_predict_age_oldc", "oldc_predict_age_oldw", "oldc_best_lambda",
                        "zebrashark.old.c.age", "zebrashark.old.w.age", zebrashark.young.age)))

save(list = ls(all.names = TRUE), file = "/scratch/sb61937/work/ZEBRA_SHARK/R/glmnet_young_old.RData")
