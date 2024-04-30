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

# create separate datasets for training set and test set
zebrashark_test <- c("GN22875", "GN22877", "GN22881", "GN22884", "GN22893", "GN22897", "GN22902", "GN22906")
zebrashark_wild <- c("GN22646", "GN22648", "GN22876", "GN22882", "GN22886", "GN22890", "GN22896", "GN22899",
                     "GN22901", "GN22903", "GN22908", "GN22911", "GN22912", "GN22914", "GN22915", "GN22916",
                     "GN22922", "GN22923")

zebrashark.percmeth.captive.test.df <- zebrashark.subset.percmeth.df[,names(zebrashark.subset.percmeth.df) %in% zebrashark_test]
zebrashark.percmeth.all.training.df <- zebrashark.subset.percmeth.df[,!names(zebrashark.subset.percmeth.df) %in% zebrashark_test]
zebrashark.percmeth.captive.training.df <- zebrashark.percmeth.all.training.df[,!names(zebrashark.percmeth.all.training.df) %in% zebrashark_wild]

head(zebrashark.percmeth.captive.test.df)
nrow(zebrashark.percmeth.captive.test.df)
ncol(zebrashark.percmeth.captive.test.df)

head(zebrashark.percmeth.captive.training.df)
nrow(zebrashark.percmeth.captive.training.df)
ncol(zebrashark.percmeth.captive.training.df)

# read in metadata for captive individuals only
zebrashark_factors_cc<-read.csv(file="/scratch/sb61937/work/ZEBRA_SHARK/R/zebrashark_factors_cc.csv")

zebrashark_factors_test<- subset(zebrashark_factors_cc, gn.number %in% zebrashark_test)
zebrashark_factors_training<- subset(zebrashark_factors_cc, !(gn.number %in% zebrashark_test))

# create dataframe only including age data and GN identifiers
zebrashark_age_test<-subset(zebrashark_factors_test, select=c(gn.number, est.age))
zebrashark_age_training<-subset(zebrashark_factors_training, select=c(gn.number, est.age))



# create transposed matrix of methylation data
meth_train_t <- as.matrix(t(zebrashark.percmeth.captive.training.df))
meth_test_t <- as.matrix(t(zebrashark.percmeth.captive.test.df))

model.cv<-cv.glmnet(meth_train_t, as.matrix(zebrashark_age_training$est.age),
                           alpha=0.5, nfolds=5, family="gaussian")
best_lambda<-model.cv$lambda.min
model<-glmnet(meth_train_t, as.matrix(zebrashark_age_training$est.age),
                     alpha=0.5, family="gaussian", nlambda=100)
en_coeff<-coef(model, s=best_lambda)

predict_age_train<-predict(model, meth_train_t, s=best_lambda)
predict_age_test<-predict(model, meth_test_t, s=best_lambda)


rm(list=setdiff(ls(), c("en_coeff", "predict_age_train", "predict_age_test", "best_lambda")))

save(list = ls(all.names = TRUE), file = "/scratch/sb61937/work/ZEBRA_SHARK/R/glmnet_single_train_test.RData")
