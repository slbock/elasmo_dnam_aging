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

########### EDIT BELOW ###########
# create list of samples - each one will be dropped in a given iteration of the leave-one-out cross-validation

oldc_uniqueSampleIDs<-c("GN22647", "GN22649", "GN22650", "GN22867", "GN22868", "GN22869",
                        "GN22870", "GN22871", "GN22873", "GN22874", "GN22875", "GN22883",
                        "GN22884", "GN22885", "GN22893", "GN22897", "GN22898", "GN22905",
                        "GN22906", "GN22907", "GN22909", "GN22910", "GN22913", "GN22917",
                        "GN22918", "GN22920", "GN22921", "GN22926")

for(i in 1:length(oldc_uniqueSampleIDs)){
  print(paste0("Processing sample: ", oldc_uniqueSampleIDs[i]))
  # remove selected sample from percent methylation matrix
  drop <- oldc_uniqueSampleIDs[i]
  # training data contains all but the sample to be left out
  meth_train <- zebrashark.subset.percmeth.old.c.df[,!(names(zebrashark.subset.percmeth.old.c.df) %in% drop)]
  # test data contains only the left out sample
  meth_test <- zebrashark.subset.percmeth.old.c.df[,(names(zebrashark.subset.percmeth.old.c.df) %in% drop)]

  # remove selected sample from age data
  age <- subset(zebrashark.old.c.age, !(gn.number == drop))
  # include only selected sample in age data
  age_test <- subset(zebrashark.old.c.age, (gn.number == drop))

  # create transposed matrix of methylation data
  meth_train_t <- as.matrix(t(meth_train))
  meth_test_t <- as.matrix(t(meth_test))

  model.cv<-cv.glmnet(meth_train_t, as.matrix(age$est.age),
                           alpha=0.5, nfolds=5, family="gaussian")
  best_lambda<-model.cv$lambda.min
  model<-glmnet(meth_train_t, as.matrix(age$est.age),
                     alpha=0.5, family="gaussian", nlambda=100)
  en_coeff<-coef(model, s=best_lambda)

  predict_age_train<-predict(model, meth_train_t, s=best_lambda)
  predict_age_test<-predict(model, meth_test_t, s=best_lambda)


  # create a data frame to hold results
  assign(paste('LOO_',oldc_uniqueSampleIDs[i],sep=''),meth_train)
  assign(paste('LOO_',oldc_uniqueSampleIDs[i],sep=''),meth_test)
  assign(paste('LOO_encoeff_',oldc_uniqueSampleIDs[i],sep=''),en_coeff)
  assign(paste('LOO_model_',oldc_uniqueSampleIDs[i],sep=''),model)
  assign(paste('LOO_modelcv_',oldc_uniqueSampleIDs[i],sep=''),model.cv)
  assign(paste('LOO_predicttest_',oldc_uniqueSampleIDs[i],sep=''),predict_age_test)
  assign(paste('LOO_predicttrain_',oldc_uniqueSampleIDs[i],sep=''),predict_age_train)

}

rm(zebrashark.subset.percmeth.df, zebrashark.subset.percmeth.old.c.df)

GN22647_coeff_df<-data.frame(LOO_encoeff_GN22647[1:15578679,])
colnames(GN22647_coeff_df)<-"coeff"
nrow(subset(GN22647_coeff_df, abs(coeff) > 0))
GN22647_coeff_df_nz<-subset(GN22647_coeff_df, abs(coeff) > 0)

GN22649_coeff_df<-data.frame(LOO_encoeff_GN22649[1:15578679,])
colnames(GN22649_coeff_df)<-"coeff"
nrow(subset(GN22649_coeff_df, abs(coeff) > 0))
GN22649_coeff_df_nz<-subset(GN22649_coeff_df, abs(coeff) > 0)

GN22650_coeff_df<-data.frame(LOO_encoeff_GN22650[1:15578679,])
colnames(GN22650_coeff_df)<-"coeff"
nrow(subset(GN22650_coeff_df, abs(coeff) > 0))
GN22650_coeff_df_nz<-subset(GN22650_coeff_df, abs(coeff) > 0)

GN22867_coeff_df<-data.frame(LOO_encoeff_GN22867[1:15578679,])
colnames(GN22867_coeff_df)<-"coeff"
nrow(subset(GN22867_coeff_df, abs(coeff) > 0))
GN22867_coeff_df_nz<-subset(GN22867_coeff_df, abs(coeff) > 0)

GN22868_coeff_df<-data.frame(LOO_encoeff_GN22868[1:15578679,])
colnames(GN22868_coeff_df)<-"coeff"
nrow(subset(GN22868_coeff_df, abs(coeff) > 0))
GN22868_coeff_df_nz<-subset(GN22868_coeff_df, abs(coeff) > 0)

GN22869_coeff_df<-data.frame(LOO_encoeff_GN22869[1:15578679,])
colnames(GN22869_coeff_df)<-"coeff"
nrow(subset(GN22869_coeff_df, abs(coeff) > 0))
GN22869_coeff_df_nz<-subset(GN22869_coeff_df, abs(coeff) > 0)

GN22870_coeff_df<-data.frame(LOO_encoeff_GN22870[1:15578679,])
colnames(GN22870_coeff_df)<-"coeff"
nrow(subset(GN22870_coeff_df, abs(coeff) > 0))
GN22870_coeff_df_nz<-subset(GN22870_coeff_df, abs(coeff) > 0)

GN22871_coeff_df<-data.frame(LOO_encoeff_GN22871[1:15578679,])
colnames(GN22871_coeff_df)<-"coeff"
nrow(subset(GN22871_coeff_df, abs(coeff) > 0))
GN22871_coeff_df_nz<-subset(GN22871_coeff_df, abs(coeff) > 0)

GN22873_coeff_df<-data.frame(LOO_encoeff_GN22873[1:15578679,])
colnames(GN22873_coeff_df)<-"coeff"
nrow(subset(GN22873_coeff_df, abs(coeff) > 0))
GN22873_coeff_df_nz<-subset(GN22873_coeff_df, abs(coeff) > 0)

GN22874_coeff_df<-data.frame(LOO_encoeff_GN22874[1:15578679,])
colnames(GN22874_coeff_df)<-"coeff"
nrow(subset(GN22874_coeff_df, abs(coeff) > 0))
GN22874_coeff_df_nz<-subset(GN22874_coeff_df, abs(coeff) > 0)

GN22875_coeff_df<-data.frame(LOO_encoeff_GN22875[1:15578679,])
colnames(GN22875_coeff_df)<-"coeff"
nrow(subset(GN22875_coeff_df, abs(coeff) > 0))
GN22875_coeff_df_nz<-subset(GN22875_coeff_df, abs(coeff) > 0)

GN22883_coeff_df<-data.frame(LOO_encoeff_GN22883[1:15578679,])
colnames(GN22883_coeff_df)<-"coeff"
nrow(subset(GN22883_coeff_df, abs(coeff) > 0))
GN22883_coeff_df_nz<-subset(GN22883_coeff_df, abs(coeff) > 0)

GN22884_coeff_df<-data.frame(LOO_encoeff_GN22884[1:15578679,])
colnames(GN22884_coeff_df)<-"coeff"
nrow(subset(GN22884_coeff_df, abs(coeff) > 0))
GN22884_coeff_df_nz<-subset(GN22884_coeff_df, abs(coeff) > 0)

GN22885_coeff_df<-data.frame(LOO_encoeff_GN22885[1:15578679,])
colnames(GN22885_coeff_df)<-"coeff"
nrow(subset(GN22885_coeff_df, abs(coeff) > 0))
GN22885_coeff_df_nz<-subset(GN22885_coeff_df, abs(coeff) > 0)

GN22893_coeff_df<-data.frame(LOO_encoeff_GN22893[1:15578679,])
colnames(GN22893_coeff_df)<-"coeff"
nrow(subset(GN22893_coeff_df, abs(coeff) > 0))
GN22893_coeff_df_nz<-subset(GN22893_coeff_df, abs(coeff) > 0)

GN22897_coeff_df<-data.frame(LOO_encoeff_GN22897[1:15578679,])
colnames(GN22897_coeff_df)<-"coeff"
nrow(subset(GN22897_coeff_df, abs(coeff) > 0))
GN22897_coeff_df_nz<-subset(GN22897_coeff_df, abs(coeff) > 0)

GN22898_coeff_df<-data.frame(LOO_encoeff_GN22898[1:15578679,])
colnames(GN22898_coeff_df)<-"coeff"
nrow(subset(GN22898_coeff_df, abs(coeff) > 0))
GN22898_coeff_df_nz<-subset(GN22898_coeff_df, abs(coeff) > 0)

GN22905_coeff_df<-data.frame(LOO_encoeff_GN22905[1:15578679,])
colnames(GN22905_coeff_df)<-"coeff"
nrow(subset(GN22905_coeff_df, abs(coeff) > 0))
GN22905_coeff_df_nz<-subset(GN22905_coeff_df, abs(coeff) > 0)

GN22906_coeff_df<-data.frame(LOO_encoeff_GN22906[1:15578679,])
colnames(GN22906_coeff_df)<-"coeff"
nrow(subset(GN22906_coeff_df, abs(coeff) > 0))
GN22906_coeff_df_nz<-subset(GN22906_coeff_df, abs(coeff) > 0)

GN22907_coeff_df<-data.frame(LOO_encoeff_GN22907[1:15578679,])
colnames(GN22907_coeff_df)<-"coeff"
nrow(subset(GN22907_coeff_df, abs(coeff) > 0))
GN22907_coeff_df_nz<-subset(GN22907_coeff_df, abs(coeff) > 0)

GN22909_coeff_df<-data.frame(LOO_encoeff_GN22909[1:15578679,])
colnames(GN22909_coeff_df)<-"coeff"
nrow(subset(GN22909_coeff_df, abs(coeff) > 0))
GN22909_coeff_df_nz<-subset(GN22909_coeff_df, abs(coeff) > 0)

GN22910_coeff_df<-data.frame(LOO_encoeff_GN22910[1:15578679,])
colnames(GN22910_coeff_df)<-"coeff"
nrow(subset(GN22910_coeff_df, abs(coeff) > 0))
GN22910_coeff_df_nz<-subset(GN22910_coeff_df, abs(coeff) > 0)

GN22913_coeff_df<-data.frame(LOO_encoeff_GN22913[1:15578679,])
colnames(GN22913_coeff_df)<-"coeff"
nrow(subset(GN22913_coeff_df, abs(coeff) > 0))
GN22913_coeff_df_nz<-subset(GN22913_coeff_df, abs(coeff) > 0)

GN22917_coeff_df<-data.frame(LOO_encoeff_GN22917[1:15578679,])
colnames(GN22917_coeff_df)<-"coeff"
nrow(subset(GN22917_coeff_df, abs(coeff) > 0))
GN22917_coeff_df_nz<-subset(GN22917_coeff_df, abs(coeff) > 0)

GN22918_coeff_df<-data.frame(LOO_encoeff_GN22918[1:15578679,])
colnames(GN22918_coeff_df)<-"coeff"
nrow(subset(GN22918_coeff_df, abs(coeff) > 0))
GN22918_coeff_df_nz<-subset(GN22918_coeff_df, abs(coeff) > 0)

GN22920_coeff_df<-data.frame(LOO_encoeff_GN22920[1:15578679,])
colnames(GN22920_coeff_df)<-"coeff"
nrow(subset(GN22920_coeff_df, abs(coeff) > 0))
GN22920_coeff_df_nz<-subset(GN22920_coeff_df, abs(coeff) > 0)

GN22921_coeff_df<-data.frame(LOO_encoeff_GN22921[1:15578679,])
colnames(GN22921_coeff_df)<-"coeff"
nrow(subset(GN22921_coeff_df, abs(coeff) > 0))
GN22921_coeff_df_nz<-subset(GN22921_coeff_df, abs(coeff) > 0)

GN22926_coeff_df<-data.frame(LOO_encoeff_GN22926[1:15578679,])
colnames(GN22926_coeff_df)<-"coeff"
nrow(subset(GN22926_coeff_df, abs(coeff) > 0))
GN22926_coeff_df_nz<-subset(GN22926_coeff_df, abs(coeff) > 0)


rm(list=setdiff(ls(), c("GN22647_coeff_df", "GN22649_coeff_df", "GN22650_coeff_df", "GN22867_coeff_df", "GN22868_coeff_df", "GN22869_coeff_df",
                         "GN22870_coeff_df", "GN22871_coeff_df", "GN22873_coeff_df", "GN22874_coeff_df", "GN22875_coeff_df", "GN22883_coeff_df",
                         "GN22884_coeff_df", "GN22885_coeff_df", "GN22893_coeff_df", "GN22897_coeff_df", "GN22898_coeff_df", "GN22905_coeff_df",
                         "GN22906_coeff_df", "GN22907_coeff_df", "GN22909_coeff_df", "GN22910_coeff_df", "GN22913_coeff_df", "GN22917_coeff_df",
                         "GN22918_coeff_df", "GN22920_coeff_df", "GN22921_coeff_df", "GN22926_coeff_df",
                         "LOO_predicttest_GN22647", "LOO_predicttest_GN22649", "LOO_predicttest_GN22650", "LOO_predicttest_GN22867", "LOO_predicttest_GN22868", "LOO_predicttest_GN22869",
                         "LOO_predicttest_GN22870", "LOO_predicttest_GN22871", "LOO_predicttest_GN22873", "LOO_predicttest_GN22874", "LOO_predicttest_GN22875", "LOO_predicttest_GN22883",
                         "LOO_predicttest_GN22884", "LOO_predicttest_GN22885", "LOO_predicttest_GN22893", "LOO_predicttest_GN22897", "LOO_predicttest_GN22898", "LOO_predicttest_GN22905",
                         "LOO_predicttest_GN22906", "LOO_predicttest_GN22907", "LOO_predicttest_GN22909", "LOO_predicttest_GN22910", "LOO_predicttest_GN22913", "LOO_predicttest_GN22917",
                         "LOO_predicttest_GN22918", "LOO_predicttest_GN22920", "LOO_predicttest_GN22921", "LOO_predicttest_GN22926",
                         "LOO_predicttrain_GN22647", "LOO_predicttrain_GN22649", "LOO_predicttrain_GN22650", "LOO_predicttrain_GN22867", "LOO_predicttrain_GN22868", "LOO_predicttrain_GN22869",
                         "LOO_predicttrain_GN22870", "LOO_predicttrain_GN22871", "LOO_predicttrain_GN22873", "LOO_predicttrain_GN22874", "LOO_predicttrain_GN22875", "LOO_predicttrain_GN22883",
                         "LOO_predicttrain_GN22884", "LOO_predicttrain_GN22885", "LOO_predicttrain_GN22893", "LOO_predicttrain_GN22897", "LOO_predicttrain_GN22898", "LOO_predicttrain_GN22905",
                         "LOO_predicttrain_GN22906", "LOO_predicttrain_GN22907", "LOO_predicttrain_GN22909", "LOO_predicttrain_GN22910", "LOO_predicttrain_GN22913", "LOO_predicttrain_GN22917",
                         "LOO_predicttrain_GN22918", "LOO_predicttrain_GN22920", "LOO_predicttrain_GN22921", "LOO_predicttrain_GN22926")))

save(list = ls(all.names = TRUE), file = "/scratch/sb61937/work/ZEBRA_SHARK/R/old_glmnet_loocv.RData")
