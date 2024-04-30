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

########### EDIT BELOW ###########
# create list of samples - each one will be dropped in a given iteration of the leave-one-out cross-validation

young_uniqueSampleIDs<-c("GN22651", "GN22652", "GN22872", "GN22877", "GN22878", "GN22879",
                         "GN22880", "GN22881", "GN22887", "GN22888", "GN22889", "GN22891",
                         "GN22892", "GN22894", "GN22895", "GN22900", "GN22902", "GN22904",
                         "GN22927", "GN22928")

for(i in 1:length(young_uniqueSampleIDs)){
  print(paste0("Processing sample: ", young_uniqueSampleIDs[i]))
  # remove selected sample from percent methylation matrix
  drop <- young_uniqueSampleIDs[i]
  # training data contains all but the sample to be left out
  meth_train <- zebrashark.subset.percmeth.young.df[,!(names(zebrashark.subset.percmeth.young.df) %in% drop)]
  # test data contains only the left out sample
  meth_test <- zebrashark.subset.percmeth.young.df[,(names(zebrashark.subset.percmeth.young.df) %in% drop)]

  # remove selected sample from age data
  age <- subset(zebrashark.young.age, !(gn.number == drop))
  # include only selected sample in age data
  age_test <- subset(zebrashark.young.age, (gn.number == drop))

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
  assign(paste('LOO_',young_uniqueSampleIDs[i],sep=''),meth_train)
  assign(paste('LOO_',young_uniqueSampleIDs[i],sep=''),meth_test)
  assign(paste('LOO_encoeff_',young_uniqueSampleIDs[i],sep=''),en_coeff)
  assign(paste('LOO_model_',young_uniqueSampleIDs[i],sep=''),model)
  assign(paste('LOO_modelcv_',young_uniqueSampleIDs[i],sep=''),model.cv)
  assign(paste('LOO_predicttest_',young_uniqueSampleIDs[i],sep=''),predict_age_test)
  assign(paste('LOO_predicttrain_',young_uniqueSampleIDs[i],sep=''),predict_age_train)

}

rm(zebrashark.subset.percmeth.df, zebrashark.subset.percmeth.young.df)

GN22651_coeff_df<-data.frame(LOO_encoeff_GN22651[1:15578679,])
colnames(GN22651_coeff_df)<-"coeff"
nrow(subset(GN22651_coeff_df, abs(coeff) > 0))
GN22651_coeff_df_nz<-subset(GN22651_coeff_df, abs(coeff) > 0)

GN22652_coeff_df<-data.frame(LOO_encoeff_GN22652[1:15578679,])
colnames(GN22652_coeff_df)<-"coeff"
nrow(subset(GN22652_coeff_df, abs(coeff) > 0))
GN22652_coeff_df_nz<-subset(GN22652_coeff_df, abs(coeff) > 0)

GN22872_coeff_df<-data.frame(LOO_encoeff_GN22872[1:15578679,])
colnames(GN22872_coeff_df)<-"coeff"
nrow(subset(GN22872_coeff_df, abs(coeff) > 0))
GN22872_coeff_df_nz<-subset(GN22872_coeff_df, abs(coeff) > 0)

GN22877_coeff_df<-data.frame(LOO_encoeff_GN22877[1:15578679,])
colnames(GN22877_coeff_df)<-"coeff"
nrow(subset(GN22877_coeff_df, abs(coeff) > 0))
GN22877_coeff_df_nz<-subset(GN22877_coeff_df, abs(coeff) > 0)

GN22878_coeff_df<-data.frame(LOO_encoeff_GN22878[1:15578679,])
colnames(GN22878_coeff_df)<-"coeff"
nrow(subset(GN22878_coeff_df, abs(coeff) > 0))
GN22878_coeff_df_nz<-subset(GN22878_coeff_df, abs(coeff) > 0)

GN22879_coeff_df<-data.frame(LOO_encoeff_GN22879[1:15578679,])
colnames(GN22879_coeff_df)<-"coeff"
nrow(subset(GN22879_coeff_df, abs(coeff) > 0))
GN22879_coeff_df_nz<-subset(GN22879_coeff_df, abs(coeff) > 0)

GN22880_coeff_df<-data.frame(LOO_encoeff_GN22880[1:15578679,])
colnames(GN22880_coeff_df)<-"coeff"
nrow(subset(GN22880_coeff_df, abs(coeff) > 0))
GN22880_coeff_df_nz<-subset(GN22880_coeff_df, abs(coeff) > 0)

GN22881_coeff_df<-data.frame(LOO_encoeff_GN22881[1:15578679,])
colnames(GN22881_coeff_df)<-"coeff"
nrow(subset(GN22881_coeff_df, abs(coeff) > 0))
GN22881_coeff_df_nz<-subset(GN22881_coeff_df, abs(coeff) > 0)

GN22887_coeff_df<-data.frame(LOO_encoeff_GN22887[1:15578679,])
colnames(GN22887_coeff_df)<-"coeff"
nrow(subset(GN22887_coeff_df, abs(coeff) > 0))
GN22887_coeff_df_nz<-subset(GN22887_coeff_df, abs(coeff) > 0)

GN22888_coeff_df<-data.frame(LOO_encoeff_GN22888[1:15578679,])
colnames(GN22888_coeff_df)<-"coeff"
nrow(subset(GN22888_coeff_df, abs(coeff) > 0))
GN22888_coeff_df_nz<-subset(GN22888_coeff_df, abs(coeff) > 0)

GN22889_coeff_df<-data.frame(LOO_encoeff_GN22889[1:15578679,])
colnames(GN22889_coeff_df)<-"coeff"
nrow(subset(GN22889_coeff_df, abs(coeff) > 0))
GN22889_coeff_df_nz<-subset(GN22889_coeff_df, abs(coeff) > 0)

GN22891_coeff_df<-data.frame(LOO_encoeff_GN22891[1:15578679,])
colnames(GN22891_coeff_df)<-"coeff"
nrow(subset(GN22891_coeff_df, abs(coeff) > 0))
GN22891_coeff_df_nz<-subset(GN22891_coeff_df, abs(coeff) > 0)

GN22892_coeff_df<-data.frame(LOO_encoeff_GN22892[1:15578679,])
colnames(GN22892_coeff_df)<-"coeff"
nrow(subset(GN22892_coeff_df, abs(coeff) > 0))
GN22892_coeff_df_nz<-subset(GN22892_coeff_df, abs(coeff) > 0)

GN22894_coeff_df<-data.frame(LOO_encoeff_GN22894[1:15578679,])
colnames(GN22894_coeff_df)<-"coeff"
nrow(subset(GN22894_coeff_df, abs(coeff) > 0))
GN22894_coeff_df_nz<-subset(GN22894_coeff_df, abs(coeff) > 0)

GN22895_coeff_df<-data.frame(LOO_encoeff_GN22895[1:15578679,])
colnames(GN22895_coeff_df)<-"coeff"
nrow(subset(GN22895_coeff_df, abs(coeff) > 0))
GN22895_coeff_df_nz<-subset(GN22895_coeff_df, abs(coeff) > 0)

GN22900_coeff_df<-data.frame(LOO_encoeff_GN22900[1:15578679,])
colnames(GN22900_coeff_df)<-"coeff"
nrow(subset(GN22900_coeff_df, abs(coeff) > 0))
GN22900_coeff_df_nz<-subset(GN22900_coeff_df, abs(coeff) > 0)

GN22902_coeff_df<-data.frame(LOO_encoeff_GN22902[1:15578679,])
colnames(GN22902_coeff_df)<-"coeff"
nrow(subset(GN22902_coeff_df, abs(coeff) > 0))
GN22902_coeff_df_nz<-subset(GN22902_coeff_df, abs(coeff) > 0)

GN22904_coeff_df<-data.frame(LOO_encoeff_GN22904[1:15578679,])
colnames(GN22904_coeff_df)<-"coeff"
nrow(subset(GN22904_coeff_df, abs(coeff) > 0))
GN22904_coeff_df_nz<-subset(GN22904_coeff_df, abs(coeff) > 0)

GN22927_coeff_df<-data.frame(LOO_encoeff_GN22927[1:15578679,])
colnames(GN22927_coeff_df)<-"coeff"
nrow(subset(GN22927_coeff_df, abs(coeff) > 0))
GN22927_coeff_df_nz<-subset(GN22927_coeff_df, abs(coeff) > 0)

GN22928_coeff_df<-data.frame(LOO_encoeff_GN22928[1:15578679,])
colnames(GN22928_coeff_df)<-"coeff"
nrow(subset(GN22928_coeff_df, abs(coeff) > 0))
GN22928_coeff_df_nz<-subset(GN22928_coeff_df, abs(coeff) > 0)


rm(list=setdiff(ls(), c("GN22651_coeff_df", "GN22652_coeff_df", "GN22872_coeff_df", "GN22877_coeff_df", "GN22878_coeff_df", "GN22879_coeff_df",
                         "GN22880_coeff_df", "GN22881_coeff_df", "GN22887_coeff_df", "GN22888_coeff_df", "GN22889_coeff_df", "GN22891_coeff_df",
                         "GN22892_coeff_df", "GN22894_coeff_df", "GN22895_coeff_df", "GN22900_coeff_df", "GN22902_coeff_df", "GN22904_coeff_df",
                         "GN22927_coeff_df", "GN22928_coeff_df",
                         "LOO_predicttest_GN22651", "LOO_predicttest_GN22652", "LOO_predicttest_GN22872", "LOO_predicttest_GN22877", "LOO_predicttest_GN22878", "LOO_predicttest_GN22879",
                         "LOO_predicttest_GN22880", "LOO_predicttest_GN22881", "LOO_predicttest_GN22887", "LOO_predicttest_GN22888", "LOO_predicttest_GN22889", "LOO_predicttest_GN22891",
                         "LOO_predicttest_GN22892", "LOO_predicttest_GN22894", "LOO_predicttest_GN22895", "LOO_predicttest_GN22900", "LOO_predicttest_GN22902", "LOO_predicttest_GN22904",
                         "LOO_predicttest_GN22927", "LOO_predicttest_GN22928",
                         "LOO_predicttrain_GN22651", "LOO_predicttrain_GN22652", "LOO_predicttrain_GN22872", "LOO_predicttrain_GN22877", "LOO_predicttrain_GN22878", "LOO_predicttrain_GN22879",
                         "LOO_predicttrain_GN22880", "LOO_predicttrain_GN22881", "LOO_predicttrain_GN22887", "LOO_predicttrain_GN22888", "LOO_predicttrain_GN22889", "LOO_predicttrain_GN22891",
                         "LOO_predicttrain_GN22892", "LOO_predicttrain_GN22894", "LOO_predicttrain_GN22895", "LOO_predicttrain_GN22900", "LOO_predicttrain_GN22902", "LOO_predicttrain_GN22904",
                         "LOO_predicttrain_GN22927", "LOO_predicttrain_GN22928")))

save(list = ls(all.names = TRUE), file = "/scratch/sb61937/work/ZEBRA_SHARK/R/young_glmnet_loocv.RData")
