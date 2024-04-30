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

# subset dataframe to only include captive individuals
# wild individuals
wild.list <- c("GN22646", "GN22648", "GN22876", "GN22882", "GN22886", "GN22890", "GN22896", "GN22899", "GN22901", "GN22903", "GN22908", "GN22911", "GN22912", "GN22914", "GN22915", "GN22916", "GN22922", "GN22923")

# wild data only
zebrashark.subset.percmeth.wild.df <- zebrashark.subset.percmeth.df[,names(zebrashark.subset.percmeth.df) %in% wild.list]
head(zebrashark.subset.percmeth.wild.df)
nrow(zebrashark.subset.percmeth.wild.df)
ncol(zebrashark.subset.percmeth.wild.df)

# captive data only
zebrashark.subset.percmeth.captive.df <- zebrashark.subset.percmeth.df[,!names(zebrashark.subset.percmeth.df) %in% wild.list]
head(zebrashark.subset.percmeth.captive.df)
nrow(zebrashark.subset.percmeth.captive.df)
ncol(zebrashark.subset.percmeth.captive.df)

# read in metadata for captive individuals only
zebrashark_factors_cc<-read.csv(file="/scratch/sb61937/work/ZEBRA_SHARK/R/zebrashark_factors_cc.csv")

# create dataframe only including age data and GN identifiers
zebrashark_age_cc<-subset(zebrashark_factors_cc, select=c(gn.number, est.age))


# create list of samples - each one will be dropped in a given iteration of the leave-one-out cross-validation

uniqueSampleIDs<-c("GN22867", "GN22868", "GN22869", "GN22870", "GN22871", "GN22872")

for(i in 1:length(uniqueSampleIDs)){
  print(paste0("Processing sample: ", uniqueSampleIDs[i]))
  # remove selected sample from percent methylation matrix
  drop <- uniqueSampleIDs[i]
  # training data contains all but the sample to be left out
  meth_train <- zebrashark.subset.percmeth.captive.df[,!(names(zebrashark.subset.percmeth.captive.df) %in% drop)]
  # test data contains only the left out sample
  meth_test <- zebrashark.subset.percmeth.captive.df[,(names(zebrashark.subset.percmeth.captive.df) %in% drop)]

  # remove selected sample from age data
  age <- subset(zebrashark_age_cc, !(gn.number == drop))
  # include only selected sample in age data
  age_test <- subset(zebrashark_age_cc, (gn.number == drop))

  # create transposed matrix of methylation data
  meth_train_t <- as.matrix(t(meth_train))
  meth_test_t <- as.matrix(t(meth_test))

  model.cv<-cv.glmnet(meth_train_t, as.matrix(log10(age$est.age)),
                           alpha=0.5, nfolds=5, family="gaussian")
  best_lambda<-model.cv$lambda.min
  model<-glmnet(meth_train_t, as.matrix(log10(age$est.age)),
                     alpha=0.5, family="gaussian", nlambda=100)
  en_coeff<-coef(model, s=best_lambda)

  predict_age_train<-predict(model, meth_train_t, s=best_lambda)
  predict_age_test<-predict(model, meth_test_t, s=best_lambda)


  # create a data frame to hold results
  assign(paste('LOO_',uniqueSampleIDs[i],sep=''),meth_train)
  assign(paste('LOO_',uniqueSampleIDs[i],sep=''),meth_test)
  assign(paste('LOO_encoeff_',uniqueSampleIDs[i],sep=''),en_coeff)
  assign(paste('LOO_model_',uniqueSampleIDs[i],sep=''),model)
  assign(paste('LOO_modelcv_',uniqueSampleIDs[i],sep=''),model.cv)
  assign(paste('LOO_predicttest_',uniqueSampleIDs[i],sep=''),predict_age_test)
  assign(paste('LOO_predicttrain_',uniqueSampleIDs[i],sep=''),predict_age_train)

}

rm(zebrashark.subset.percmeth.df, zebrashark.subset.percmeth.captive.df)

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

GN22872_coeff_df<-data.frame(LOO_encoeff_GN22872[1:15578679,])
colnames(GN22872_coeff_df)<-"coeff"
nrow(subset(GN22872_coeff_df, abs(coeff) > 0))
GN22872_coeff_df_nz<-subset(GN22872_coeff_df, abs(coeff) > 0)


rm(list=setdiff(ls(), c("GN22867_coeff_df_nz", "GN22868_coeff_df_nz", "GN22869_coeff_df_nz",
                        "GN22870_coeff_df_nz", "GN22871_coeff_df_nz", "GN22872_coeff_df_nz",
                        "LOO_predicttrain_GN22867", "LOO_predicttrain_GN22868", "LOO_predicttrain_GN22869",
                        "LOO_predicttrain_GN22870", "LOO_predicttrain_GN22871", "LOO_predicttrain_GN22872",
                        "LOO_predicttest_GN22867", "LOO_predicttest_GN22868", "LOO_predicttest_GN22869",
                        "LOO_predicttest_GN22870", "LOO_predicttest_GN22871", "LOO_predicttest_GN22872")))

save(list = ls(all.names = TRUE), file = "/scratch/sb61937/work/ZEBRA_SHARK/R/log_glmnet_2_53S_short.RData")
