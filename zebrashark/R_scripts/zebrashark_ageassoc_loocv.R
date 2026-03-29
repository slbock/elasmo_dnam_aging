library(impute)
library(dplyr)
library(caret)
library(methylKit)
library(GenomicRanges)
library(genomation)
library(WGCNA)
library(reshape2)
library(glmnet)

## Set this option so genomic coordinates do not get printed as scientific notation
options(scipen=50)

# load percent methylation data for covered, filtered loci (n = 950,534 CpGs)
zebrashark.subset.percmeth.df<-read.csv(file="/dir/work/ZEBRA_SHARK/R/zebrashark.5x.percmeth.removinvar.csv")
head(zebrashark.subset.percmeth.df)
rownames(zebrashark.subset.percmeth.df)<-zebrashark.subset.percmeth.df$X
zebrashark.subset.percmeth.df$loc<-zebrashark.subset.percmeth.df$X


# subset dataframe to only include captive individuals
# wild individuals to be removed: GN22646, GN22648, GN22876, GN22882, GN22886, GN22890
zebrashark.subset.percmeth.captive.df<-subset(zebrashark.subset.percmeth.df, select = -c(X, GN22646, GN22648, GN22876, GN22882, GN22886, GN22890))
head(zebrashark.subset.percmeth.captive.df)
nrow(zebrashark.subset.percmeth.captive.df)
ncol(zebrashark.subset.percmeth.captive.df)

# read in metadata for captive individuals only
zebrashark_factors_cc<-read.csv(file="/dir/work/ZEBRA_SHARK/R/zebrashark_factors_cc.csv")

# create dataframe only including age data and GN identifiers
zebrashark_age_cc<-subset(zebrashark_factors_cc, select=c(gn.number, est.age))


# create list of samples - each one will be dropped in a given iteration of the leave-one-out cross-validation

uniqueSampleIDs<-c("GN22647", "GN22649", "GN22650", "GN22651", "GN22652",
                   "GN22867", "GN22868", "GN22869", "GN22870", "GN22871",
                   "GN22872", "GN22873", "GN22874", "GN22875", "GN22877",
                   "GN22878", "GN22879", "GN22880", "GN22881", "GN22883",
                   "GN22884", "GN22885", "GN22887", "GN22888", "GN22889",
                   "GN22891", "GN22892", "GN22893", "GN22894")

for(i in 1:length(uniqueSampleIDs)){
  print(paste0("Processing sample: ", uniqueSampleIDs[i]))
  # remove selected sample from percent methylation matrix
  drop <- uniqueSampleIDs[i]
  # training data contains all but the sample to be left out
  meth_train_loc <- zebrashark.subset.percmeth.captive.df[,!(names(zebrashark.subset.percmeth.captive.df) %in% drop)]
  # test data contains only the left out sample
  meth_test_loc_1 <- zebrashark.subset.percmeth.captive.df[,(names(zebrashark.subset.percmeth.captive.df) %in% drop)]
  meth_test_loc<-data.frame(meth=meth_test_loc_1, loc=zebrashark.subset.percmeth.captive.df$loc)

  # remove selected sample from age data
  age <- subset(zebrashark_age_cc, !(gn.number == drop))
  # include only selected sample in age data
  age_test <- subset(zebrashark_age_cc, (gn.number == drop))

  # create transposed matrix of methylation data
  meth_train_loc_t <- as.matrix(t(meth_train_loc[,1:28]))
  meth_test_loc_t <- as.matrix(t(meth_test_loc[,1]))

  # test for spearman correlations
  pearson<-corAndPvalue(meth_train_loc_t, age$est.age, method="pearson")

  # add spearman correlations to original methylation dataset
  meth_train_loc<-cbind(meth_train_loc, cor=pearson$cor)
  meth_test_loc<-cbind(meth_test_loc, cor=pearson$cor)

  # subset methylation data to only include age-associated sites
  meth_train_loc_ageassoc<-subset(meth_train_loc, abs(cor) > 0.5)
  meth_test_loc_ageassoc<-subset(meth_test_loc, abs(cor) > 0.5)

  # create transposed methylation matrix with age-associated loci only
  meth_train_loc_ageassoc_t<-as.matrix(t(meth_train_loc_ageassoc[,1:28]))
  meth_test_loc_ageassoc_t<-as.matrix(t(meth_test_loc_ageassoc[,1]))

  model.cv<-cv.glmnet(meth_train_loc_ageassoc_t, as.matrix(age$est.age),
                           alpha=0.5, nfolds=5, family="gaussian")
  best_lambda<-model.cv$lambda.min
  model<-glmnet(meth_train_loc_ageassoc_t, as.matrix(age$est.age),
                     alpha=0.5, family="gaussian", nlambda=100)
  en_coeff<-coef(model, s=best_lambda)

  predict_age_train<-predict(model, meth_train_loc_ageassoc_t, s=best_lambda)
  predict_age_test<-predict(model, meth_test_loc_ageassoc_t, s=best_lambda)


  # create a data frame to hold results
  assign(paste('LOO_',uniqueSampleIDs[i],sep=''),meth_train_loc)
  assign(paste('LOO_test_',uniqueSampleIDs[i],sep=''),meth_test_loc)
  assign(paste('LOO_encoeff_',uniqueSampleIDs[i],sep=''),en_coeff)
  assign(paste('LOO_model_',uniqueSampleIDs[i],sep=''),model)
  assign(paste('LOO_modelcv_',uniqueSampleIDs[i],sep=''),model.cv)
  assign(paste('LOO_predicttest_',uniqueSampleIDs[i],sep=''),predict_age_test)
  assign(paste('LOO_predicttrain_',uniqueSampleIDs[i],sep=''),predict_age_train)

}

warnings()
rm(list=setdiff(ls(), c("LOO_GN22647", "LOO_test_GN22647", "LOO_encoeff_GN22647", "LOO_model_GN22647", "LOO_modelcv_GN22647", "LOO_predicttest_GN22647", "LOO_predicttrain_GN22647",
                        "LOO_GN22649", "LOO_test_GN22649", "LOO_encoeff_GN22649", "LOO_model_GN22649", "LOO_modelcv_GN22649", "LOO_predicttest_GN22649", "LOO_predicttrain_GN22649",
                        "LOO_GN22650", "LOO_test_GN22650", "LOO_encoeff_GN22650", "LOO_model_GN22650", "LOO_modelcv_GN22650", "LOO_predicttest_GN22650", "LOO_predicttrain_GN22650",
                        "LOO_GN22651", "LOO_test_GN22651", "LOO_encoeff_GN22651", "LOO_model_GN22651", "LOO_modelcv_GN22651", "LOO_predicttest_GN22651", "LOO_predicttrain_GN22651",
                        "LOO_GN22652", "LOO_test_GN22652", "LOO_encoeff_GN22652", "LOO_model_GN22652", "LOO_modelcv_GN22652", "LOO_predicttest_GN22652", "LOO_predicttrain_GN22652",
                        "LOO_GN22867", "LOO_test_GN22867", "LOO_encoeff_GN22867", "LOO_model_GN22867", "LOO_modelcv_GN22867", "LOO_predicttest_GN22867", "LOO_predicttrain_GN22867",
                        "LOO_GN22868", "LOO_test_GN22868", "LOO_encoeff_GN22868", "LOO_model_GN22868", "LOO_modelcv_GN22868", "LOO_predicttest_GN22868", "LOO_predicttrain_GN22868",
                        "LOO_GN22869", "LOO_test_GN22869", "LOO_encoeff_GN22869", "LOO_model_GN22869", "LOO_modelcv_GN22869", "LOO_predicttest_GN22869", "LOO_predicttrain_GN22869",
                        "LOO_GN22870", "LOO_test_GN22870", "LOO_encoeff_GN22870", "LOO_model_GN22870", "LOO_modelcv_GN22870", "LOO_predicttest_GN22870", "LOO_predicttrain_GN22870",
                        "LOO_GN22871", "LOO_test_GN22871", "LOO_encoeff_GN22871", "LOO_model_GN22871", "LOO_modelcv_GN22871", "LOO_predicttest_GN22871", "LOO_predicttrain_GN22871",
                        "LOO_GN22872", "LOO_test_GN22872", "LOO_encoeff_GN22872", "LOO_model_GN22872", "LOO_modelcv_GN22872", "LOO_predicttest_GN22872", "LOO_predicttrain_GN22872",
                        "LOO_GN22873", "LOO_test_GN22873", "LOO_encoeff_GN22873", "LOO_model_GN22873", "LOO_modelcv_GN22873", "LOO_predicttest_GN22873", "LOO_predicttrain_GN22873",
                        "LOO_GN22874", "LOO_test_GN22874", "LOO_encoeff_GN22874", "LOO_model_GN22874", "LOO_modelcv_GN22874", "LOO_predicttest_GN22874", "LOO_predicttrain_GN22874",
                        "LOO_GN22875", "LOO_test_GN22875", "LOO_encoeff_GN22875", "LOO_model_GN22875", "LOO_modelcv_GN22875", "LOO_predicttest_GN22875", "LOO_predicttrain_GN22875",
                        "LOO_GN22877", "LOO_test_GN22877", "LOO_encoeff_GN22877", "LOO_model_GN22877", "LOO_modelcv_GN22877", "LOO_predicttest_GN22877", "LOO_predicttrain_GN22877",
                        "LOO_GN22878", "LOO_test_GN22878", "LOO_encoeff_GN22878", "LOO_model_GN22878", "LOO_modelcv_GN22878", "LOO_predicttest_GN22878", "LOO_predicttrain_GN22878",
                        "LOO_GN22879", "LOO_test_GN22879", "LOO_encoeff_GN22879", "LOO_model_GN22879", "LOO_modelcv_GN22879", "LOO_predicttest_GN22879", "LOO_predicttrain_GN22879",
                        "LOO_GN22880", "LOO_test_GN22880", "LOO_encoeff_GN22880", "LOO_model_GN22880", "LOO_modelcv_GN22880", "LOO_predicttest_GN22880", "LOO_predicttrain_GN22880",
                        "LOO_GN22881", "LOO_test_GN22881", "LOO_encoeff_GN22881", "LOO_model_GN22881", "LOO_modelcv_GN22881", "LOO_predicttest_GN22881", "LOO_predicttrain_GN22881",
                        "LOO_GN22883", "LOO_test_GN22883", "LOO_encoeff_GN22883", "LOO_model_GN22883", "LOO_modelcv_GN22883", "LOO_predicttest_GN22883", "LOO_predicttrain_GN22883",
                        "LOO_GN22884", "LOO_test_GN22884", "LOO_encoeff_GN22884", "LOO_model_GN22884", "LOO_modelcv_GN22884", "LOO_predicttest_GN22884", "LOO_predicttrain_GN22884",
                        "LOO_GN22885", "LOO_test_GN22885", "LOO_encoeff_GN22885", "LOO_model_GN22885", "LOO_modelcv_GN22885", "LOO_predicttest_GN22885", "LOO_predicttrain_GN22885",
                        "LOO_GN22887", "LOO_test_GN22887", "LOO_encoeff_GN22887", "LOO_model_GN22887", "LOO_modelcv_GN22887", "LOO_predicttest_GN22887", "LOO_predicttrain_GN22887",
                        "LOO_GN22888", "LOO_test_GN22888", "LOO_encoeff_GN22888", "LOO_model_GN22888", "LOO_modelcv_GN22888", "LOO_predicttest_GN22888", "LOO_predicttrain_GN22888",
                        "LOO_GN22889", "LOO_test_GN22889", "LOO_encoeff_GN22889", "LOO_model_GN22889", "LOO_modelcv_GN22889", "LOO_predicttest_GN22889", "LOO_predicttrain_GN22889",
                        "LOO_GN22891", "LOO_test_GN22891", "LOO_encoeff_GN22891", "LOO_model_GN22891", "LOO_modelcv_GN22891", "LOO_predicttest_GN22891", "LOO_predicttrain_GN22891",
                        "LOO_GN22892", "LOO_test_GN22892", "LOO_encoeff_GN22892", "LOO_model_GN22892", "LOO_modelcv_GN22892", "LOO_predicttest_GN22892", "LOO_predicttrain_GN22892",
                        "LOO_GN22893", "LOO_test_GN22893", "LOO_encoeff_GN22893", "LOO_model_GN22893", "LOO_modelcv_GN22893", "LOO_predicttest_GN22893", "LOO_predicttrain_GN22893",
                        "LOO_GN22894", "LOO_test_GN22894", "LOO_encoeff_GN22894", "LOO_model_GN22894", "LOO_modelcv_GN22894", "LOO_predicttest_GN22894", "LOO_predicttrain_GN22894")))

GN22647_coeff_df<-data.frame(LOO_encoeff_GN22647[1:nrow(LOO_encoeff_GN22647),])
colnames(GN22647_coeff_df)<-"coeff"
nrow(subset(GN22647_coeff_df, abs(coeff) > 0))

GN22649_coeff_df<-data.frame(LOO_encoeff_GN22649[1:nrow(LOO_encoeff_GN22649),])
colnames(GN22649_coeff_df)<-"coeff"
nrow(subset(GN22649_coeff_df, abs(coeff) > 0))

GN22650_coeff_df<-data.frame(LOO_encoeff_GN22650[1:nrow(LOO_encoeff_GN22650),])
colnames(GN22650_coeff_df)<-"coeff"
nrow(subset(GN22650_coeff_df, abs(coeff) > 0))

GN22651_coeff_df<-data.frame(LOO_encoeff_GN22651[1:nrow(LOO_encoeff_GN22651),])
colnames(GN22651_coeff_df)<-"coeff"
nrow(subset(GN22651_coeff_df, abs(coeff) > 0))

GN22652_coeff_df<-data.frame(LOO_encoeff_GN22652[1:nrow(LOO_encoeff_GN22652),])
colnames(GN22652_coeff_df)<-"coeff"
nrow(subset(GN22652_coeff_df, abs(coeff) > 0))

GN22867_coeff_df<-data.frame(LOO_encoeff_GN22867[1:nrow(LOO_encoeff_GN22867),])
colnames(GN22867_coeff_df)<-"coeff"
nrow(subset(GN22867_coeff_df, abs(coeff) > 0))

GN22868_coeff_df<-data.frame(LOO_encoeff_GN22868[1:nrow(LOO_encoeff_GN22868),])
colnames(GN22868_coeff_df)<-"coeff"
nrow(subset(GN22868_coeff_df, abs(coeff) > 0))

GN22869_coeff_df<-data.frame(LOO_encoeff_GN22869[1:nrow(LOO_encoeff_GN22869),])
colnames(GN22869_coeff_df)<-"coeff"
nrow(subset(GN22869_coeff_df, abs(coeff) > 0))

GN22870_coeff_df<-data.frame(LOO_encoeff_GN22870[1:nrow(LOO_encoeff_GN22870),])
colnames(GN22870_coeff_df)<-"coeff"
nrow(subset(GN22870_coeff_df, abs(coeff) > 0))

GN22871_coeff_df<-data.frame(LOO_encoeff_GN22871[1:nrow(LOO_encoeff_GN22871),])
colnames(GN22871_coeff_df)<-"coeff"
nrow(subset(GN22871_coeff_df, abs(coeff) > 0))

GN22872_coeff_df<-data.frame(LOO_encoeff_GN22872[1:nrow(LOO_encoeff_GN22872),])
colnames(GN22872_coeff_df)<-"coeff"
nrow(subset(GN22872_coeff_df, abs(coeff) > 0))

GN22873_coeff_df<-data.frame(LOO_encoeff_GN22873[1:nrow(LOO_encoeff_GN22873),])
colnames(GN22873_coeff_df)<-"coeff"
nrow(subset(GN22873_coeff_df, abs(coeff) > 0))

GN22874_coeff_df<-data.frame(LOO_encoeff_GN22874[1:nrow(LOO_encoeff_GN22874),])
colnames(GN22874_coeff_df)<-"coeff"
nrow(subset(GN22874_coeff_df, abs(coeff) > 0))

GN22875_coeff_df<-data.frame(LOO_encoeff_GN22875[1:nrow(LOO_encoeff_GN22875),])
colnames(GN22875_coeff_df)<-"coeff"
nrow(subset(GN22875_coeff_df, abs(coeff) > 0))

GN22877_coeff_df<-data.frame(LOO_encoeff_GN22877[1:nrow(LOO_encoeff_GN22877),])
colnames(GN22877_coeff_df)<-"coeff"
nrow(subset(GN22877_coeff_df, abs(coeff) > 0))

GN22878_coeff_df<-data.frame(LOO_encoeff_GN22878[1:nrow(LOO_encoeff_GN22878),])
colnames(GN22878_coeff_df)<-"coeff"
nrow(subset(GN22878_coeff_df, abs(coeff) > 0))

GN22879_coeff_df<-data.frame(LOO_encoeff_GN22879[1:nrow(LOO_encoeff_GN22879),])
colnames(GN22879_coeff_df)<-"coeff"
nrow(subset(GN22879_coeff_df, abs(coeff) > 0))

GN22880_coeff_df<-data.frame(LOO_encoeff_GN22880[1:nrow(LOO_encoeff_GN22880),])
colnames(GN22880_coeff_df)<-"coeff"
nrow(subset(GN22880_coeff_df, abs(coeff) > 0))

GN22881_coeff_df<-data.frame(LOO_encoeff_GN22881[1:nrow(LOO_encoeff_GN22881),])
colnames(GN22881_coeff_df)<-"coeff"
nrow(subset(GN22881_coeff_df, abs(coeff) > 0))

GN22883_coeff_df<-data.frame(LOO_encoeff_GN22883[1:nrow(LOO_encoeff_GN22883),])
colnames(GN22883_coeff_df)<-"coeff"
nrow(subset(GN22883_coeff_df, abs(coeff) > 0))

GN22884_coeff_df<-data.frame(LOO_encoeff_GN22884[1:nrow(LOO_encoeff_GN22884),])
colnames(GN22884_coeff_df)<-"coeff"
nrow(subset(GN22884_coeff_df, abs(coeff) > 0))

GN22885_coeff_df<-data.frame(LOO_encoeff_GN22885[1:nrow(LOO_encoeff_GN22885),])
colnames(GN22885_coeff_df)<-"coeff"
nrow(subset(GN22885_coeff_df, abs(coeff) > 0))

GN22887_coeff_df<-data.frame(LOO_encoeff_GN22887[1:nrow(LOO_encoeff_GN22887),])
colnames(GN22887_coeff_df)<-"coeff"
nrow(subset(GN22887_coeff_df, abs(coeff) > 0))

GN22888_coeff_df<-data.frame(LOO_encoeff_GN22888[1:nrow(LOO_encoeff_GN22888),])
colnames(GN22888_coeff_df)<-"coeff"
nrow(subset(GN22888_coeff_df, abs(coeff) > 0))

GN22889_coeff_df<-data.frame(LOO_encoeff_GN22889[1:nrow(LOO_encoeff_GN22889),])
colnames(GN22889_coeff_df)<-"coeff"
nrow(subset(GN22889_coeff_df, abs(coeff) > 0))

GN22891_coeff_df<-data.frame(LOO_encoeff_GN22891[1:nrow(LOO_encoeff_GN22891),])
colnames(GN22891_coeff_df)<-"coeff"
nrow(subset(GN22891_coeff_df, abs(coeff) > 0))

GN22892_coeff_df<-data.frame(LOO_encoeff_GN22892[1:nrow(LOO_encoeff_GN22892),])
colnames(GN22892_coeff_df)<-"coeff"
nrow(subset(GN22892_coeff_df, abs(coeff) > 0))

GN22893_coeff_df<-data.frame(LOO_encoeff_GN22893[1:nrow(LOO_encoeff_GN22893),])
colnames(GN22893_coeff_df)<-"coeff"
nrow(subset(GN22893_coeff_df, abs(coeff) > 0))

GN22894_coeff_df<-data.frame(LOO_encoeff_GN22894[1:nrow(LOO_encoeff_GN22894),])
colnames(GN22894_coeff_df)<-"coeff"
nrow(subset(GN22894_coeff_df, abs(coeff) > 0))



save(list = ls(all.names = TRUE), file = "/dir/work/ZEBRA_SHARK/R/glmnet_all_pearson_ageassoc.RData")
