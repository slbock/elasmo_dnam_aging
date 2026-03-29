library(impute)
library(dplyr)
library(caret)
library(methylKit)
library(GenomicRanges)
library(genomation)
library(WGCNA)
library(reshape2)
library(glmnet)

load("/dir/work/ZEBRA_SHARK/R/filteredByCov_ouput/zebrashark.5x.percmeth.53S.impute.RData")

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
zebrashark_factors_cc<-read.csv(file="/dir/work/ZEBRA_SHARK/R/zebrashark_factors_cc.csv")

# read in metadata for wild individuals only
zebrashark_factors_w<-read.csv(file="/dir/work/ZEBRA_SHARK/R/zebrashark_factors_wild.csv")

# create dataframe only including age data and GN identifiers
zebrashark_age_cc<-subset(zebrashark_factors_cc, select=c(gn.number, est.age))

# create dataframe only including age data and GN identifiers
zebrashark_age_w<-subset(zebrashark_factors_w, select=c(gn.number, est.age))

# create transposed matrix of methylation data
meth_captive_loc_t <- as.matrix(t(zebrashark.subset.percmeth.captive.df[,1:48]))
meth_wild_loc_t <- as.matrix(t(zebrashark.subset.percmeth.wild.df[,1:18]))

# test for spearman correlations
spearman_captive<-corAndPvalue(meth_captive_loc_t, zebrashark_factors_cc$est.age, method="spearman")
pearson_captive<-corAndPvalue(meth_captive_loc_t, zebrashark_factors_cc$est.age, method="pearson")

# add spearman correlations to original methylation dataset
zebrashark.subset.percmeth.captive.df<-cbind(zebrashark.subset.percmeth.captive.df, cor=spearman_captive$cor, pcor=pearson_captive$cor)
write.csv(zebrashark.subset.percmeth.captive.df, file="/dir/work/ZEBRA_SHARK/R/zebrashark_subset_percmeth_captive.csv")

zebrashark.subset.percmeth.wild.df<-cbind(zebrashark.subset.percmeth.wild.df, cor=spearman_captive$cor, pcor=pearson_captive$cor)
write.csv(zebrashark.subset.percmeth.wild.df, file="/dir/work/ZEBRA_SHARK/R/zebrashark_subset_percmeth_wild.csv")

# subset methylation data to only include age-associated sites
meth_captive_loc_ageassoc<-subset(zebrashark.subset.percmeth.captive.df, abs(pcor) > 0.5)
meth_wild_loc_ageassoc<-subset(zebrashark.subset.percmeth.wild.df, abs(pcor) > 0.5)

# create transposed methylation matrix with age-associated loci only
meth_captive_loc_ageassoc_t<-as.matrix(t(meth_captive_loc_ageassoc[,1:48]))
meth_wild_loc_ageassoc_t<-as.matrix(t(meth_wild_loc_ageassoc[,1:18]))








predict_age_captive<-predict(model, meth_captive_loc_ageassoc_t, s=best_lambda)
predict_age_wild<-predict(model, meth_wild_loc_ageassoc_t, s=best_lambda)

save(list = ls(all.names = TRUE), file = "/dir/work/ZEBRA_SHARK/R/glmnet_all_pearson_wild.RData")
