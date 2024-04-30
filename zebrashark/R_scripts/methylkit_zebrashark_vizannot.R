library(impute)
library(dplyr)
library(caret)
library(methylKit)
library(GenomicRanges)
library(genomation)
library(WGCNA)
library(ggplot2)
library(MASS)
library(car)

# load filtered subsets

zebrashark.5x.percmeth.53S <- read.csv("/scratch/sb61937/work/ZEBRA_SHARK/R/filteredByCov_ouput/zebrashark.5x.percmeth.53S.csv")
zebrashark_factors <- read.csv("/scratch/sb61937/work/ZEBRA_SHARK/R/zebrashark_factors.csv")

# make genomic locus the rowname
rownames(zebrashark.5x.percmeth.53S) <- zebrashark.5x.percmeth.53S[,1]
# remove genomic locus as variable
percmeth_53S <- zebrashark.5x.percmeth.53S[,-1]
head(percmeth_53S)
str(percmeth_53S)

# get percent methylation of most significant age-correlated sites
age_cp_pearson_sig <- read.csv("/scratch/sb61937/work/ZEBRA_SHARK/R/age.cp.pearson.sig.csv")
colnames(age_cp_pearson_sig) <- c("locus", "pval", "fdr", "cor")

age_cp_pearson_sig_meth <- merge(age_cp_pearson_sig, percmeth_53S, by.x="locus", by.y="row.names")
write.csv(age_cp_pearson_sig_meth, file="/scratch/sb61937/work/ZEBRA_SHARK/R/age.cp.pearson.sig.meth.csv")

age_cp_spearman_sig <- read.csv("/scratch/sb61937/work/ZEBRA_SHARK/R/age.cp.spearman.sig.csv")
colnames(age_cp_spearman_sig) <- c("locus", "pval", "fdr", "cor")

age_cp_spearman_sig_meth <- merge(age_cp_spearman_sig, percmeth_53S, by.x="locus", by.y="row.names")
write.csv(age_cp_spearman_sig_meth, file="/scratch/sb61937/work/ZEBRA_SHARK/R/age.cp.spearman.sig.meth.csv")

#percmeth_53S_t <-t(percmeth_53S)
#str(percmeth_53S_t)

#_______________________________________________________________________________
#which(apply(percmeth_53S_t, 2, var)==0)

#percmeth_53S_t_sc <-percmeth_53S_t[ , which(apply(percmeth_53S_t, 2, var) != 0)]


# create PCA based on methylation of all filtered sites
#pca_all <- prcomp(t(na.omit(percmeth_53S)))
#summary(pca_all)
#est_age <- zebrashark_factors$est.age
#type_c <- zebrashark_factors$type.c

#pc1 <- data.frame(pca_all$x[, 1])
#pc2 <- data.frame(pca_all$x[, 2])

#pc1.2.df <- cbind(pc1, pc2, est_age, type_c)
#colnames(pc1.2.df) <- c("pc1", "pc2", "est_age", "type_c")

#save(pc1.2.df, file="/scratch/sb61937/work/ZEBRA_SHARK/R/pc1_2_df.Rdata")
#save(pca_all, file="/scratch/sb61937/work/ZEBRA_SHARK/R/pca_all.Rdata")


# create PCA based on methylation of all filtered sites but with zero variance sites removed
# PCA scaled
#pca_all_sc <- prcomp(na.omit(percmeth_53S_t_sc), scale.=TRUE)
#summary(pca_all_sc)
#est_age <- zebrashark_factors$est.age
#type_c <- zebrashark_factors$type.c

#pc1_sc <- data.frame(pca_all_sc$x[, 1])
#pc2_sc <- data.frame(pca_all_sc$x[, 2])

#pc1.2.sc.df <- cbind(pc1_sc, pc2_sc, est_age, type_c)
#colnames(pc1.2.sc.df) <- c("pc1", "pc2", "est_age", "type_c")

#save(pc1.2.sc.df, file="/scratch/sb61937/work/ZEBRA_SHARK/R/pc1_2_sc_df.Rdata")
#save(pca_all_sc, file="/scratch/sb61937/work/ZEBRA_SHARK/R/pca_all_sc.Rdata")
