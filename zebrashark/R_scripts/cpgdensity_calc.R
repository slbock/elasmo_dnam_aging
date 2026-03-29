library(impute)
library(dplyr)
library(caret)
library(methylKit)
library(GenomicRanges)
library(genomation)
library(WGCNA)
library(car)
library(ggplot2)


windows_200bp<-read.table("/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/windows.bed")
CpG_sSteTig4<-read.table("/dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Genome/CpG.sSteFas4_shorter.bed")

windows_200bp_gr<- makeGRangesFromDataFrame(windows_200bp, keep.extra.columns=FALSE, ignore.strand=TRUE, seqinfo=NULL,
                                                           seqnames.field="V1",
                                                           start.field="V2",
                                                           end.field="V3",
                                                           starts.in.df.are.0based=TRUE)

CpGs_sSteTig4_gr<- makeGRangesFromDataFrame(CpG_sSteTig4, keep.extra.columns=FALSE, ignore.strand=TRUE, seqinfo=NULL,
                                           seqnames.field="V1",
                                           start.field="V2",
                                           end.field="V3",
                                           starts.in.df.are.0based=TRUE)

# line below counts how many CpGs are in each window
windows_200bp_gr_overlap<-data.frame(countOverlaps(windows_200bp_gr, CpGs_sSteTig4_gr))
# name count column
colnames(windows_200bp_gr_overlap)<-"number.CpGs"
# add count column to windows_200bp dataframe
windows_200bp_count<-cbind(windows_200bp, windows_200bp_gr_overlap)
max(windows_200bp_count$number.CpGs)
mean(windows_200bp_count$number.CpGs)
median(windows_200bp_count$number.CpGs)
sum(windows_200bp_count$number.CpGs)

# assign a bin variable that categorizes windows by CpGs/200bp
# define bins encompassing values of 0 to 100, broken up by 5 (note this is low because we're working with subset data)
windows_200bp_count <- windows_200bp_count %>%
                      mutate(CpG_bin = cut(number.CpGs, breaks=c(seq(-5, 100, 5))))

write.csv(windows_200bp_count, file="/dir/work/ZEBRA_SHARK/R/windows_200bp_count.csv")

cpg_density_window<-ggplot(data=windows_200bp_count, aes(x=CpG_bin)) + geom_bar()
save(cpg_density_window, file="cpg_density_window.RData")

# now I need to map covered, filtered CpGs to windows and
## characterize these covered CpGs by their local CpG density of the window they overlap

# test file of covered CpGs
cov_filter_cpgs<-read.csv("/dir/work/ZEBRA_SHARK/R/zebrashark.5x.all.removinvar.csv")

cov_filter_cpgs_gr<- makeGRangesFromDataFrame(cov_filter_cpgs, keep.extra.columns=FALSE, ignore.strand=TRUE, seqinfo=NULL,
                                        seqnames.field="seqnames",
                                        start.field="starts",
                                        end.field="ends",
                                        starts.in.df.are.0based=FALSE)

cov_filter_cpgs_overlap <- findOverlaps(query=windows_200bp_gr, subject=cov_filter_cpgs_gr)
cov_filter_cpgs_overlap
cov_filter_cpgs_overlap_df<-data.frame(cov_filter_cpgs_overlap)

# create a row number variable for the windows dataframe
windows_200bp_count$row.n<-seq.int(nrow(windows_200bp_count))

# now add a variable to the subjects dataframe with the covered CpGs that has the queryHits number from findOverlaps
cov_filter_cpgs$row.n<-seq.int(nrow(cov_filter_cpgs))
cov_filter_cpgs$queryHits <-cov_filter_cpgs_overlap_df$queryHits[match(cov_filter_cpgs$row.n, cov_filter_cpgs_overlap_df$subjectHits)]

# now that you have the row number of your corresponding windows, you can merge your covered cpgs dataframe with the windows information
cov_filter_cpgs<-merge(cov_filter_cpgs, windows_200bp_count, by.x="queryHits", by.y="row.n", all.x=TRUE, all.y=FALSE)

write.csv(cov_filter_cpgs, file="/dir/work/ZEBRA_SHARK/R/zebrashark.5x.all.removinvar.density")
