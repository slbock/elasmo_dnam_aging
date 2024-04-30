library(dplyr)
library(caret)
library(methylKit)
library(GenomicRanges)
library(genomation)
library(WGCNA)
library(ggplot2)


load("/scratch/sb61937/work/ZEBRA_SHARK/R/zebrashark.5x.meth.destr.all.RData")

# remove invariable
zebrashark.5x.all.removinvar<-read.csv(file="/scratch/sb61937/work/ZEBRA_SHARK/R/zebrashark.5x.all.removinvar.csv")

# now that the dataset is filtered the way I want, I now need to generate a methylBase object from this so I can run differential methylation analysis
# I think the easiest way to do this is generate a GRanges object from the dataframe, then use the selectByOverlap on the methylBase object

# Need to create GRanges from dataframe
zebrashark.5x.all.removinvar.gr<- makeGRangesFromDataFrame(zebrashark.5x.all.removinvar, keep.extra.columns=FALSE, ignore.strand=TRUE, seqinfo=NULL,
                                                         seqnames.field="seqnames",
                                                         start.field="starts",
                                                         end.field="ends",
                                                         starts.in.df.are.0based=FALSE)

zebrashark.5x.meth.destr.subset<-selectByOverlap(zebrashark.5x.meth.destr.all, zebrashark.5x.all.removinvar.gr)

# re-name methylbase object (more simple)
zebrashark.subset<-zebrashark.5x.meth.destr.subset
dim(zebrashark.subset)

## Do not need to generate percent methylation matrix - we already have that

## Determine overlap status for all covered/filtered CpGs with genes and CpG islands, shores, shelves, or open sea

# Gene context
# I think I need to create a GRanges object for the genomic coordinates of each gene part separately
# Specifically, have a object for each of the following: promoter, exon, and intron
# then I need to probably use the countOverlaps function for each one of these with the full CpG dataset
# if this works correctly, there should not be any CpGs with a 1 for more than one feature - b/c theoretically, a single CpG can't be simultaneously in a promoter and an intron??
# I can then use these overlap counts to create a single column in which each CpG is assigned a gene context - promotor, intron, exon, intergenic
# I should also probably run the function annotateWithGeneParts to verify the percentages match, and this use of Genomic Ranges is doing what I think it's doing
# input annotation files do not have underscores in the scaffold names - I processed these previously

annot<-readTranscriptFeatures("/scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Genome/sSteTig4_annot_edited.bed")
exons<-annot$exons
introns<-annot$introns
promoters<-annot$promoters

## Remove underscores in seqnames and seqlevels
seqlevels(zebrashark.5x.all.removinvar.gr)<-sub('NW_', 'NW', seqlevels(zebrashark.5x.all.removinvar.gr))
seqlevels(zebrashark.5x.all.removinvar.gr)<-sub('NC_', 'NC', seqlevels(zebrashark.5x.all.removinvar.gr))

### Determine overlap with promoters
zebrashark.subset.promoter.overlap<-countOverlaps(zebrashark.5x.all.removinvar.gr, promoters, ignore.strand=TRUE)
zebrashark.subset.promoter.overlap.df<-data.frame(zebrashark.subset.promoter.overlap)

# create new column in promoter overlap dataframe that is binary
zebrashark.subset.promoter.overlap.df<- zebrashark.subset.promoter.overlap.df %>%
  mutate(bin = case_when(zebrashark.subset.promoter.overlap.df > 0 ~ 1,
                         zebrashark.subset.promoter.overlap.df == 0 ~ 0))

# number of CpGs overlapping a promoter
sum(zebrashark.subset.promoter.overlap.df$bin)

# percentage of CpGs overlapping a promoter
sum(zebrashark.subset.promoter.overlap.df$bin)/nrow(zebrashark.subset.promoter.overlap.df)

write.csv(zebrashark.subset.promoter.overlap.df, file="zebrashark.promoter.overlap.csv", row.names=TRUE)

### Determine overlap with exons
zebrashark.subset.exon.overlap<-countOverlaps(zebrashark.5x.all.removinvar.gr, exons, ignore.strand=TRUE)
zebrashark.subset.exon.overlap.df<-data.frame(zebrashark.subset.exon.overlap)

# create new column in exon overlap dataframe that is binary
zebrashark.subset.exon.overlap.df<- zebrashark.subset.exon.overlap.df %>%
  mutate(bin = case_when(zebrashark.subset.exon.overlap.df > 0 ~ 1,
                         zebrashark.subset.exon.overlap.df == 0 ~ 0))

# number of CpGs overlapping an exon
sum(zebrashark.subset.exon.overlap.df$bin)

# percentage of CpGs overlapping a exon
sum(zebrashark.subset.exon.overlap.df$bin)/nrow(zebrashark.subset.exon.overlap.df)

write.csv(zebrashark.subset.exon.overlap.df, file="zebrashark.exon.overlap.csv", row.names=TRUE)


### Determine overlap with introns
zebrashark.subset.intron.overlap<-countOverlaps(zebrashark.5x.all.removinvar.gr, introns, ignore.strand=TRUE)
zebrashark.subset.intron.overlap.df<-data.frame(zebrashark.subset.intron.overlap)

# create new column in intron overlap dataframe that is binary
zebrashark.subset.intron.overlap.df<- zebrashark.subset.intron.overlap.df %>%
  mutate(bin = case_when(zebrashark.subset.intron.overlap.df > 0 ~ 1,
                         zebrashark.subset.intron.overlap.df == 0 ~ 0))

# number of CpGs overlapping an intron
sum(zebrashark.subset.intron.overlap.df$bin)

# percentage of CpGs overlapping a intron
sum(zebrashark.subset.intron.overlap.df$bin)/nrow(zebrashark.subset.intron.overlap.df)

write.csv(zebrashark.subset.intron.overlap.df, file="zebrashark.intron.overlap.csv", row.names=TRUE)

# validate overlap results with annotateWithGeneParts function in genomation
annotateWithGeneParts(zebrashark.5x.all.removinvar.gr, annot,
                      strand = FALSE, intersect.chr = FALSE)

### Determine overlap with CpG islands

# read in bed files with coordinates for islands
cpg.islands <- read.table("/scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Genome/sSteTig4_cpgisland_edited.bed")
cpg.islands$name <-"island"
# name the columns
names(cpg.islands) <- c("seqnames", "start", "end", "name")
#convert zero based start in bedfile to 1 based start
cpg.islands$start=cpg.islands$start+1
#convert to granges
cpg.islands.gr<- makeGRangesFromDataFrame(cpg.islands, keep.extra.columns = FALSE,
                         ignore.strand = TRUE, seqinfo = NULL,
                         seqnames.field = "seqnames",
                         start.field = "start",
                         end.field ="end",
                         starts.in.df.are.0based = FALSE)
# make sure the underscores in seqnames and seqlevels are removed
head(cpg.islands.gr)

# determine overlap with CpG islands
zebrashark.subset.island.overlap<-countOverlaps(zebrashark.5x.all.removinvar.gr, cpg.islands.gr, ignore.strand=TRUE)
zebrashark.subset.island.overlap.df<-data.frame(zebrashark.subset.island.overlap)

# create new column in island overlap dataframe that is binary
zebrashark.subset.island.overlap.df<- zebrashark.subset.island.overlap.df %>%
  mutate(bin = case_when(zebrashark.subset.island.overlap.df > 0 ~ 1,
                         zebrashark.subset.island.overlap.df == 0 ~ 0))

# number of CpGs overlapping an island
sum(zebrashark.subset.island.overlap.df$bin)

# percentage of CpGs overlapping an island
sum(zebrashark.subset.island.overlap.df$bin)/nrow(zebrashark.subset.island.overlap.df)

write.csv(zebrashark.subset.island.overlap.df, file="zebrashark.island.overlap.csv", row.names=FALSE)

# validate results
annotateWithFeature(zebrashark.5x.all.removinvar.gr, cpg.islands.gr, strand = FALSE, extend = 0, intersect.chr = FALSE)


# read in bed files with coordinates for shores
cpg.shores <- read.table("/scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Genome/edited_sSteTig4_cpg_shores.bed")
cpg.shores$name <-"shores"

# name the columns
names(cpg.shores) <- c("seqnames", "start", "end", "name")
#convert zero based start in bedfile to 1 based start
cpg.shores$start=cpg.shores$start+1
#convert to granges
cpg.shores.gr<- makeGRangesFromDataFrame(cpg.shores, keep.extra.columns = FALSE,
                                          ignore.strand = TRUE, seqinfo = NULL,
                                          seqnames.field = "seqnames",
                                          start.field = "start",
                                          end.field ="end",
                                          starts.in.df.are.0based = FALSE)
# make sure the underscores in seqnames and seqlevels are removed
head(cpg.shores.gr)

# determine overlap with CpG shores
zebrashark.subset.shore.overlap<-countOverlaps(zebrashark.5x.all.removinvar.gr, cpg.shores.gr, ignore.strand=TRUE)
zebrashark.subset.shore.overlap.df<-data.frame(zebrashark.subset.shore.overlap)

# create new column in shore overlap dataframe that is binary
zebrashark.subset.shore.overlap.df<- zebrashark.subset.shore.overlap.df %>%
  mutate(bin = case_when(zebrashark.subset.shore.overlap.df > 0 ~ 1,
                         zebrashark.subset.shore.overlap.df == 0 ~ 0))

# number of CpGs overlapping a shore
sum(zebrashark.subset.shore.overlap.df$bin)

# percentage of CpGs overlapping a shore
sum(zebrashark.subset.shore.overlap.df$bin)/nrow(zebrashark.subset.shore.overlap.df)

write.csv(zebrashark.subset.shore.overlap.df, file="zebrashark.shore.overlap.csv", row.names=FALSE)

# validate results
annotateWithFeature(zebrashark.5x.all.removinvar.gr, cpg.shores.gr, strand = FALSE, extend = 0, intersect.chr = FALSE)


# read in bed files with coordinates for shelves
cpg.shelves <- read.table("/scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Genome/edited_sSteTig4_cpg_shelves.bed")
cpg.shelves$name <-"shelves"

# name the columns
names(cpg.shelves) <- c("seqnames", "start", "end", "name")
#convert zero based start in bedfile to 1 based start
cpg.shelves$start=cpg.shelves$start+1
# convert to granges
cpg.shelves.gr<- makeGRangesFromDataFrame(cpg.shelves, keep.extra.columns = FALSE,
                                         ignore.strand = TRUE, seqinfo = NULL,
                                         seqnames.field = "seqnames",
                                         start.field = "start",
                                         end.field ="end",
                                         starts.in.df.are.0based = FALSE)
# make sure the underscores in seqnames and seqlevels are removed
head(cpg.shelves.gr)


# determine overlap with CpG shelves
zebrashark.subset.shelf.overlap<-countOverlaps(zebrashark.5x.all.removinvar.gr, cpg.shelves.gr, ignore.strand=TRUE)
zebrashark.subset.shelf.overlap.df<-data.frame(zebrashark.subset.shelf.overlap)

# create new column in shelf overlap dataframe that is binary
zebrashark.subset.shelf.overlap.df<- zebrashark.subset.shelf.overlap.df %>%
  mutate(bin = case_when(zebrashark.subset.shelf.overlap.df > 0 ~ 1,
                         zebrashark.subset.shelf.overlap.df == 0 ~ 0))

# number of CpGs overlapping a shelf
sum(zebrashark.subset.shelf.overlap.df$bin)

# percentage of CpGs overlapping a shelf
sum(zebrashark.subset.shelf.overlap.df$bin)/nrow(zebrashark.subset.shelf.overlap.df)

write.csv(zebrashark.subset.shelf.overlap.df, file="zebrashark.shelf.overlap.csv", row.names=FALSE)

# validate results
annotateWithFeature(zebrashark.5x.all.removinvar.gr, cpg.shelves.gr, strand = FALSE, extend = 0, intersect.chr = FALSE)
