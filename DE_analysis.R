

  #########################################################################
  #                                                                       #                                 
  #        RNA-Seq differential analysis expression with DESeq2           #                
  #                                                                       #
  #                                                                       #
  #########################################################################


#install & load DESeq2

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

# set working directory
setwd("~/<path>/bioinfo_exc")

# read count data 
counts <- read.table("data/counts.txt", sep = "\t", header = TRUE, row.names = 1)

# have a look
head(counts)

# are there missing values?
sum(is.na(counts))

# how many samples/genes?
dim(counts)
# we have 57992 genes and 178 samples

# calculate count per million (CPM) in order to normlize counts by library size 
cpm <- apply(counts,2, function(x) (x/sum(x))*1e6)

