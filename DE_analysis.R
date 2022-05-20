

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

# read sample annotation (metadata)
annot <- read.table("data/sample-annotation.txt", sep = "\t", header = T)

# have a look
head(annot)
summary(annot)

#how many samples in each type?
table(annot$type)

'''
lesional   normal 
95       83
'''

#how many samples? (are the same as number of column in counts?)
dim(annot)[1]

#are the names the same in counts and annotation?
sum(colnames(counts) %in% annot$sample_id)


#make a vector of sample_id for each condition (normal/)
samples_normal <- annot$sample_id[annot$type == 'normal'] 
samples_lesional <- annot$sample_id[annot$type == 'lesional']

# for each condition mark genes that have CPM > 1 in at least 75% of samples
# names of genes pass the filter in the normal group
genes2keep_normal <- names(which(rowSums(cpm[, samples_normal] > 1) > length(samples_normal)*0.75))

# how many are there?
length(genes2keep_normal) 
# 17880 genes

# names of genes pass the filter in the lesional group
genes2keep_lesional <- names(which(rowSums(cpm[, samples_lesional] > 1) > length(samples_lesional)*0.75))

# how many are there?
length(genes2keep_lesional) 
# 17248 genes

# combine the 2 names and remove duplications
genes2keep <- unique(c(genes2keep_normal, genes2keep_lesional))

# subset the count table according to the filtering condition
cpm <- cpm[genes2keep ,]

# how many genes we have left?
dim(cpm)[1]
#18286
