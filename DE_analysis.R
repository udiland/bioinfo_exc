

#########################################################################
#                                                                       #                                 
#        RNA-Seq differential analysis expression with DESeq2           #                
#                                                                       #
#                                                                       #
#########################################################################


# load packges
library(ggplot2)
library(pheatmap)

#install & load DESeq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

# set working directory
setwd("~/<path>/bioinfo_exc")

# read count data (place column number 1 as rownames)
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

# save count as RDS file
saveRDS(cpm, file = "noraml_lesion_CPM.rds")


# read sample annotation (metadata)
annot <- read.table("data/sample-annotation.txt", sep = "\t", header = T, row.names = 1)

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
sum(colnames(counts) %in% rownames(annot))



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

# save filtered count as RDS file
saveRDS(cpm, file = "noraml_lesion_CPM_filtered.rds")

# check the samples integrety and variation

# plot log10(cpm) distribution for each sample in normal group
plot(density(log10(cpm[,samples_normal[1]] + 0.001)))

for (i in samples_normal[2:length(samples_normal)]){
  lines(density(log10(cpm[, i] + 0.001)))
}
# the log10(cpm) distribution looks more or less similar in the normal group

# plot log10(cpm) distribution for each sample in lesional group
plot(density(log10(cpm[,samples_lesional[1]] + 0.001)))
for (i in samples_lesional[2:length(samples_lesional)]){
  lines(density(log10(cpm[, i] + 0.001)))
}
# from the plot its looks like one sample values distribute different than the rest 
which(colMeans(log10(cpm[,samples_lesional] + 0.001)) < 1)
#sample SRR1146078 has values that are a bit lower than most of samples in the lesion group

# change cpm matrix to dataframe for convience
cpm <- as.data.frame(cpm)
# look only at the 
plot(density(log10(cpm$SRR1146078)))

# it is not very crutial because the shape of the distribution looks normal, but I
# will remove this sample (I will still have 94 samples left)

keep = !"SRR1146078" == names(cpm)
cpm <- cpm[, keep]

# also remove tha sample from the annotation table
# copy row names to new column (otherwise they will be lost)
annot$names <- rownames(annot)
# filter rows 
keep <- !(rownames(annot) == "SRR1146078") 

annot <- annot[keep,]

# remove 'names' column
annot <- data.frame(type = annot[, c("type")], row.names = annot$names)

# define 'type' column as factor and specify what levels is the 'base' (what we will compare to)
annot$type <- factor(annot$type, levels= c("normal" ,"lesional"))

#check that the order of sample names is the same in the column of the CPM and rows of annotation table
all(colnames(cpm) == rownames(annot))

# get the count data based on the filtartion of CPM
counts <- counts[rownames(cpm), colnames(cpm)]

# build deseq2 data set object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = annot,
                              design = ~type)

# print summary (and sanity check)
dds

# run DE analysis
dds <- DESeq(dds)

# take the results
res <- results(dds)
res






sampleDists <- dist(t(cpm[, samples_lesional]))
sampleDistMatrix <- as.matrix( sampleDists )
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists)

