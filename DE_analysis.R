

#########################################################################
#                                                                       #                                 
#        RNA-Seq differential analysis expression with DESeq2           #                
#                                                                       #
#                                                                       #
#########################################################################


# load packges
library(ggplot2)
library(pheatmap)
library(factoextra)


#install & load DESeq2 (cannot be ibstalled with conda on a windows machine)
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
samples_normal <- rownames(annot)[annot$type == 'normal']
samples_lesional <- rownames(annot)[annot$type == 'lesional']

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
# transform data to log2CPM (add 1 to aviod taking log of 0)
log2cpm <- log2(cpm + 1)

# plot log2(cpm) distribution for each sample in normal group
plot(density(log2cpm[,samples_normal[1]]))

for (i in samples_normal[2:length(samples_normal)]){
  lines(density(log2cpm[, i]))
}

# the log2cpm distribution of each sample looks more or less similar in the normal group

# plot log2cpm distribution for each sample in lesional group
plot(density(log2cpm[,samples_lesional[1]]))
for (i in samples_lesional[2:length(samples_lesional)]){
  lines(density(log2cpm[, i]))
}
# from the plot its looks like one sample values distribute different than the rest 
hist((colMeans(log2cpm[,samples_lesional])))

which((colMeans(log2cpm[,samples_lesional])) < 3.5)

#sample SRR1146078 has values that are a bit lower than most of samples in the lesion group

# change cpm matrix to dataframe for convience
log2cpm <- as.data.frame(log2cpm)
# look only at the 
plot(density(log2cpm$SRR1146078))

# it is not very crutial because the shape of the distribution looks the same, but I
# will remove this sample (I will still have 94 samples left)

keep = !"SRR1146078" == names(log2cpm)
log2cpm <- log2cpm[, keep]

# also remove tha sample from the annotation table
# copy row names to new column (otherwise they will be lost)
annot$names <- rownames(annot)
# filter rows 
keep <- !(rownames(annot) == "SRR1146078") 

annot <- annot[keep,]

# remove 'names' column
annot <- data.frame(type = annot[, c("type")], row.names = annot$names)



#### use top variable/expressed genes to measure variantion between samples ################

#Compute the mean, variance and cv for each gene and sort them in decreasing order
means <- rowMeans(log2cpm)
vars <- apply(log2cpm,1,var)

par(mfrow=c(1,2),mar = c(5,3,2,1))

#Sort, select and plot the top 500 highly variable genes from the data
vars <- sort( vars , decreasing=T)
top_var <- names(vars) [1:500]
boxplot(t(log2cpm[ top_var[1:20],]),  las=2 , col="grey" , main="log2CPM (top var 20 genes)" ,cex=.2)

#Sort, select and plot the top 500 highly expressed genes from the data
means <- sort( means , decreasing=T)
top_means <- names(means) [1:500]
boxplot(t(log2cpm[ top_means[1:20],]),  las=2 , col="grey" , main="log2CPM (top express 20 genes)" ,cex=.2)

# how many in both groups
sum(top_means %in% top_var)
sum(top_var %in% top_means)

# PCA analysis of log2CPM data
# PC of top 500 most vairable genes
PC <-  prcomp(t(log2cpm[top_var ,]), center = TRUE, scale. = TRUE)

par(mfrow=c(1,2),mar = c(5,3,2,1))

plot(PC$x[,1] , PC$x[,2], cex=2, col=factor(annot$type), xlab="PC1", ylab="PC2", 
     pch=16, main="PCA of top 500 most varaible genes", las=1)
text(PC$x[,1] , PC$x[,2], cex=.7, labels = paste0(rownames(annot)), pos = 3)

plot(PC$x[,3] , PC$x[,4], cex=2, col=factor(annot$type), xlab="PC3", ylab="PC4",
     pch=16, main="PCA of top 500 most varaible genes", las=1)
text(PC$x[,3] , PC$x[,4], cex=.7, labels = paste0(rownames(annot)), pos = 3)

# from the PCA analysis we can see that PC1 is good at differentiating between the two types of samples. however,
# its looks like there is a sample from one type with expression pattern more similar to samples from the other type 

# chack sample grouping with top expressed genes

PC <-  prcomp(t(log2cpm[top_means ,]), center = TRUE, scale. = TRUE)

par(mfrow=c(1,2),mar = c(5,3,2,1))

plot(PC$x[,1] , PC$x[,2], cex=2, col=factor(annot$type), xlab="PC1", ylab="PC2", 
     pch=16, main="PCA of top 500 most varaible genes", las=1)
text(PC$x[,1] , PC$x[,2], cex=.7, labels = paste0(rownames(annot)), pos = 3)

plot(PC$x[,3] , PC$x[,4], cex=2, col=factor(annot$type), xlab="PC3", ylab="PC4",
     pch=16, main="PCA of top 500 most varaible genes", las=1)
text(PC$x[,3] , PC$x[,4], cex=.7, labels = paste0(rownames(annot)), pos = 3)

# can genearly see the same pattern but here some of the samples from the different groups are closer


# use k-means to cluster the two groups
set.seed(123)

# K-means on samples with top 500 expression. with k=2
top500_exp_k2 <- kmeans(t(log2cpm[top_means,]), 2)

#plot the groups
fviz_cluster(list(data = t(log2cpm[top_means ,]), cluster = top500_exp_k2$cluster),
             ellipse.type = "norm", geom = "point", stand = FALSE,
             palette = "jco", ggtheme = theme_classic())

# we can see nice seperation between the groups

# collect samples that were asigned to cluster 1 and 2
clust1 <- names(which(top500_exp_k2$cluster == 1))

clust2 <- names(which(top500_exp_k2$cluster == 2))


# find which cluster coresponed to 'lesional' and to 'normal'
lesional_samples <- rownames(annot)[annot$type == 'lesional']

normal_samples <- rownames(annot)[annot$type == 'normal']

sum(lesional_samples %in% clust1) #93
sum(lesional_samples %in% clust2) #1


# from the above I can conclude that k-mean cluster 1 is the 'lesional' group and the 'normal'
# group coresponed to cluster 2
sum(normal_samples %in% clust2)

# so the one sample in clust2 that come from 'lesional' should be removed.
lesional_samples[lesional_samples %in% clust2]

# sample "SRR1146216" shold be removed!
keep = !"SRR1146216" == names(log2cpm)
log2cpm <- log2cpm[, keep]

# also remove tha sample from the annotation table
# copy row names to new column (otherwise they will be lost)
annot$names <- rownames(annot)
# filter rows 
keep <- !(rownames(annot) == "SRR1146216") 

annot <- annot[keep,]

# remove 'names' column
annot <- data.frame(type = annot[, c("type")], row.names = annot$names)



# define 'type' column as factor and specify what levels is the 'base' (what we will compare to)
annot$type <- factor(annot$type, levels= c("normal" ,"lesional"))

#check that the order of sample names is the same in the column of the CPM and rows of annotation table
all(colnames(log2cpm) == rownames(annot))

# get the count data based on the filtartion of CPM
counts <- counts[rownames(log2cpm), colnames(log2cpm)]

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

