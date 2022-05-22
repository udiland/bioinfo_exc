

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
library(EnhancedVolcano)
set.seed(123)

#install & load DESeq2 (cannot be ibstalled with conda on a windows machine)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

# set working directory to where the you have the reopsitory
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
annot_samples <- read.table("data/sample-annotation.txt", sep = "\t", header = T, row.names = 1)

# have a look
head(annot_samples)
summary(annot_samples)

#how many samples in each type?
table(annot_samples$type)

'''
lesional   normal 
95       83
'''

#how many samples? (are the same as number of column in counts?)
dim(annot_samples)[1]

#are the names the same in counts and annotation?
sum(colnames(counts) %in% rownames(annot_samples))


#make a vector of sample_id for each condition (normal/)
samples_normal <- rownames(annot_samples)[annot_samples$type == 'normal']
samples_lesional <- rownames(annot_samples)[annot_samples$type == 'lesional']

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
annot_samples$names <- rownames(annot_samples)
# filter rows 
keep <- !(rownames(annot_samples) == "SRR1146078") 
annot_samples <- annot_samples[keep,]

# remove 'names' column
annot_samples <- data.frame(type = annot_samples[, c("type")], row.names = annot_samples$names)



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

plot(PC$x[,1] , PC$x[,2], cex=2, col=factor(annot_samples$type), xlab="PC1", ylab="PC2", 
     pch=16, main="PCA of top 500 most varaible genes", las=1)
# add sample labels (optional)
text(PC$x[,1] , PC$x[,2], cex=.7, labels = paste0(rownames(annot_samples)), pos = 3)

plot(PC$x[,3] , PC$x[,4], cex=2, col=factor(annot_samples$type), xlab="PC3", ylab="PC4",
     pch=16, main="PCA of top 500 most varaible genes", las=1)
text(PC$x[,3] , PC$x[,4], cex=.7, labels = paste0(rownames(annot_samples)), pos = 3)
'''
 from the PCA analysis we can see that PC1 is good at differentiating between the two types of samples while
 PC 2,3,4 do not seperate acording to type. 
 however,its looks like there is a sample from one type with expression pattern more similar to samples from the other type 
'''

# chack sample grouping with top expressed genes

PC <-  prcomp(t(log2cpm[top_means ,]), center = TRUE, scale. = TRUE)

plot(PC$x[,1] , PC$x[,2], cex=2, col=factor(annot_samples$type), xlab="PC1", ylab="PC2", 
     pch=16, main="PCA of top 500 most varaible genes", las=1)
text(PC$x[,1] , PC$x[,2], cex=.7, labels = paste0(rownames(annot_samples)), pos = 3)

plot(PC$x[,3] , PC$x[,4], cex=2, col=factor(annot_samples$type), xlab="PC3", ylab="PC4",
     pch=16, main="PCA of top 500 most varaible genes", las=1)
text(PC$x[,3] , PC$x[,4], cex=.7, labels = paste0(rownames(annot_samples)), pos = 3)

# can genearly see the same pattern but here some of the samples from the different groups are closer in PC1


## use k-means to cluster the two groups ###################################################################
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
lesional_samples <- rownames(annot_samples)[annot_samples$type == 'lesional']

normal_samples <- rownames(annot_samples)[annot_samples$type == 'normal']

sum(lesional_samples %in% clust1) #1  
sum(lesional_samples %in% clust2) #93
# 'lesional' -> clust2, 'normal' -> clust1

# from the above I can conclude that k-mean cluster 2 is the 'lesional' group and the 'normal'
# group coresponed to cluster 1 
#ATTENTION! the cluster annotation is random so the cluster names
#(i.e. cluster1, cluster2) can switch when re-running the code
sum(lesional_samples %in% clust1)

# so the one sample in clust1 that come from 'lesional' should be removed.
lesional_samples[lesional_samples %in% clust1]

# sample "SRR1146216" shold be removed!
keep = !"SRR1146216" == names(log2cpm)
log2cpm <- log2cpm[, keep]

# also remove tha sample from the annotation table
# copy row names to new column (otherwise they will be lost)
annot_samples$names <- rownames(annot_samples)
# filter rows 
keep <- !(rownames(annot_samples) == "SRR1146216") 

annot_samples <- annot_samples[keep,]

# remove 'names' column
annot_samples <- data.frame(type = annot_samples[, c("type")], row.names = annot_samples$names)



# define 'type' column as factor and specify what levels is the 'base' (what we will compare to)
annot_samples$type <- factor(annot_samples$type, levels= c("normal" ,"lesional"))

#check that the order of sample names is the same in the column of the CPM and rows of annotation table
all(colnames(log2cpm) == rownames(annot_samples))

# subset the count data based on the filtartion of CPM
counts <- counts[rownames(log2cpm), colnames(log2cpm)]

# build deseq2 data set object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = annot_samples,
                              design = ~type)

# print summary (and sanity check)
dds

# run DE analysis
dds <- DESeq(dds)

# take the results
res <- results(dds)
res


# read gene annotations
genes_annot <- read.table("data/gene-annotation.txt", sep="\t", header = T, quote = "", fill = F)

# make a df with results and annotations (only for genes that hase annotations)
res_annotated <- merge(as.data.frame(res), genes_annot, by.x=0, by.y=1, all.x=T)
head(res_annotated, 2)

# filter columns and write results
colnames(res_annotated)

res_annotated <- res_annotated[, c(1,3,7,8,9,10)]

colnames(res_annotated)[1] <- "ENSMBL.ID"

write.table(res_annotated, "DESeq_results.tsv", sep = '\t', row.names = F)

# select top 100 most significant and annotaed genes
res_annotated_only <- merge(as.data.frame(res), genes_annot, by.x=0, by.y=1)

# order the table by p-adjusted
res_annotated_only <- res_annotated_only[order(res_annotated_only$padj) ,]

#select top 100
res_annotated_only <- res_annotated_only[1:100 ,]

# merge with log2cpm
top100_sig_cpm <- merge(res_annotated_only, log2cpm, by.x=1, by.y=0, all.x=T)

# remove unnessesary columns
colnames(top100_sig_cpm)

top100_sig_cpm <- top100_sig_cpm[, c(1,11:length(colnames(top100_sig_cpm)))]

rownames(top100_sig_cpm) <- top100_sig_cpm$Row.names

top100_sig_cpm <- top100_sig_cpm[, c(2:length(colnames(top100_sig_cpm)))]

# z-score normalize the log2cpm
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

top100_sig_cpm_norm <- t(apply(top100_sig_cpm, 1, cal_z_score))

# create heatmap
HM <- pheatmap(top100_sig_cpm_norm, annotation_col = annot_samples, show_rownames = F, show_colnames = F)

pdf("Top100_sig_heatmap.pdf")
HM
dev.off()


# plot volcano
#check what is the minimum padj value of top 100 geens in order to color the plot with top 100 significant genes 
res <- as.data.frame(res)
min_padj <- max(sort(res$padj)[1:100])

# volcano plot with top 100 sig geens in red (FC cutoff 1)
vp <- EnhancedVolcano(res, x = "log2FoldChange", y = "padj", lab="", pCutoff = min_padj)

pdf("VolcanoPlot_top_100.pdf")
vp
dev.off()

###### END PART1 ####################################################################################

'''
In this section I will make an interactive graph the show a scatter plot with the mean log2cpm of the two conditions.
significant genes (padj < 0.05) that are up-regulated (FC >1) will be colored blue
and significant genes (padj < 0.05) that are down-regulated (FC < -1) will be colored red
it uses plotly to create html file that will show gene information while hoovering on a point
and by clicking on a point (gene) it will direct the user to the gene page on ENSMBL website.
The plot can be zoom in/out and saved as a picture
'''

# make a table of gene information, differentail expression and log2cpm
expression_cpm <- merge(res_annotated, log2cpm, by.x=1, by.y=0)
#sanity check
dim(expression_cpm)
dim(res_annotated)
dim(log2cpm)
head(expression_cpm)

#calculate log2cpm mean for each type
#add samples names as different column
annot_samples$samples <- rownames(annot_samples)
# create a vector of samples name for normal and lesional
normal_samples <- annot_samples$samples[annot_samples$type == 'normal']
lesional_samples <- annot_samples$samples[annot_samples$type == 'lesional']

# add to main expression_cpm table the mean log2cpm columns for each type
expression_cpm$normal_log2cpm_mean <- rowMeans(expression_cpm[, normal_samples])
expression_cpm$lesional_log2cpm_mean <- rowMeans(expression_cpm[, lesional_samples])

#tag genes as up-regulated, down regulated or not changed by significance
# with padj < 0.05 and |log2FoldChange| > 1
expression_cpm$group <- ifelse(expression_cpm$padj < 0.05 & expression_cpm$log2FoldChange >= 1, 1, 0)
expression_cpm$group <- ifelse(expression_cpm$padj < 0.05 & expression_cpm$log2FoldChange <= -1, -1, expression_cpm$group)
expression_cpm$group <- as.factor(expression_cpm$group)

# how many genes in each group?
table(expression_cpm$group)
'''
   -1     0     1 
 1162 16005  1119 
1162 genes are downregulated in lesional compare noraml
1119 genes are upregulated
16005 do not change significatly
'''

# add text annotation
expression_cpm$text <- paste("ENSMBL.ID: ", expression_cpm$ENSMBL.ID,
                             "\nlog2FoldChange: ", round(expression_cpm$log2FoldChange, 4),
                             "\nSYMBOL: ", expression_cpm$SYMBOL, 
                             "\nGENENAME: ", expression_cpm$GENENAME,
                             "\nPV_adj:", expression_cpm$padj, sep = "")
#add link 
expression_cpm$link <- paste("https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", expression_cpm$ENSMBL.ID, sep="")

# plot log2cpm with significance
plot <- ggplot(data = expression_cpm, aes(x = normal_log2cpm_mean, y = lesional_log2cpm_mean, color=group, text=text)) + 
  geom_point() + theme_classic() + 
  scale_color_manual(values = c("0" = "black", "-1" = "tomato2","1"= "skyblue4")) +
  theme(legend.position = "none", axis.line = element_blank()) + xlab("Normal [log2 CPM]") + 
  ylab("Lesional [log2 CPM]")
        


p2 <- ggplotly(plot, tooltip = c('text'))


# add link to each color group
for(i in 1:length(p2$x$data)){
  if(p2$x$data[[i]]$name %in% levels(expression_cpm$group)){
    p2$x$data[[i]]$customdata = expression_cpm$link
  }
}


js <- "
function(el) {
  el.on('plotly_click', function(d) {
    var url = d.points[0].customdata;
    //url
    window.open(url);
  });
}"


p3 <- ggplotly(onRender(p2, js))


saveWidget(p3, "log2cpm_normal_vs_lesional.html")
