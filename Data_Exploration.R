# Exploratory Gene Expression Data Analysis in R
# Mina Peyton
# 16Jul2020

# Load sample information file
Sinfo <- read.csv("SampleInfo.csv", stringsAsFactors = FALSE)
Sinfo

# Load count table of aligned reads
Raw <- read.table(file = "RSEM_read_counts.txt", sep = "\t", 
                   stringsAsFactors = FALSE,
                   header = TRUE, row.names = 1)

names(Raw) <- c("M1_ctrl", "M1_PA", "M2_ctrl",	"M2_PA", 
                "M3_ctrl", "M3_PA", "F1_ctrl", "F1_PA", 
                "F2_ctrl", "F2_PA", "F3_ctrl", "F3_PA")

# Count distribution
hist(log2(rowSums(Raw)))

# Only include genes with at least 12 total counts (1 count/sample on average)
Exp <- Raw[rowSums(Raw) > 11, ]
dim(Exp)

# Boxplot of raw counts
boxplot(Exp, las = 2)

# Quick solution - add 1 count to each cell of the dataframe, then take log2
Exp_log2 <- log2(Exp + 1)
boxplot(Exp_log2, las = 2)

# Calculate median gene expression value for each sample
SampleMedians <- apply(Exp_log2, 2, median)
SampleMedians

# Look for systematic differences in counts by treatment or sex
boxplot(SampleMedians ~ Sinfo$Treatment)
t.test(SampleMedians ~ Sinfo$Treatment)

boxplot(SampleMedians ~ Sinfo$Sex)
t.test(SampleMedians ~ Sinfo$Sex)

# Dist calculates distances between rows, use t() to transpose
RawDist <- dist(t(Exp_log2), method = "euclidean")
plot(hclust(RawDist, method = "average"))

# Simple normalization
GrandMedian <- mean(SampleMedians)
CorrectionFactors <- GrandMedian - SampleMedians 
CorrectionFactors

# Create an empty matrix with the same dimensions and names as Exp_log2
ExpNorm <- matrix(nrow = length(rownames(Exp_log2)),
                  ncol = length(colnames(Exp_log2)),
                  dimnames = list(rownames(Exp_log2), colnames(Exp_log2)))

# Loop through each column (sample) and fill in medians adjusted with sample correction factor
for(col in colnames(ExpNorm)){
  ExpNorm[ , col] <- Exp_log2[ , col] + CorrectionFactors[col]
}

# Boxplot of normalized counts
boxplot(ExpNorm, ylab = "log2 counts", 
        main = "Normalized Counts", las = 2)

# Cluster dendrogram of normalized data
NormDist <- dist(t(ExpNorm), method = "euclidean")
plot(hclust(NormDist, method = "average"))

# PCA plot
PCA <- prcomp(t(NormDist))
plot(PCA$x[ , 1], PCA$x[ , 2], pch = 19)

# Use colors to highlight groups
plot(PCA$x[ , 1], PCA$x[ , 2], pch = 19, col = as.factor(Sinfo$Donor))
plot(PCA$x[ , 1], PCA$x[ , 2], pch = 19, col = as.factor(Sinfo$Sex))
plot(PCA$x[ , 1], PCA$x[ , 2], pch = 19, col = as.factor(Sinfo$Treatment))

# Identify top 50 genes with highest SD across all samples
NormSD <- apply(ExpNorm, 1, sd)
TopSD <- names(NormSD[order(NormSD, decreasing = TRUE)[1:50]])

# Make a heatmap of top 50 SD genes
library("gplots")
heatmap.2(as.matrix(ExpNorm[TopSD, ]), 
scale = "row",
trace = "none",
col = rev(redblue(100)))
