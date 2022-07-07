# Differential gene expression analysis with edgeR
# Mina Peyton
# 16July2020

# Load edgeR package
library(edgeR)

## Load data and annotations 
Raw <- read.table(file = "RSEM_read_counts.txt", 
                   sep = "\t", 
                   stringsAsFactors = FALSE,
                   header = TRUE, 
                   row.names = 1)

AnnotInfo <- read.csv("AnnotTable.csv", stringsAsFactors = FALSE)

## Create edgeR DGE object using count data and annotation info
Groups <- factor(rep(c(1, 2), 6), labels = c("ctrl", "PA"))
DGE <- DGEList(counts = Raw, group = Groups, genes = AnnotInfo)
DGE$samples

## Assign cell donor and treatment groups
colnames(DGE$counts)
Donor <- factor(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6))
Treatment <- factor(rep(c("ctrl", "PA"), 6))
data.frame(Sample = colnames(DGE), Donor, Treatment)

## Create design matrix
design <- model.matrix(~ Donor + Treatment)
rownames(design) <- colnames(DGE)
design

## Remove low expression genes & normalize 
# identify well-expressed genes:
keep <- filterByExpr(DGE, design, min.count = 10, min.total.count = 15)
DGE <- DGE[keep, ] # keep only well-expressed genes

DGE <- calcNormFactors(DGE) 
DGE$samples 

## Estimate dispersion
# calculate overall dispersion
DGE <- estimateGLMCommonDisp(DGE, design, verbose = TRUE) 

# calculate dispersion trend based on gene abundance
DGE <- estimateGLMTrendedDisp(DGE, design) 

# calculate separate dispersion for each gene
DGE <- estimateGLMTagwiseDisp(DGE, design) 

## Visualizing dispersion
plotBCV(DGE)

# Generate normalized log2 transformed CPMs:
CPM <- cpm(DGE, normalized.lib.sizes = TRUE, log = TRUE)
CPM <- as.data.frame(CPM)

# Visualize normalized data
boxplot(CPM, las = 2, ylab = "log2 CPM", main = "Normalized Data")

## GLM likelihood ratio test for identifying DE genes
## Fit a gene-wise negative binomial generalized linear model
fit <- glmFit(DGE, design)
colnames(fit$coefficients)

## Conduct likelihood ratio test to identify DE genes of PA treatment
lrt <- glmLRT(fit, coef = "TreatmentPA") #specify the comparison of interest
topTags(lrt)

## Summary of genes that are up- or down-regulated
de <- decideTestsDGE(lrt, adjust.method = "fdr")
summary(de)

## Export results
Results <- as.data.frame(topTags(lrt, n = dim(DGE)[1]))
write.csv(Results, file = "edgeR_Results.csv", row.names = FALSE)

## MA Plot with all DE genes in red
detags <- rownames(DGE)[as.logical(de)]
plotSmear(lrt, de.tags = detags)
abline(h = c(-1, 1), col = "blue")

# Number of genes with p < 0.05 (uncorrected) and absolute log2FC > 1 (2-fold change):
DE <- Results[Results$PValue < 0.05 & abs(Results$logFC) > 1, ]
dim(DE)[1]

## MA Plot with 738 DE genes in red
detags2 <- DE$ENSEMBL
plotSmear(lrt, de.tags = detags2, cex = 0.5, main = "MA Plot edgeR")
abline(h = c(-1, 1), col = "blue")
legend("bottom", legend = "DE genes (FC > 2, p < 0.05)", pch = 16, 
       col = "red", bty = "n")

## Volcano plot with 738 DE genes in red  
plot(-log10(Results$PValue) ~ Results$logFC, 
     xlab = "log2 fold change", ylab = "-log10 p-value", 
     main = "Effect of Pseudomonas exposure", xlim = c(-10, 10))

points(-log10(Results$PValue[Results$ENSEMBL %in% detags2]) ~ 
        Results$logFC[Results$ENSEMBL %in% detags2], 
        col = "red", pch = 16)

## Heatmap of 738 DE genes (FC > 2, p < 0.05)
library(gplots)
heatmap.2(as.matrix(CPM[rownames(CPM) %in% detags2, ]), 
          scale = "row", trace = "none")


 