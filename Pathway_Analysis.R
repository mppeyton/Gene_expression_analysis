# Gene Expression Data Pathway Analysis in R
# Mina Peyton
# 16Jul2020

library(edgeR)
library(gplots)
library(org.Hs.eg.db)
library(GO.db)

GoResults<-goana(lrt,
                    geneid=lrt$genes$ENTREZID,
                    species="Hs",
                    FDR = .05)

str(GoanaResults)
head(GoanaResults,25)
dim(GoanaResults) 

# adjust the p-values with FDR method

GoResults$FDR.Up<-p.adjust(GoanaResults$P.Up, method = "fdr")
GoResults$FDR.Down<-p.adjust(GoanaResults$P.Down, method = "fdr")

head(GoResults)

## online pathway analysis tools
EdgeRdata <- read.csv("edgeR_Results.csv",1)
head(EdgeRdata)

#ID column of all expressed genes
write.table(EdgeRdata$ENSEMBL,file="Background_GeneIDs.txt",
            sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE) 

#ID column of significant genes
#define DE genes with FDR<0.05, like the cutoff we used for goana()
DE_FDR0.05 <- EdgeRdata[which(EdgeRdata$FDR<0.05),"ENSEMBL"]
write.table(DE_FDR0.05,file="DE_FDR05_GeneIDs.txt",
            sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE) 


# PANTHER 
PANTHER<-EdgeRdata[,c("ENSEMBL","logFC")]
head(PANTHER)
write.table(PANTHER,file="PANTHER_IDandFC.txt", 
            col.names=FALSE,row.names=FALSE,quote=FALSE,sep ="\t")
