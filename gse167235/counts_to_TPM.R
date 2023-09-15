library(stringr)
library(sva)
library(data.table)
library(readxl)
library(DESeq2)
library(DESeq)
library(pamr)
library(ggpubr)
library(Seurat)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(ggplot2)
library(gplots)
library(pca3d)
library(rgl)
library(scatterplot3d)
library(FactoMineR)
library(ggfortify)
library(useful)
library(tidyverse)
library(kableExtra)
library(xfun)
library(psych)
library(limma)
library(calibrate)
library(pheatmap)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
library(biomaRt)
library(GenomicFeatures)
library(xlsx)
setwd("~/gse167235")
ensembl_list<-read.table("tmp01",header=FALSE)
dim(ensembl_list)
human<-useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id","start_position","end_position"), filters="ensembl_gene_id", values=ensembl_list, mart=human)

gene_coords$size=gene_coords$end_position - gene_coords$start_position
head(gene_coords)

a<-"gse167235_counts"
counts01<-data.frame(fread(a),check.names=FALSE,row.names=1)
counts<-as.matrix(counts01)
head(counts)[,1:9]

featureLength01<-read.table("feature_length01",head=TRUE)
featureLength02<-t(featureLength01)
featureLength<-as.vector(featureLength02)

Counts_to_tpm <- function(counts, featureLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- featureLength
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

tmp_counts<-Counts_to_tpm(counts,featureLength)
write.table(tmp_counts,"tmp_counts")


