# ER pathway of DEGs

install.packages("BiocManager")
BiocManager::install(c("pathview", "org.At.tair.db", "AnnotationDbi"))


setwd("C:\\Users\\p0095402\\Desktop\\DESeq2\\DEseq2_shoot")

# Load required libraries
library(pathview)
library(org.At.tair.db)
library(AnnotationDbi)
library(dplyr)

# Load DEGs data
deg_data <- read.csv("DESeq2_significant_2_LNP1OEvsRNAi_shoots_results.csv", header = TRUE)  # Adjust sep if needed

# Ensure correct column names
colnames(deg_data) <- c("X", "BaseMean", "LogFC", "lfcSE", "Stat", "PValue", "Padj")

# Convert TAIR Gene IDs (AT1G...) to Entrez IDs
deg_data$Entrez <- mapIds(org.At.tair.db, 
                          keys = deg_data$X, 
                          column = "ENTREZID", 
                          keytype = "TAIR", 
                          multiVals = "first")

# Remove NA values (genes that couldn’t be mapped)
deg_data <- na.omit(deg_data)

# Prepare named vector for Pathview
gene_data <- deg_data$LogFC
names(gene_data) <- deg_data$Entrez  # Now using Entrez IDs

# Run Pathview with the converted Entrez IDs
pathview(gene.data = gene_data,
         pathway.id = "ath04141",
         species = "ath",
         out.suffix = "Arabidopsis_DEGs",
         kegg.native = TRUE,
         low = "green", high = "red", mid = "white")  # Green = Down, Red = Up
