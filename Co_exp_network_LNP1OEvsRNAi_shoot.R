# network analysis Gene Co-expression Network Analysis: WGCNA 26.09.24

# Arabidopsis DEGs networking with related to heat stress

install.packages(c("DESeq2", "tidyverse"))

# Load the packages
library(DESeq2)
library(WGCNA)
library(ggplot2)


#Set a working directory:
setwd("C:\\Users\\p0095402\\Desktop\\DESeq2\\Co-Expression_analysis_shoots")

#Load count data
raw_data <- read.table("Combined_counts.txt", header=TRUE, row.names=1, sep="\t")

# Define sample metadata for WT and RNAi shoots (adjust according to your actual sample names)
counts_data <- raw_data[, c("LNP1_OE_1S.bam", "LNP1_OE_2S.bam", "LNP1_OE_3S.bam",
                          "RNAi_1S.bam", "RNAi_2S.bam", "RNAi_3S.bam")]


# Define sample conditions
conditions <- c("LNP1_OE", "LNP1_OE", "LNP1_OE", "RNAi", "RNAi", "RNAi")
col_data <- data.frame(condition = factor(conditions))
rownames(col_data) <- colnames(counts_data)  # Ensure rownames match

# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts_data), 
                              colData = col_data, 
                              design = ~ condition)

# Extract raw counts from DESeqDataSet
raw_counts <- counts(dds, normalized = FALSE)

# Convert raw counts to long format for easier plotting with ggplot

library(tidyr)
library(ggplot2)

raw_counts_df <- as.data.frame(raw_counts)
raw_counts_long <- pivot_longer(raw_counts_df, cols = everything(), names_to = "Sample", values_to = "Expression")

# Density plot
ggplot(raw_counts_long, aes(x = log2(Expression + 1), fill = Sample)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Raw Counts (log2-transformed)", 
       x = "log2(counts + 1)", y = "Density") +
  theme_minimal()



# Run DESeq normalization
dds <- DESeq(dds)

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Convert to data frame for easier handling
normalized_counts_df <- as.data.frame(normalized_counts)

# Define new sample names
new_sample_names <- c("LNP1 OE 1", "LNP1 OE 2", "LNP1 OE 3", "RNAi 1", "RNAi 2", "RNAi 3")

# Ensure the length of new_sample_names matches the number of columns in normalized_counts
colnames(normalized_counts_df) <- new_sample_names



# Save normalized counts from DESeq2
write.csv(normalized_counts_df, "C:\\Users\\p0095402\\Desktop\\DESeq2\\Co-Expression_analysis_shoots\\normalized_WT_LNP1OE_RNAi_shoots.csv", row.names = TRUE)

# visualize the normalized data
# Density Plot for DESeq2 Normalized Counts
normalized_counts_df_long <- normalized_counts_df %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Expression")

ggplot(normalized_counts_df_long, aes(x = Expression, fill = Sample)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of Normalized Counts (DESeq2)") +
  theme_minimal()


# Install the WGCNA package if you haven't yet
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  install.packages("WGCNA")
}
library(WGCNA)


# WGCNA requires clean data 
# Check for missing values
gsg <- goodSamplesGenes(normalized_counts_df, verbose = 3)
gsg$allOK  # Should return TRUE if there are no problems with the data

# Filter out low-expression genes (optional but recommended for WGCNA)
filtered_counts <- normalized_counts_df[rowMeans(normalized_counts_df) > 1, ]
view(normalized_counts_df) 


# Transpose the data because WGCNA expects genes in columns and samples in rows
datExpr <- t(filtered_counts)

# Pick a soft-thresholding power
powers <- c(1:20) # Range of powers to test
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results to visualize the scale-free topology fit
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.9, col = "red")

# Set the soft-thresholding power (based on the previous step, say power = 6)
softPower <- 6

# Construct adjacency matrix
adjacency <- adjacency(datExpr, power = softPower)

# Turn adjacency into a topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)

# Convert TOM into a dissimilarity matrix
dissTOM <- 1 - TOM

# Hierarchical clustering of genes based on the dissimilarity matrix
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Plot the clustering tree (dendrogram)
plot(geneTree, xlab = "", sub = "", main = "Gene Dendrogram", labels = FALSE, hang = 0.01)

# Open a new plotting device with specified width and height
dev.new(width = 10, height = 8)

# Plot the gene dendrogram again
plot(geneTree, xlab = "", sub = "", main = "Gene Dendrogram", labels = FALSE, hang = 0.01)



# Dynamic tree cut to identify modules
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)

# Convert module labels into colors for visualization
dynamicColors <- labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)


# Export network to visualize in Cytoscape
exportNetworkToCytoscape(TOM, edgeFile = "CytoscapeEdges.txt", nodeFile = "CytoscapeNodes.txt",
                         weighted = TRUE, threshold = 0.02, nodeNames = colnames(datExpr),
                         nodeAttr = dynamicColors)


# Sample names corresponding to the replicates in your expression matrix
sampleLabels <- c("LNP1 OE 1", "LNP1 OE 2", "LNP1 OE 3", "RNAi 1", "RNAi 2", "RNAi 3")

# Create a data frame with the experimental conditions for each sample
sampleTraits <- data.frame(
  Condition = c(rep("OE", 3), rep("RNAi", 3))
)

# Make sure the row names of the trait data match the sample labels
rownames(sampleTraits) <- sampleLabels

# Print out the sample trait data frame
print(sampleTraits)

# Check for missing values (NA)
sum(is.na(datExpr))  # This will tell you how many NA values are present

# Check for infinite values
sum(is.infinite(datExpr))  # This will tell you how many infinite values are present

# Remove rows with any NA values
data <- na.omit(datExpr)

# Or remove columns with any NA values
data <- data[, colSums(is.na(datExpr)) == 0]


# Check that sample names match between the trait data and expression data
all(rownames(sampleTraits) == colnames(datExpr))  # Should return TRUE



# Remove rows with any NA values
data <- na.omit(MEs)

# Or remove columns with any NA values
data <- data[, colSums(is.na(MEs)) == 0]


# Correlate module eigengenes with traits
moduleTraitCor <- cor(MEs, sampleTraits, use = "p")

# Convert the condition column to a factor (if it's not already)
sampleTraits$Condition <- as.factor(sampleTraits$Condition)

# Use model.matrix to create dummy variables for each condition
numericTraits <- model.matrix(~Condition - 1, data = sampleTraits)

# Inspect the numeric traits matrix
print(numericTraits)

# Perform the correlation between module eigengenes and numeric traits
moduleTraitCor <- cor(MEs, numericTraits, use = "p")

# Inspect the correlation matrix
print(moduleTraitCor)

# Set margins to something smaller (default is c(5, 4, 4, 2) + 0.1)
par(mar = c(3, 3, 2, 1))  # This reduces the margin size

# Now run your plot command
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(numericTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  main = "Module-Trait Relationships"
)

# Save the plot to a PNG file
png("module_trait_relationships.png", width = 1200, height = 900)

# Generate the heatmap
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(numericTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  main = "Module-Trait Relationships"
)

# Close the PNG device
dev.off()

# let identify modules and number of genes

# Assuming you've already created the gene dendrogram (geneTree) and calculated the dissimilarity (TOM)

# Dynamic Tree Cut to define modules
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                             deepSplit = 2, pamRespectsDendro = FALSE, 
                             minClusterSize = 30)

# Convert numeric labels into colors for visualization
moduleColors <- labels2colors(dynamicMods)

# View module colors for the first few genes
table(moduleColors)  # This shows the number of genes in each module

# Count the number of genes in each module
moduleGeneCount <- table(moduleColors)

# Display the result
print(moduleGeneCount)

# Create a data frame mapping genes to their modules
geneInfo <- data.frame(Gene = rownames(datExpr), Module = moduleColors)

# View the first few entries
head(geneInfo)

# Extract genes belonging to a specific module (e.g., "blue")
blueGenes <- geneInfo[geneInfo$Module == "blue", "Gene"]

# View the first few genes in the blue module
head(blueGenes)

# Create a barplot of the number of genes in each module
barplot(moduleGeneCount, main = "Number of Genes in Each Module", 
        xlab = "Module", ylab = "Gene Count", col = names(moduleGeneCount))

dev.off()






































































































































