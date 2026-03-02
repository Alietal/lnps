# RNAi vs WT shoots
library(DESeq2)
library(pheatmap)
library(ggplot2)

setwd("C:\\Users\\p0095402\\Desktop\\DESeq2\\DEseq2_shoot")

# Load count data
counts <- read.table("Combined_counts.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


#save file as excel to look individual LNPs
install.packages("writexl")
library(writexl)

# Save the data frame to an Excel file
write_xlsx(counts, "raw_data_shoots.xlsx")

# Define sample metadata for WT and RNAi roots (adjust according to your actual sample names)
counts_data <- counts[, c("WT_1S.bam", "WT_2S.bam", "WT_3S.bam",
  "LNP1_OE_1S.bam", "LNP1_OE_2S.bam", "LNP1_OE_3S.bam")]


# Create a metadata frame to describe your samples and conditions.
coldata <- data.frame(
  row.names = colnames(counts_data),
  condition = rep(c("WT-S", "LNP1_OE-S"), each = 3))
View(coldata)
# Convert the condition column in your coldata data frame to factors before creating the DESeqDataSet object.
coldata$condition <- as.factor(coldata$condition)


# prepare the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
results <- results(dds)
View(dds)

# View top differentially expressed genes (adjust as needed)
top_genes <- head(results, n=10)
top_genes

#MA plot 

plotMA(results, main="WT vs LNP1 OE shoot", ylim=c(-2, 2))

# Volcano Plot
volcanoData <- as.data.frame(results)
volcanoData$padj[is.na(volcanoData$padj)] <- 1  # Replace NA with 1 for plotting
volcanoData

# Generate Volcano Plot
library(EnhancedVolcano)
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "WT vs OE shoot")

ggplot(volcanoData, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha=0.5) +
  scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) +
  labs(title="LNP1 OE vs RNAi shoot", x="Log2 Fold Change", y="-Log10 Adjusted P-value") +
  theme_minimal()

# Heatmap for top differentially expressed genes
top_genes <- rownames(head(results[order(results$padj), ], 50))  # Adjust number as needed
normalized_counts <- counts(dds, normalized=TRUE)
selected_counts <- normalized_counts[top_genes, ]

pheatmap(log2(selected_counts + 1), cluster_rows=TRUE, cluster_cols=TRUE, 
         annotation_col=results, main="Heatmap of Top Differentially Expressed Genes")

# PCA Plot
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")

# Expression Profile for Individual Genes
gene <- "AT1G23810"  # Replace with your gene of interest
plotCounts(dds, gene=gene, intgroup="condition", main=paste("Expression of", gene))

# Boxplot for a specific gene
data <- plotCounts(dds, gene=gene, intgroup="condition", returnData=TRUE)

ggplot(data, aes(x=condition, y=count)) +
  geom_boxplot(aes(fill=condition)) +
  theme_minimal() +
  ggtitle(paste("Expression of", gene))


top_genes <- head(order(results$padj), 50)  # Select top 30 genes
normalized_counts <- counts(dds, normalized = TRUE)
pheatmap(normalized_counts[top_genes, ], cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = TRUE)

# Define new sample names
new_sample_names <- c("LNP1 OE 1", "LNP1 OE 2", "LNP1 OE 3", "RNAi 1", "RNAi 2", "RNAi 3")

# Ensure the length of new_sample_names matches the number of columns in normalized_counts
colnames(normalized_counts) <- new_sample_names

# Generate the heatmap
pheatmap(normalized_counts[top_genes, ], cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = FALSE)





# filter signficant results
significant_results <- results[which(results$padj < 0.05 & abs(results$log2FoldChange) > 1), ]
head(significant_results)

# Remove rows with NA values
significant_results <- na.omit(significant_results)

# Display the first few rows of the cleaned dataframe
head(significant_results)

# Remove rows with NA values
significant_results <- significant_results[complete.cases(significant_results), ]

# Display the first few rows of the cleaned dataframe
head(significant_results)


# Step 2: Filter significant results
significant_results <- results[which(results$padj < 0.05 & abs(results$log2FoldChange) > 1), ]

# Step 3: Remove rows with NA values
significant_results <- na.omit(significant_results)

# Step 4: Save the cleaned dataframe back to a CSV file
write.csv(significant_results, "DESeq2_significant_LNP1OEvsRNAi_shoots_results.csv", row.names = TRUE)

# Display the first few rows of the cleaned dataframe
head(significant_results)


# Significant results
write.csv(as.data.frame(significant_results), file = "DESeq2_significant_LNP1OEvsRNAi_shoots_results.csv")



# normalize value for heatmap

# Perform variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)

# Extract the transformed values for the top genes
vst_counts <- assay(vsd)[top_genes, ]
# Define new sample names
new_sample_names <- c("LNP1 OE 1", "LNP1 OE 2", "LNP1 OE 3", "RNAi 1", "RNAi 2", "RNAi 3")

# Ensure the length of new_sample_names matches the number of columns in normalized_counts
colnames(vst_counts) <- new_sample_names


# Generate the heatmap with VST transformed data
pheatmap(vst_counts, cluster_rows = FALSE, show_rownames = TRUE, 
         cluster_cols = FALSE, main="LNP1 OE vs RNAi shoot")

# Apply log2 transformation to normalized counts
log2_counts <- log2(normalized_counts + 1)  # Adding 1 to avoid log2(0)

# Generate the heatmap with log2 transformed data
pheatmap(log2_counts[top_genes, ], cluster_rows = FALSE, show_rownames = TRUE, 
         cluster_cols = FALSE, main="LNP1 OE vs RNAi shoot")


# Load necessary libraries
library(ggplot2)

library(dplyr)

# Load the data from your CSV file
df <- read.csv("DESeq2_significant_LNP1OEvsRNAi_shoots_results.csv", header = TRUE)

# Ensure that the 'gene' column exists
df <- df %>% rename(gene = "X")

# Add regulation column if not present
df <- df %>%
  mutate(regulation = case_when(
    log2FoldChange > 1 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -1 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Not significant"
  ))

# Filter for only upregulated and downregulated genes
significant_genes <- df %>%
  filter(regulation %in% c("Upregulated", "Downregulated"))

# Select top 10 upregulated and downregulated genes
top_genes <- significant_genes %>%
  group_by(regulation) %>%
  top_n(20, wt = -padj)

# Create bar plot
ggplot(top_genes, aes(x = reorder(gene, log2FoldChange), y = log2FoldChange, fill = regulation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal() +
  labs(title = "Top Differentially Expressed Genes LNP1 OE vs RNAi shoot", x = "Genes", y = "Log2 Fold Change") +
  theme(legend.title = element_blank())



# DESeq2_significant_RNAivsWT_roots_results

# Define significance and log2 fold change thresholds
padj_threshold <- 0.05
log2fc_threshold <- 1

# Filter significant results
significant_results <- results[which(results$padj < padj_threshold & abs(results$log2FoldChange) > log2fc_threshold), ]

# Count upregulated and downregulated genes
upregulated_genes <- sum(significant_results$log2FoldChange > log2fc_threshold, na.rm = TRUE)
downregulated_genes <- sum(significant_results$log2FoldChange < -log2fc_threshold, na.rm = TRUE)

# Display the counts
cat("Number of upregulated genes:", upregulated_genes, "\n")
cat("Number of downregulated genes:", downregulated_genes, "\n")

# Up and down regulate Significant results
write.csv(as.data.frame(significant_genes), file = "Up_down_regulate_LNP1OEvsRNAi_shoots_results.csv")


results$significance <- ifelse(results$padj < padj_threshold & abs(results$log2FoldChange) > log2fc_threshold,
                               ifelse(results$log2FoldChange > log2fc_threshold, "Upregulated", "Downregulated"), "Not Significant")

volcano_plot <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "LNP1 OE vs RNAi shoot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Gene Regulation") +
  theme(plot.title = element_text(hjust = 0.5))

# Print the plot
print(volcano_plot)

# Create a dataframe for plotting
gene_counts <- data.frame(
  Regulation = c("Upregulated", "Downregulated"),
  Count = c(upregulated_genes, downregulated_genes)
)

# Create the bar graph
bar_graph <- ggplot(gene_counts, aes(x = Regulation, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal() +
  labs(title = "LNP1 OE vs RNAi shoot",
       x = "Gene Regulation",
       y = "Number of Genes") +
  theme(plot.title = element_text(hjust = 0.5))
print(bar_graph)

# heatmap

# Get the list of significant genes
significant_genes <- rownames(significant_results)

# Extract the normalized counts for these genes
normalized_counts <- counts(dds, normalized = TRUE)

# Subset the normalized counts to include only significant genes
significant_counts <- normalized_counts[significant_genes, ]

# Log-transform the counts to improve visualization
log_counts <- log2(significant_counts + 1)

# Create the heatmap
pheatmap(log_counts, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = FALSE, 
         show_colnames = TRUE, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Heatmap of Upregulated and Downregulated Genes")

pheatmap(log2_counts[significant_genes, ], cluster_rows = FALSE, show_rownames = TRUE, 
         cluster_cols = FALSE, main="LNP1 OE vs RNAi shoot")



# Load necessary libraries
library(DESeq2)
library(pheatmap)

# Assuming you have already performed the DESeq2 analysis and have the results
# For example:
# dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)
# dds <- DESeq(dds)
# results <- results(dds)

# Define significance and log2 fold change thresholds
padj_threshold <- 0.05
log2fc_threshold <- 1

# Filter significant results
significant_results <- results[which(results$padj < padj_threshold & abs(results$log2FoldChange) > log2fc_threshold), ]

# Sort results by log2FoldChange
sorted_results <- significant_results[order(significant_results$log2FoldChange, decreasing = TRUE), ]

# Select top 20 upregulated and top 20 downregulated genes
top_genes <- c(rownames(sorted_results)[1:100], rownames(sorted_results)[(nrow(sorted_results)-99):nrow(sorted_results)])

# Extract the normalized counts for the selected genes
normalized_counts <- counts(dds, normalized = TRUE)

# Subset the normalized counts to include only selected genes
selected_counts <- normalized_counts[top_genes, ]

# Define new sample names
new_sample_names <- c("LNP1 OE 1", "LNP1 OE 2", "LNP1 OE 3", "RNAi 1", "RNAi 2", "RNAi 3")

# Ensure the length of new_sample_names matches the number of columns in normalized_counts
colnames(selected_counts) <- new_sample_names


# Log-transform the counts to improve visualization
log_counts <- log2(selected_counts + 1)

# Create the heatmap
pheatmap(log_counts, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top 100 Upregulated and Downregulated Genes")

pheatmap(log_counts, cluster_rows = FALSE, show_rownames = TRUE, 
         cluster_cols = FALSE, main="LNP1 OE vs RNAi shoot")

write.csv(as.data.frame(log_counts), file = "log_top100_genes_LNP1OEvsRNAi_shoots_results.csv")


