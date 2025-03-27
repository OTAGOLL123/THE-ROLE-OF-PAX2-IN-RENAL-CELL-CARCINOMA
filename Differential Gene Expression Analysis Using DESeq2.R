# Load required package
library(DESeq2)

# Define input file paths
count_matrix_file <- "./data/raw_counts.csv"       # Replace with your actual file path
sample_metadata_file <- "./data/sample_metadata.csv"

# Load count matrix and sample metadata
count_data <- read.csv(count_matrix_file, row.names = 1)
sample_data <- read.csv(sample_metadata_file, row.names = 1)

# Check for matching sample names
if (!all(colnames(count_data) == rownames(sample_data))) {
  stop("Sample names in count matrix and metadata do not match. Please check your input files.")
}

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_data,
  design = ~ Condition  # Replace 'Condition' with the appropriate column name in your metadata
)

# Optional: Filter out lowly expressed genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2 pipeline
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Filter significantly differentially expressed genes
sig_res <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]

# Write output
write.csv(as.data.frame(res), "./results/deseq2_results_all.csv")
write.csv(as.data.frame(sig_res), "./results/deseq2_significant_genes.csv")

# Print summary
summary(res)
