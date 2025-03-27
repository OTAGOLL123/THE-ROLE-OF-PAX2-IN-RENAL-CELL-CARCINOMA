# Load required library
library(ggplot2)

# Define file paths
expression_file <- "./data/expression_matrix.csv"
group_file <- "./data/sample_metadata.csv"
output_file <- "./results/pca_plot.pdf"

# Load expression matrix and sample metadata
expr_data <- read.csv(expression_file, row.names = 1, header = TRUE, check.names = FALSE)
sample_info <- read.csv(group_file, row.names = 1, header = TRUE)

# Log2 transform the expression data
expr_data <- log2(expr_data + 1)

# Transpose expression matrix (samples as rows)
expr_data <- t(expr_data)

# Perform PCA
pca <- prcomp(expr_data, scale. = TRUE)

# Extract first two principal components
pca_scores <- as.data.frame(pca$x[, 1:2])
colnames(pca_scores) <- c("PC1", "PC2")

# Add group information
pca_scores$Group <- sample_info$Group  # Replace 'Group' with actual column name if different

# Calculate explained variance
explained_var <- round(100 * summary(pca)$importance[2, 1:2], 2)
pc1_label <- paste0("PC1 (", explained_var[1], "%)")
pc2_label <- paste0("PC2 (", explained_var[2], "%)")

# Create PCA plot
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(title = "PCA of Gene Expression", x = pc1_label, y = pc2_label) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# Save plot
ggsave(output_file, plot = pca_plot, width = 8, height = 6)
