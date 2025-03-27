# Load required package
library(ConsensusClusterPlus)

# Set input and working directory
expression_file <- "./data/expression_matrix.txt"
output_dir <- "./results/consensus_clustering"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

# Read expression data
expr_matrix <- read.table(expression_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
expr_matrix <- as.matrix(expr_matrix)

# Run consensus clustering
max_k <- 9
clustering_results <- ConsensusClusterPlus(
  d = expr_matrix,
  maxK = max_k,
  reps = 50,
  pItem = 0.8,
  pFeature = 1,
  title = output_dir,
  clusterAlg = "pam",
  distance = "euclidean",
  seed = 123456,
  plot = "png"
)

# Select desired number of clusters (adjust based on evaluation)
selected_k <- 2
cluster_assignment <- clustering_results[[selected_k]][["consensusClass"]]
cluster_df <- data.frame(Cluster = cluster_assignment)

# Rename clusters with letters for clarity (optional)
label_map <- LETTERS[1:max_k]
unique_clusters <- levels(factor(cluster_df$Cluster))
cluster_df$Cluster <- label_map[match(cluster_df$Cluster, unique_clusters)]

# Save cluster assignments
cluster_output <- rbind(ID = colnames(cluster_df), cluster_df)
write.table(cluster_output, file = "./cluster_assignments.txt", sep = "\t", quote = FALSE, col.names = FALSE)
