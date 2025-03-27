# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Input: Vector of differentially expressed gene IDs (Entrez)
gene_list <- c("1234", "5678", "91011")  # Replace with actual Entrez IDs

# Perform GO enrichment analysis
go_results <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",               # Options: "BP", "CC", or "MF"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# View top results
head(go_results)

# Visualization: Top enriched GO terms
dotplot(go_results, showCategory = 10) +
  ggtitle("Top 10 Enriched GO Terms")

barplot(go_results, showCategory = 10, title = "Top 10 Enriched GO Terms")

emapplot(go_results)  # Enrichment map
cnetplot(go_results, showCategory = 5, circular = TRUE, colorEdge = TRUE)  # GO term-gene network

# Save results to file
write.csv(as.data.frame(go_results), file = "./results/go_enrichment_results.csv")
