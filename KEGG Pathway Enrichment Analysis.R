# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Set enrichment thresholds
pval_cutoff <- 0.05
qval_cutoff <- 0.05

# Choose coloring variable for plots
color_by <- if (qval_cutoff <= 0.05) "qvalue" else "pvalue"

# Input: vector of gene symbols (replace with your input source)
gene_symbols <- unique(as.vector(input_data[, 2]))  # Replace 'input_data' with your actual data frame

# Convert gene symbols to Entrez IDs
entrez_ids <- mget(gene_symbols, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrez_ids <- as.character(entrez_ids)
conversion_table <- data.frame(Symbol = gene_symbols, EntrezID = entrez_ids)
valid_genes <- entrez_ids[entrez_ids != "NA"]

# Perform KEGG enrichment analysis
kegg_result <- enrichKEGG(
  gene = valid_genes,
  organism = "hsa",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

# Format and filter results
kegg_df <- as.data.frame(kegg_result)
kegg_df$geneID <- sapply(kegg_df$geneID, function(gene_str) {
  ids <- strsplit(gene_str, "/")[[1]]
  paste(conversion_table$Symbol[match(ids, conversion_table$EntrezID)], collapse = "/")
})
kegg_df <- kegg_df[kegg_df$pvalue < pval_cutoff & kegg_df$qvalue < qval_cutoff, ]

# Save filtered results
write.table(kegg_df, file = "./results/kegg_enrichment_filtered.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Determine number of pathways to display
top_n <- min(30, nrow(kegg_df))

# Generate bar plot
pdf(file = "./results/kegg_barplot.pdf", width = 8, height = 7)
barplot(kegg_result, drop = TRUE, showCategory = top_n, label_format = 100, color = color_by)
dev.off()

# Generate dot plot
pdf(file = "./results/kegg_dotplot.pdf", width = 8, height = 7)
dotplot(kegg_result, showCategory = top_n, orderBy = "GeneRatio", label_format = 100, color = color_by)
dev.off()

# Generate gene-pathway network plot
pdf(file = "./results/kegg_cnetplot.pdf", width = 8, height = 5.25)
cnet_plot <- cnetplot(kegg_result, circular = TRUE, showCategory = 5, colorEdge = TRUE)
print(cnet_plot)
dev.off()
