##############################################
# Volcano Plot of Differentially Expressed Genes
##############################################

library(EnhancedVolcano)

# Load differential expression results
deg_data <- read.csv("./results/deseq2_results_all.csv", row.names = 1)

# Generate volcano plot
EnhancedVolcano(deg_data,
  lab = rownames(deg_data),                  # Gene labels
  x = 'log2FoldChange',                      # X-axis: log2 fold change
  y = 'padj',                                # Y-axis: adjusted p-value
  title = 'Volcano Plot of Differentially Expressed Genes',
  xlab = 'Log2 Fold Change',
  ylab = '-Log10 Adjusted P-value',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 3,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 0.8,
  legendPosition = 'right',
  legendLabSize = 10,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  gridlines.major = FALSE,
  gridlines.minor = FALSE
)
