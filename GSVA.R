# Load necessary libraries
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)

# Input file paths
expFile <- "merge.txt"
groupFile <- "Cluster.txt"
gmtFile <- "c2.cp.kegg.v7.2.symbols.gmt"

# Read expression data
rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, -1]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

# Load gene sets
geneSets <- getGmt(gmtFile, geneIdType = SymbolIdentifier())

# GSVA analysis
gsvaResult <- gsva(data,
                   geneSets,
                   min.sz = 10,
                   max.sz = 500,
                   verbose = TRUE,
                   parallel.sz = 1)

# Save GSVA results
gsvaOut <- rbind(id = colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file = "gsvaOut.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# Read group information
groupInfo <- read.table(groupFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Match samples
gsvaResult <- t(gsvaResult)
commonSamples <- intersect(rownames(gsvaResult), rownames(groupInfo))
gsvaResult <- gsvaResult[commonSamples, , drop = FALSE]
groupInfo <- groupInfo[commonSamples, , drop = FALSE]
mergedData <- cbind(gsvaResult, groupInfo)
projectID <- gsub("(.*?)\\_.*", "\\1", rownames(mergedData))
mergedData <- cbind(mergedData, Project = projectID)

# Differential analysis
adjPFilter <- 0.05
groupLabels <- as.vector(mergedData$Group)
comparisons <- combn(levels(factor(groupLabels)), 2)

for (i in 1:ncol(comparisons)) {
  treat <- mergedData[mergedData$Group == comparisons[2, i], ]
  control <- mergedData[mergedData$Group == comparisons[1, i], ]
  data <- rbind(control, treat)

  type <- as.vector(data$Group)
  ann <- data[, c(ncol(data), ncol(data) - 1)]
  dataMatrix <- t(data[, -c(ncol(data) - 1, ncol(data))])
  design <- model.matrix(~0 + factor(type))
  colnames(design) <- levels(factor(type))

  fit <- lmFit(dataMatrix, design)
  contrast <- paste0(comparisons[2, i], "-", comparisons[1, i])
  cont.matrix <- makeContrasts(contrast, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)

  allDiff <- topTable(fit2, adjust = "fdr", number = 200000)
  allDiffOut <- rbind(id = colnames(allDiff), allDiff)
  write.table(allDiffOut, file = paste0(contrast, ".all.txt"), sep = "\t", quote = FALSE, col.names = FALSE)

  diffSig <- allDiff[with(allDiff, abs(logFC) > 0.1 & adj.P.Val < adjPFilter), ]
  diffSigOut <- rbind(id = colnames(diffSig), diffSig)
  write.table(diffSigOut, file = paste0(contrast, ".diff.txt"), sep = "\t", quote = FALSE, col.names = FALSE)

  # Heatmap visualization
  colors <- c("#0066FF", "#FF9900", "#FF0000", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")
  ann_colors <- list()
  groupColor <- colors[1:length(levels(factor(groupLabels)))]
  names(groupColor) <- levels(factor(groupLabels))
  ann_colors[["Group"]] <- groupColor[c(comparisons[1, i], comparisons[2, i])]

  termNum <- 20
  diffTerms <- rownames(diffSig)
  if (length(diffTerms) < termNum) termNum <- length(diffTerms)
  selectedTerms <- diffTerms[1:termNum]
  heatmapData <- dataMatrix[selectedTerms, ]

  pdf(file = paste0(contrast, ".heatmap.pdf"), height = 6, width = 10)
  pheatmap(heatmapData,
           annotation = ann,
           annotation_colors = ann_colors,
           color = colorRampPalette(c(rep("blue", 2), "white", rep("red", 2)))(50),
           cluster_cols = FALSE,
           show_colnames = FALSE,
           gaps_col = as.vector(cumsum(table(type))),
           scale = "row",
           fontsize = 10,
           fontsize_row = 7,
           fontsize_col = 10)
  dev.off()
}
