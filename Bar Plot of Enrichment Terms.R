# Load required package
library(ggplot2)

# Define file paths
input_file <- "./data/enrichment_results.txt"
output_file <- "./results/enrichment_barplot.pdf"

# Read enrichment data
enrichment_data <- read.table(input_file, header = TRUE, sep = "\t", check.names = FALSE)

# Order terms by p-value (descending)
ordered_terms <- enrichment_data[order(enrichment_data$pvalue, decreasing = TRUE), "Term"]
enrichment_data$Term <- factor(enrichment_data$Term, levels = ordered_terms)

# Create bar plot
barplot <- ggplot(data = enrichment_data) +
  geom_bar(aes(x = Term, y = Count, fill = pvalue), stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "pink", high = "lightblue") +
  xlab("Term") +
  ylab("Gene Count") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black")
  )

# Save plot
ggsave(output_file, plot = barplot, width = 7, height = 5)
