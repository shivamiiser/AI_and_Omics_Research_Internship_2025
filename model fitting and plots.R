# =====================================================================
#         Final Corrected Script for GSE21359 Analysis
# =====================================================================

# Load Libraries
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)

# --- TROUBLESHOOTING: Run this once to close any stuck devices ---
# Run this a few times until you get an error message.
dev.off()
dev.off()

# Load your pre-processed R objects
load("Results/preprocessed_data.RData")

# --- PROBE-TO-GENE MAPPING ---
# CORRECTED: Use 'filtered_data' which was loaded from the RData file
probe_ids <- rownames(filtered_data)

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"
)

processed_data_df <- filtered_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::mutate(SYMBOL = gene_symbols[PROBEID]) %>%
  dplyr::relocate(SYMBOL, .after = PROBEID)

processed_data_df <- processed_data_df %>%
  dplyr::filter(!is.na(SYMBOL))

expr_only <- processed_data_df %>%
  dplyr::select(-PROBEID, -SYMBOL)

averaged_data <- limma::avereps(expr_only, ID = processed_data_df$SYMBOL)
data <- data.matrix(as.data.frame(averaged_data))

# --- DIFFERENTIAL GENE EXPRESSION ANALYSIS ---
smoking_status <- phenotype_data$`smoking status:ch1`
is_smoker <- grepl("smoker|COPD", smoking_status, ignore.case = TRUE)
groups <- factor(ifelse(is_smoker, "Smoker", "Non-Smoker"))

# Synchronize data and groups
is_na_group <- is.na(groups)
if (sum(is_na_group) > 0) {
  groups <- groups[!is_na_group]
  data <- data[, !is_na_group]
}

# Create design matrix and fit the model
design <- model.matrix(~0 + groups)
colnames(design) <- make.names(levels(groups))
fit_1 <- lmFit(data, design)

# Create contrast matrix
contrast_matrix <- makeContrasts(Smoker_vs_NonSmoker = Smoker - Non.Smoker,
                                 levels = design)
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)
fit_2 <- eBayes(fit_contrast)

# Get DEG results
deg_results <- topTable(fit_2, coef = "Smoker_vs_NonSmoker", number = Inf)

deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated", "No")
))

# --- DATA VISUALIZATION ---

# --- CORRECTED: This block now SAVES the volcano plot to a file ---
# It will NOT appear on your screen. Look for the file in your "Results" folder.
png("Results/volcano_plot_smokers.png", width = 2000, height = 1500, res = 300)

# When saving, you don't need print()
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Smoker vs. Non-Smoker",
       x = "log2 Fold Change", y = "-log10(Adjusted P-value)")

# --- THIS IS THE MISSING LINE ---
dev.off()

print("Volcano plot has been SAVED to your 'Results' folder.")

# To DISPLAY the plot on screen, run this block separately:
print(
  ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No" = "grey")) +
    theme_minimal() +
    labs(title = "Volcano Plot: Smoker vs. Non-Smoker",
         x = "log2 Fold Change", y = "-log10(Adjusted P-value)")
)


# --- Heatmap of Top 25 DEGs ---

# Make sure you have the 'deg_updown' object from the previous step
# Select the top 25 genes with the smallest adjusted p-values
top_25_genes <- head(rownames(deg_updown[order(deg_updown$adj.P.Val), ]), 25)

# Subset the expression data for these top genes
heatmap_data <- data[top_25_genes, ]

# Generate unique column names for the heatmap (e.g., Smoker_1, Smoker_2)
group_char <- as.character(groups)
heatmap_names <- ave(group_char, group_char, FUN = function(x) paste0(x, "_", seq_along(x)))
colnames(heatmap_data) <- heatmap_names

# Save the heatmap to your "Results" folder
png("Results/heatmap_top25_DEGs.png", width = 2000, height = 1500, res = 300)

pheatmap(
  heatmap_data,
  scale = "row", # Use "row" to scale genes and see relative changes (z-score)
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Top 25 Differentially Expressed Genes (Smoker vs. Non-Smoker)"
)

dev.off()

print("Top 25 DEG heatmap has been SAVED to your 'Results' folder.")


# This will give you the exact counts for your summary
print(paste("Number of Upregulated Genes:", nrow(upregulated)))
print(paste("Number of Downregulated Genes:", nrow(downregulated)))