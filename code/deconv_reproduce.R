# ------------------------------------------------------------------------------
# CIBERSORTx LM22 Whole Blood Analysis — Job2 (Flexible Group Variable)
# Author: Shaurita D. Hutchins
# Description: Comprehensive immune deconvolution and visualization workflow
# Input: CIBERSORTx_Job2_Results.csv
# ID map: id_map.csv (PaperID <-> CGDSID)Y
# Metadata: Metadata_2024_11_20.csv (with ID column)
# Run details:
# - LM22 leukocyte signature (547 genes)
# - B-mode batch correction
# - Quantile normalization disabled
# - 1000 permutations
# - R version 4.5.1 (2025-06-13)
# ------------------------------------------------------------------------------

# 1. Load libraries ------------------------------------------------------------
library(tidyverse)
library(rstatix)
library(MASS)
library(FSA)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)

# 2. Parameters ----------------------------------------------------------------
group_var <- "Affected"  # Change to "Category" as needed
meta_treatment_cols <- c("GC", "HCQ", "AZA", "MTX", "MMF", "CYC")  # treatments if present

# 3. Load data -----------------------------------------------------------------
fractions <- read_csv("data/CIBERSORTx_Job2_Results.csv", show_col_types = FALSE)
colnames(fractions)[1] <- "CGDSID"
celltype_cols <- setdiff(colnames(fractions), c("CGDSID", "P-value", "Correlation", "RMSE"))
fractions <- fractions %>% mutate(across(all_of(celltype_cols), as.numeric))

# 4. Map to PaperID ------------------------------------------------------------
id_map <- read_csv("data/id_map.csv", show_col_types = FALSE) %>%
  dplyr::select(PaperID, CGDSID) %>%
  dplyr::filter(!is.na(PaperID) & !is.na(CGDSID)) %>%
  dplyr::rename(SampleID = PaperID)

fractions <- fractions %>%
  left_join(id_map, by = "CGDSID") %>%
  relocate(SampleID, .before = CGDSID)

# 5. Merge metadata ------------------------------------------------------------
meta <- read_csv("data/Metadata_2024_11_20.csv", show_col_types = FALSE)
df <- fractions %>% left_join(meta, by = c("SampleID" = "ID"))

if (!group_var %in% names(df)) stop(paste("Grouping variable", group_var, "not found in metadata."))
if (any(is.na(df[[group_var]]))) warning("Some samples missing group information.")

# ------------------------------------------------------------------------------
# PART I — Group Comparisons
# ------------------------------------------------------------------------------

# 6. Summary statistics --------------------------------------------------------
summary_stats <- df %>%
  pivot_longer(all_of(celltype_cols), names_to = "CellType", values_to = "Proportion") %>%
  group_by(CellType, .data[[group_var]]) %>%
  summarise(
    mean = mean(Proportion, na.rm = TRUE),
    sd = sd(Proportion, na.rm = TRUE),
    n = sum(!is.na(Proportion)),
    .groups = "drop"
  )
write_csv(summary_stats, paste0("Job2_CellTypeSummary_By_", group_var, ".csv"))

# 7. Non-parametric comparisons ------------------------------------------------
message(paste("Running Mann–Whitney tests for", group_var, "…"))

df_test <- df %>% filter(!is.na(.data[[group_var]]))
test_results <- tibble(CellType = character(), p = numeric())

for (ct in celltype_cols) {
  tmp <- df_test %>%
    dplyr::select(all_of(group_var), all_of(ct)) %>%
    dplyr::filter(!is.na(.data[[group_var]]), !is.na(.data[[ct]]))

  # make sure both groups are represented
  if (length(unique(tmp[[group_var]])) == 2 && all(table(tmp[[group_var]]) > 1)) {
    pval <- tryCatch({
      stats::wilcox.test(tmp[[ct]] ~ tmp[[group_var]])$p.value
    }, error = function(e) NA_real_)
    test_results <- dplyr::bind_rows(test_results,
                                     tibble(CellType = ct, p = pval))
  }
}
test_results <- test_results %>%
  dplyr::mutate(p.adj = p.adjust(p, method = "BH"),
         significance = case_when(
           is.na(p.adj) ~ "",
           p.adj < 0.001 ~ "***",
           p.adj < 0.01  ~ "**",
           p.adj < 0.05  ~ "*",
           TRUE ~ "ns"
         ))
write_csv(test_results, paste0("Job2_GroupComparisonResults_By_", group_var, ".csv"))

# 8. Ordered logistic regression (example) -------------------------------------
df_quart <- df_test %>%
  mutate(across(all_of(celltype_cols), ~ ntile(., 4), .names = "{.col}_quart"))

if ("Neutrophils_quart" %in% colnames(df_quart)) {
  # build formula safely
  formula_str <- paste("as.factor(Neutrophils_quart) ~", group_var, "+ Sex")
  model <- polr(as.formula(formula_str), data = df_quart, Hess = TRUE)

  print(summary(model))
  print(exp(coef(model)))  # odds ratios
}


# 9. Boxplots with significance ------------------------------------------------
plot_df <- df_test %>%
  pivot_longer(all_of(celltype_cols), names_to = "CellType", values_to = "Proportion")

sig_df <- test_results %>%
  dplyr::mutate(group1 = levels(as.factor(df_test[[group_var]]))[1],
         group2 = levels(as.factor(df_test[[group_var]]))[2],
         y.position = 1.05 * max(plot_df$Proportion, na.rm = TRUE)) %>%
  dplyr::select(CellType, group1, group2, p.adj, significance, y.position) %>%
  dplyr::rename(p = p.adj)

p <- ggboxplot(plot_df, x = group_var, y = "Proportion",
               fill = group_var, color = group_var,
               add = "jitter", palette = "Dark2", outlier.shape = NA) +
  facet_wrap(~ CellType, scales = "free_y", ncol = 4) +
  stat_pvalue_manual(sig_df, label = "significance", tip.length = 0.01,
                     bracket.size = 0.3, hide.ns = TRUE, position = position_dodge(0.8)) +
  labs(title = paste("CIBERSORTx LM22 Cell Proportions by", group_var),
       y = "Estimated proportion", x = NULL) +
  theme_bw() + theme(strip.text = element_text(size = 9, face = "bold"),
                     axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "bottom")
ggsave(paste0("Job2_Boxplots_By_", group_var, ".png"), p, width = 10, height = 8, dpi = 300)

# ------------------------------------------------------------------------------
# PART II — Composition and Correlation Analyses
# ------------------------------------------------------------------------------

# 10. Stacked bar plot ---------------------------------------------------------
palette22 <- colorRampPalette(brewer.pal(12, "Set3"))(length(celltype_cols))
plot_df2 <- df %>%
  pivot_longer(all_of(celltype_cols), names_to = "CellType", values_to = "Proportion") %>%
  group_by(SampleID, .data[[group_var]], CellType) %>%
  summarise(Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop") %>%
  arrange(.data[[group_var]], SampleID) %>%
  dplyr::mutate(SampleID = factor(SampleID, levels = unique(SampleID)))

g <- ggplot(plot_df2, aes(x = SampleID, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.8, color = "white", size = 0.1) +
  scale_fill_manual(values = palette22) +
  labs(title = paste("CIBERSORTx LM22 Immune Composition by", group_var),
       y = "Relative proportion", x = "Sample ID", fill = "Cell type") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, size = 7, color = "grey20"),
        panel.grid.major.x = element_blank(),
        legend.position = "right")
ggsave(paste0("Job2_StackedBar_By_", group_var, ".png"), g, width = 11, height = 7, dpi = 300)

# 11. PCA of immune composition ------------------------------------------------
# Remove constant (zero-variance) columns
pca_data <- df %>% dplyr::select(all_of(celltype_cols)) %>% drop_na()
constant_cols <- names(which(apply(pca_data, 2, sd, na.rm = TRUE) == 0))

if (length(constant_cols) > 0) {
  message("Removing constant columns from PCA: ", paste(constant_cols, collapse = ", "))
  pca_data <- pca_data %>% dplyr::select(-all_of(constant_cols))
}

# Run PCA only if enough variation remains
if (nrow(pca_data) > 2 && ncol(pca_data) > 1) {
  pca_res <- prcomp(pca_data, scale. = TRUE)

  # Map PCA coordinates back to metadata
  pca_df <- as.data.frame(pca_res$x[, 1:2]) %>%
    mutate(
      SampleID = df$SampleID[complete.cases(pca_data)],
      Group = df[[group_var]][complete.cases(pca_data)]
    )

  # Plot PCA with sample labels
  pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = Group, label = SampleID)) +
    geom_point(size = 3, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20, segment.color = "grey50") +
    theme_bw() +
    labs(
      title = paste("PCA of CIBERSORTx Immune Cell Proportions by", group_var),
      subtitle = "Each point represents one individual sample",
      x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "% variance)"),
      y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "% variance)")
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  ggsave(paste0("Job2_PCA_By_", group_var, "_Labeled.png"),
         pca_plot, width = 8, height = 6, dpi = 300)
} else {
  warning("Not enough variable cell types for PCA.")
}



# 12. Correlation matrix (Spearman) -------------------------------------------
library(pheatmap)

# Select numeric immune cell columns only
df_corr <- df %>%
  dplyr::select(all_of(celltype_cols)) %>%
  dplyr::select(where(is.numeric))

# Drop columns that are constant or mostly missing
df_corr <- df_corr[, apply(df_corr, 2, function(x) var(x, na.rm = TRUE) > 0.0 & sum(!is.na(x)) > 2)]

# Compute Spearman correlation
corr_matrix <- suppressWarnings(cor(df_corr, method = "spearman", use = "pairwise.complete.obs"))

# Replace NA correlations (if any) with 0 to allow clustering
corr_matrix[is.na(corr_matrix)] <- 0

# Generate heatmap
pheatmap(
  corr_matrix,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 8,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = paste("Spearman Correlation Between Immune Cell Types (", group_var, ")", sep = "")
)

ggplot2::ggsave(paste0("Job2_CorrelationHeatmap_By_", group_var, ".png"),
                width = 7, height = 6, dpi = 300)


# 13. Z-score pseudoheatmap by treatment --------------------------------------
# Define which metadata columns to use for annotation (edit as needed)
# Ensure pheatmap is loaded
library(pheatmap)

# Step 1: Reconfirm z_mat is numeric
stopifnot(is.numeric(z_mat))

# Step 2: Rebuild annotation cleanly
annotation_df <- df %>%
  dplyr::select(SampleID, Category) %>%
  dplyr::filter(!is.na(Category)) %>%
  dplyr::distinct(SampleID, .keep_all = TRUE)

# Set rownames properly
rownames(annotation_df) <- annotation_df$SampleID
annotation_df <- annotation_df[, "Category", drop = FALSE]

# Convert to a pure data.frame with factors (not tibbles)
annotation_df <- as.data.frame(annotation_df)
annotation_df$Category <- as.factor(annotation_df$Category)

# Step 3: Align rows between annotation and matrix
common_samples <- intersect(rownames(z_mat), rownames(annotation_df))
z_mat <- z_mat[common_samples, , drop = FALSE]
annotation_df <- annotation_df[common_samples, , drop = FALSE]

# Step 4: Double check types
cat("Annotation structure:\n")
str(annotation_df)

# Step 5: Plot heatmap
pheatmap(
  t(z_mat),
  annotation_col = annotation_df,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50),
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_colnames = TRUE,
  show_rownames = TRUE,
  fontsize_row = 8,
  fontsize_col = 7,
  border_color = NA,
  main = "Z-score of Immune Cell Types by Category"
)

ggplot2::ggsave("Job2_CategoryPseudoheatmap.png",
                width = 8, height = 6, dpi = 300)
