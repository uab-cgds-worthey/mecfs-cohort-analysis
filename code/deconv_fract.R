# ------------------------------------------------------------------------------
# CIBERSORTx LM22 Analysis — Job2 (with ID mapping)
# Author: Shaurita D. Hutchins
# Description: Analyze immune cell-type proportions from CIBERSORTx LM22 output
# Input: CIBERSORTx_Job2_Results.csv
# ID map: id_map.csv (PaperID <-> CGDSID)
# Metadata: Metadata_2024_11_20.csv (with ID column)
# Run details:
# - LM22 leukocyte signature (547 genes)
# - B-mode batch correction
# - Quantile normalization disabled
# - 1000 permutations
# ------------------------------------------------------------------------------

# 1. Load libraries ------------------------------------------------------------
library(tidyverse)
library(rstatix)
library(MASS)   # for polr (ordered logistic regression)
library(FSA)    # for Dunn’s test

# 2. Load CIBERSORTx results ---------------------------------------------------
fractions <- read_csv("data/CIBERSORTx_Job2_Results.csv", show_col_types = FALSE)

# Ensure first column is CGDSID
colnames(fractions)[1] <- "CGDSID"

# Remove diagnostic columns
celltype_cols <- setdiff(colnames(fractions), c("CGDSID", "P-value", "Correlation", "RMSE"))

fractions <- fractions %>%
  mutate(across(all_of(celltype_cols), as.numeric))

# 3. Map to PaperID ------------------------------------------------------------
id_map <- read_csv("data/id_map.csv", show_col_types = FALSE) %>%
  dplyr::select(PaperID, CGDSID) %>%
  dplyr::filter(!is.na(PaperID) & !is.na(CGDSID)) %>%
  dplyr::rename(SampleID = PaperID)

# Join mapping to replace CGDSID with PaperID
fractions <- fractions %>%
  left_join(id_map, by = "CGDSID") %>%
  relocate(SampleID, .before = CGDSID)

# Check for unmapped IDs
if (any(is.na(fractions$SampleID))) {
  warning("Some CGDSIDs were not found in id_map.csv")
  print(fractions %>% filter(is.na(SampleID)) %>% select(CGDSID))
}

# 4. Merge phenotype metadata --------------------------------------------------
# Metadata file must have an "ID" column (matches PaperID)
meta <- read_csv("data/Metadata_2024_11_20.csv", show_col_types = FALSE)

df <- fractions %>%
  left_join(meta, by = c("SampleID" = "ID"))

if (any(is.na(df$Affected))) warning("Some samples missing group information.")

# 5. Summary statistics --------------------------------------------------------
summary_stats <- df %>%
  pivot_longer(all_of(celltype_cols), names_to = "CellType", values_to = "Proportion") %>%
  group_by(CellType, Affected) %>%
  summarise(
    mean = mean(Proportion, na.rm = TRUE),
    sd   = sd(Proportion, na.rm = TRUE),
    n    = sum(!is.na(Proportion)),
    .groups = "drop"
  )

write_csv(summary_stats, "Job2_CellTypeSummary_ByGroup.csv")

# 6. Non-parametric comparisons ------------------------------------------------
# Filter out groups with fewer than 2 samples or all-zero data
# 6. Non-parametric comparisons ------------------------------------------------
message("Running non-parametric group comparisons (Mann–Whitney U)…")

df_test <- df %>% dplyr::filter(!is.na(Affected))

# Build an empty results tibble
test_results <- tibble::tibble(CellType = character(), p = numeric())

for (ct in celltype_cols) {
  tmp <- df_test %>%
    dplyr::select(Affected, !!sym(ct)) %>%
    dplyr::filter(!is.na(Affected), !is.na(!!sym(ct)))

  # make sure both groups are represented
  if (length(unique(tmp$Affected)) == 2 &&
      all(table(tmp$Affected) > 1)) {

    pval <- tryCatch({
      stats::wilcox.test(tmp[[ct]] ~ tmp$Affected)$p.value
    }, error = function(e) NA_real_)

    test_results <- dplyr::bind_rows(
      test_results,
      tibble::tibble(CellType = ct, p = pval)
    )
  }
}

# Multiple-testing correction and significance labels
test_results <- test_results %>%
  dplyr::mutate(
    p.adj = p.adjust(p, method = "BH"),
    significance = dplyr::case_when(
      is.na(p.adj) ~ "",
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

readr::write_csv(test_results, "Job2_GroupComparisonResults.csv")
message("Saved: Job2_GroupComparisonResults.csv")


# 7. Ordered logistic regression ----------------------------------------------
# Convert proportions to quartiles
df_quart <- df_test %>%
  mutate(across(all_of(celltype_cols), ~ ntile(., 4), .names = "{.col}_quart"))

# Example model for one representative cell type
model <- polr(as.factor(Neutrophils_quart) ~ Affected + Sex, data = df_quart, Hess = TRUE)
summary(model)
exp(coef(model))  # odds ratios

# 8. Visualization with significance -------------------------------------------
library(ggpubr)

# Prepare data for plotting
plot_df <- df_test %>%
  tidyr::pivot_longer(all_of(celltype_cols),
                      names_to = "CellType",
                      values_to = "Proportion")

# Prepare significance data for ggpubr
# Convert significance table into the format ggpubr expects
sig_df <- test_results %>%
  dplyr::mutate(
    group1 = "Unaffected",
    group2 = "Affected",
    y.position = 1.05 * max(plot_df$Proportion, na.rm = TRUE)  # position for stars
  ) %>%
  dplyr::select(CellType, group1, group2, p.adj, significance, y.position) %>%
  dplyr::rename(p = p.adj)

# Merge and plot
p <- ggpubr::ggboxplot(
  plot_df,
  x = "Affected",
  y = "Proportion",
  fill = "Affected",
  palette = c("#4575b4", "#d73027"),   # blue/red palette for clarity
  outlier.shape = NA
) +
  facet_wrap(~ CellType, scales = "free_y", ncol = 4) +
  ggpubr::stat_pvalue_manual(
    sig_df,
    label = "significance",
    tip.length = 0.01,
    bracket.size = 0.3,
    hide.ns = TRUE
  ) +
  labs(
    title = "CIBERSORTx LM22 Immune Cell Proportions by Affected Status",
    subtitle = "Boxplots with Mann–Whitney significance annotations",
    y = "Estimated Proportion",
    x = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# Display
print(p)

# Save to file
ggplot2::ggsave("Job2_CellType_Boxplots_with_Significance.png", p, width = 10, height = 8, dpi = 300)

# 9. Save cleaned data ---------------------------------------------------------
write_csv(df, "Job2_CellTypeProportions_Merged.csv")

# 10. Stacked bar

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Create 22-color palette
palette22 <- colorRampPalette(brewer.pal(12, "Set3"))(22)

# Prepare data
plot_df <- df %>%
  pivot_longer(all_of(celltype_cols), names_to = "CellType", values_to = "Proportion") %>%
  group_by(SampleID, Affected, CellType) %>%
  summarise(Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop")

# Order samples: unaffected first, then affected
plot_df <- plot_df %>%
  arrange(Affected, SampleID) %>%
  mutate(SampleID = factor(SampleID, levels = unique(SampleID)))

# ---------------------- Unified grouped bar chart ----------------------
ggplot(plot_df, aes(x = SampleID, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.8, color = "white", size = 0.1) +
  scale_fill_manual(values = palette22) +
  labs(
    title = "CIBERSORTx LM22 Immune Composition per Sample",
    subtitle = "Samples grouped by clinical status (Affected vs Unaffected)",
    y = "Relative proportion",
    x = "Sample ID",
    fill = "Cell type"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7, color = "grey20"),
    axis.title.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8)
  ) +
  # Add divider line between groups
  geom_vline(
    xintercept = length(unique(plot_df$SampleID[plot_df$Affected == "Unaffected"])) + 0.5,
    color = "black",
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = length(unique(plot_df$SampleID[plot_df$Affected == "Unaffected"])) / 2,
    y = 1.05,
    label = "Unaffected",
    fontface = "bold",
    size = 3.5
  ) +
  annotate(
    "text",
    x = length(unique(plot_df$SampleID[plot_df$Affected == "Unaffected"])) +
      length(unique(plot_df$SampleID[plot_df$Affected == "Affected"])) / 2,
    y = 1.05,
    label = "Affected",
    fontface = "bold",
    size = 3.5
  )
