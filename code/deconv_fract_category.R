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
library(ggpubr)

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

fractions <- fractions %>%
  left_join(id_map, by = "CGDSID") %>%
  relocate(SampleID, .before = CGDSID)

if (any(is.na(fractions$SampleID))) {
  warning("Some CGDSIDs were not found in id_map.csv")
  print(fractions %>% filter(is.na(SampleID)) %>% select(CGDSID))
}

# 4. Merge phenotype metadata --------------------------------------------------
meta <- read_csv("data/Metadata_2024_11_20.csv", show_col_types = FALSE)

df <- fractions %>%
  left_join(meta, by = c("SampleID" = "ID"))

if (any(is.na(df$Category))) warning("Some samples missing Category information.")

# 5. Summary statistics --------------------------------------------------------
summary_stats <- df %>%
  pivot_longer(all_of(celltype_cols), names_to = "CellType", values_to = "Proportion") %>%
  group_by(CellType, Category) %>%
  summarise(
    mean = mean(Proportion, na.rm = TRUE),
    sd   = sd(Proportion, na.rm = TRUE),
    n    = sum(!is.na(Proportion)),
    .groups = "drop"
  )

write_csv(summary_stats, "Job2_CellTypeSummary_ByCategory.csv")

# 6. Non-parametric comparisons ------------------------------------------------
message("Running non-parametric group comparisons (Kruskal–Wallis or Mann–Whitney)…")

df_test <- df %>% filter(!is.na(Category))

test_results <- tibble(CellType = character(), p = numeric())

for (ct in celltype_cols) {
  tmp <- df_test %>%
    dplyr::select(Category, !!sym(ct)) %>%
    dplyr::filter(!is.na(Category), !is.na(!!sym(ct)))

  if (length(unique(tmp$Category)) == 2) {
    pval <- tryCatch({
      stats::wilcox.test(tmp[[ct]] ~ tmp$Category)$p.value
    }, error = function(e) NA_real_)
  } else if (length(unique(tmp$Category)) > 2) {
    pval <- tryCatch({
      stats::kruskal.test(tmp[[ct]] ~ tmp$Category)$p.value
    }, error = function(e) NA_real_)
  } else {
    pval <- NA_real_
  }

  test_results <- bind_rows(test_results, tibble(CellType = ct, p = pval))
}

test_results <- test_results %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    significance = case_when(
      is.na(p.adj) ~ "",
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

write_csv(test_results, "Job2_CategoryComparisonResults.csv")
message("Saved: Job2_CategoryComparisonResults.csv")

# 7. Ordered logistic regression ----------------------------------------------
df_quart <- df_test %>%
  mutate(across(all_of(celltype_cols), ~ ntile(., 4), .names = "{.col}_quart"))

if (length(unique(df_quart$Category)) >= 2 && length(unique(df_quart$Sex)) >= 2) {
  model <- polr(as.factor(Neutrophils_quart) ~ Category + Sex, data = df_quart, Hess = TRUE)
  print(summary(model))
  print(exp(coef(model)))
}

# 8. Visualization with significance ------------------------------------------
plot_df <- df_test %>%
  pivot_longer(all_of(celltype_cols), names_to = "CellType", values_to = "Proportion")

sig_df <- test_results %>%
  dplyr::mutate(
    group1 = levels(factor(df_test$Category))[1],
    group2 = levels(factor(df_test$Category))[2],
    y.position = 1.05 * max(plot_df$Proportion, na.rm = TRUE)
  ) %>%
  dplyr::select(CellType, group1, group2, p.adj, significance, y.position) %>%
  dplyr::rename(p = p.adj)

p <- ggpubr::ggboxplot(
  plot_df,
  x = "Category",
  y = "Proportion",
  fill = "Category",
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
    title = "CIBERSORTx LM22 Immune Cell Proportions by Category",
    subtitle = "Boxplots with Mann–Whitney / Kruskal–Wallis significance annotations",
    y = "Estimated Proportion",
    x = NULL
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave("Job2_CellType_Boxplots_byCategory.png", p, width = 10, height = 8, dpi = 300)

# 9. Stacked bar ---------------------------------------------------------------
palette22 <- colorRampPalette(brewer.pal(12, "Set3"))(22)

plot_df <- df %>%
  pivot_longer(all_of(celltype_cols), names_to = "CellType", values_to = "Proportion") %>%
  group_by(SampleID, Category, CellType) %>%
  summarise(Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop") %>%
  arrange(Category, SampleID) %>%
  mutate(SampleID = factor(SampleID, levels = unique(SampleID)))

ggplot(plot_df, aes(x = SampleID, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.8, color = "white", size = 0.1) +
  scale_fill_manual(values = palette22) +
  labs(
    title = "CIBERSORTx LM22 Immune Composition per Sample",
    subtitle = "Samples grouped by Category",
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
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8)
  ) +
  # dynamic dividers and labels per Category
  {
    cat_counts <- plot_df %>%
      distinct(SampleID, Category) %>%
      count(Category) %>%
      mutate(x_mid = cumsum(n) - n / 2, x_div = cumsum(n) + 0.5)

    list(
      geom_vline(
        data = cat_counts[-nrow(cat_counts), ],
        aes(xintercept = x_div),
        color = "black",
        linewidth = 0.6,
        linetype = "dashed"
      ),
      geom_text(
        data = cat_counts,
        aes(x = x_mid, y = 1.05, label = Category),
        fontface = "bold",
        size = 3.5,
        inherit.aes = FALSE
      )
    )
  }

# 10. Save merged data ---------------------------------------------------------
write_csv(df, "Job2_CellTypeProportions_Merged.csv")

# ------------------------------------------------------------------------------
# END OF SCRIPT
# ------------------------------------------------------------------------------
