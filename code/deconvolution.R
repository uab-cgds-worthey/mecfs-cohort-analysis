# ------------------------------------------------------------------------------
# CIBERSORTx Group Mode Analysis — Job1
# Author: Shaurita D. Hutchins
# Description: Summarize cell-type proportions and representative gene profiles
# ------------------------------------------------------------------------------

# 1. Load libraries ------------------------------------------------------------
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rstatix)

# 2. Define file paths ---------------------------------------------------------
prefix <- "data/CIBERSORTx_Job1_output/CIBERSORTxGEP_Job1_"

fractions_file <- paste0(prefix, "Fractions.txt")
geps_file <- paste0(prefix, "GEPs.txt")
qvals <- read.table(qvals_file, sep = "\t", header = TRUE, row.names = 1,
                    check.names = FALSE, stringsAsFactors = FALSE)

# 3. Load and clean fractions --------------------------------------------------
fractions <- read.table(fractions_file, sep = "\t", header = TRUE, check.names = FALSE)

# Rename the first column as sample ID
colnames(fractions)[1] <- "SampleID"

# Drop the diagnostic metric columns
fraction_cols <- setdiff(colnames(fractions),
                         c("SampleID", "P-value", "Correlation", "RMSE"))

# Convert to numeric for downstream calculations
fractions <- fractions %>%
  mutate(across(all_of(fraction_cols), as.numeric))

# Check data structure
str(fractions)
dim(fractions)

# 4. Load GEPs and Q-values ----------------------------------------------------
geps <- read.delim(geps_file, row.names = 1, check.names = FALSE)
qvals <- read.delim(qvals_file, row.names = 1, check.names = FALSE)

dim(geps)
dim(qvals)

# 5. Summarize estimated cell-type proportions ---------------------------------
fractions_summary <- fractions %>%
  pivot_longer(all_of(fraction_cols), names_to = "CellType", values_to = "Proportion") %>%
  group_by(CellType) %>%
  summarise(mean_prop = mean(Proportion, na.rm = TRUE),
            sd_prop = sd(Proportion, na.rm = TRUE)) %>%
  arrange(desc(mean_prop))

# Plot mean proportions
ggplot(fractions_summary, aes(x = reorder(CellType, -mean_prop), y = mean_prop, fill = CellType)) +
  geom_col() +
  labs(x = "Cell type", y = "Mean estimated proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 6. Identify significant genes per cell type ----------------------------------
# Move rownames to a column before pivoting
sig_genes <- qvals %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "CellType", values_to = "qval") %>%
  filter(qval < 0.05) %>%
  left_join(
    geps %>%
      rownames_to_column("Gene") %>%
      pivot_longer(-Gene, names_to = "CellType", values_to = "Expr"),
    by = c("Gene", "CellType")
  )

# Top 50 expressed genes per cell type
top_genes <- sig_genes %>%
  group_by(CellType) %>%
  top_n(50, Expr) %>%
  distinct(Gene, CellType, .keep_all = TRUE)

# 7. Heatmap of top expressed genes --------------------------------------------
gene_subset <- unique(top_genes$Gene)

# Ensure we only use genes actually present in geps
gene_subset <- intersect(rownames(geps), gene_subset)

if (length(gene_subset) == 0) {
  warning("No matching genes found between GEPs and top_genes; check gene ID formats.")
} else {
  message(length(gene_subset), " genes found for heatmap plotting.")

  mat <- geps[rownames(geps) %in% gene_subset, ]

  # Replace infinite or NA values with 0 before scaling
  mat[!is.finite(mat)] <- 0

  pheatmap(
    mat,
    scale = "row",
    clustering_distance_cols = "correlation",
    show_rownames = FALSE,
    main = "Top expressed genes across cell types"
  )
}


# 8. Optional: GO enrichment per cell type -------------------------------------
enrichment_results <- list()

for (ct in unique(top_genes$CellType)) {
  genes_ct <- top_genes %>% filter(CellType == ct) %>% pull(Gene)
  ego <- enrichGO(
    gene = genes_ct,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  enrichment_results[[ct]] <- ego
  message("Top GO terms for ", ct)
  print(head(ego@result[, c("Description", "p.adjust")], 5))
}

# Example enrichment visualization
if (length(enrichment_results) > 0) {
  dotplot(enrichment_results[[1]], showCategory = 10)
}

# 9. Save key outputs ----------------------------------------------------------
write_csv(fractions_summary, "Job1_CellTypeProportionSummary.csv")
write_csv(top_genes, "Job1_TopGenes_PerCellType.csv")
