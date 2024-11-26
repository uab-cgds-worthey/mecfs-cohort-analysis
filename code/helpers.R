library(biomaRt)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)

check_order <- function(sample_metadata, counts) {
  # Check if data is ordered properly
  if (!all(rownames(sample_metadata) %in% colnames(counts)) ||
    !all(rownames(sample_metadata) == colnames(counts))) {
    # Order sample_metadata
    sample_metadata <- arrange(sample_metadata, ID)

    # Check again after ordering
    if (!all(rownames(sample_metadata) %in% colnames(counts)) ||
      !all(rownames(sample_metadata) == colnames(counts))) {
      cat("Error: Data could not be ordered properly.")
      return(FALSE)
    }
  }
  return(sample_metadata)
}


retrieve_gene_info <- function(values, filters) {
  # Create biomart object
  biomart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

  # Define attributes to retrieve
  attributes <- c(
    "ensembl_gene_id_version",
    "entrezgene_id",
    "ensembl_gene_id",
    "description",
    "genecards",
    "hgnc_symbol"
  )

  # Retrieve gene information from biomart
  gene_info <- getBM(
    attributes = attributes,
    filters = filters,
    values = values,
    mart = biomart
  )

  return(gene_info)
}


process_and_save_results <- function(data_df, output_file) {
  sorted_data <- data_df[order(data_df$padj), ]
  write.csv(sorted_data, file = output_file)
  sorted_data_df <- as.data.frame(sorted_data)
  return(sorted_data_df)
}


rename_counts_columns <- function(counts_df, metadata_df, id_column, rna_samples_id_column) {
  # Create a named vector for mapping
  name_mapping <- setNames(metadata_df[[id_column]], metadata_df[[rna_samples_id_column]])

  # Check if all column names in counts are present in the name mapping
  if (!all(colnames(counts_df) %in% names(name_mapping))) {
    stop("Not all columns in counts can be mapped to new names.")
  }

  # Use the mapping to rename columns
  counts_df <- counts_df[, names(name_mapping)] # Subset the counts_df to columns that are in name_mapping
  colnames(counts_df) <- name_mapping[colnames(counts_df)]

  return(counts_df)
}


# Function to filter genes by exclusion or inclusion
filter_genes <- function(filtered_data, genes_to_filter, include = FALSE) {
  if (include) {
    return(filtered_data %>% filter(gene_name %in% genes_to_filter))
  } else {
    return(filtered_data %>% filter(!gene_name %in% genes_to_filter))
  }
}

# Function to calculate statistics for each gene
calculate_gene_stats <- function(dds, ensembl_id, intgroup) {
  d <- plotCounts(dds, gene = ensembl_id, intgroup = intgroup, returnData = TRUE)
  d_stats <- d %>%
    group_by(Affected) %>%
    summarise(
      mean_count = mean(count),
      sd_count = sd(count)
    )
  d <- merge(d, d_stats, by = "Affected")
  return(d)
}

# Function to determine regulation status
determine_regulation_status <- function(data) {
  sapply(1:nrow(data), function(i) {
    z_score <- (data$count[i] - data$mean_count[i]) / data$sd_count[i]
    if (z_score > 1) {
      return("Upregulated")
    } else if (z_score < -1) {
      return("Downregulated")
    } else {
      return("Neutral")
    }
  })
}

# Function to create and save individual gene plots
create_gene_plot <- function(dp, gene_name, output_path) {
  plot <- ggboxplot(dp,
                    x = "Affected", y = "count", add = "jitter",
                    color = "Affected", palette = c("red", "navy"),
                    title = gene_name
  ) +
    geom_text_repel(aes(label = rownames(dp))) +
    scale_color_manual(
      name = "Affected Status",
      labels = c("Affected", "Unaffected"),
      values = c("red", "navy")
    )
  ggsave(
    filename = paste0(gene_name, "-subsetted-plot-counts.png"),
    path = output_path,
    plot = plot, dpi = 450
  )
  return(plot)
}

# Function to create and save a faceted plot
create_faceted_plot <- function(data, output_path) {
  plot <- ggplot(data, aes(x = Affected, y = count, color = Affected)) +
    geom_boxplot(alpha = 0.35) +
    geom_jitter(width = 0.2, size = .45) +
    geom_text_repel(aes(label = SampleID), size = 2,
                    max.overlaps = Inf, force = 10, color = "black",
                    segment.color = "gray", min.segment.length = 0,
                    point.padding = 0.1) +
    facet_wrap(~ Gene, scales = "free_y") +
    theme_bw() +
    labs(title = "Gene Expression Counts by Affected Status", x = "Affected Status", y = "Count") +
    scale_color_manual(values = c("red", "navy")) +
    theme(
      legend.position = "bottom",
      plot.margin = margin(10, 10, 10, 10),
      strip.text = element_text(size = 4),
      plot.title = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 4),
      panel.spacing = unit(.35, "lines")
    ) +
    coord_cartesian(clip = "off")
  ggsave(
    filename = "faceted_genes_of_interest_plot_counts.png",
    path = output_path,
    plot = plot, dpi = 2200
  )
  return(plot)
}

# Function to generate heatmap
create_heatmap <- function(mat, gene_info, vsd_coldata, output_path, ann_colors, clusters = 3) {
  rownames(mat) <- gene_info$gene_name[match(rownames(mat), gene_info$Ensembl_ID)]
  df_sub <- as.data.frame(vsd_coldata[, c("Affected", "OverallCategory")])
  png(file.path(output_path, "genes_of_interest_heatmap.png"), width = 12, height = 10, units = "in", res = 1200)
  draw(ComplexHeatmap::pheatmap(
    mat = mat,
    color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
    annotation_col = df_sub,
    annotation_colors = ann_colors,
    fontsize = 6,
    angle_col = c("0"),
    cellheight = 10,
    cellwidth = 24,
    border_color = "black"
  ))
  dev.off()
}

# Main function to process and plot genes of interest
process_genes_of_interest <- function(filtered_data, dds, vsd, genes_to_filter, output_path) {
  gene_results_list <- list()
  gi_all_counts <- data.frame()

  # Filter genes
  filtered_without_genes <- filter_genes(filtered_data, genes_to_filter, include = FALSE)

  for (ensembl_id in filtered_without_genes$Ensembl_ID) {
    d <- calculate_gene_stats(dds, ensembl_id, "Affected")
    regulation_status <- determine_regulation_status(d)

    gene_name <- filtered_without_genes$gene_name[filtered_without_genes$Ensembl_ID == ensembl_id]
    temp_gene_df <- data.frame(Patient = rownames(d), !!gene_name := regulation_status, stringsAsFactors = FALSE)

    if (length(gene_results_list) == 0) {
      gene_results_list[[1]] <- temp_gene_df
    } else {
      gene_results_list[[1]] <- merge(gene_results_list[[1]], temp_gene_df, by = "Patient", all = TRUE)
    }

    dp <- d
    dp$Gene <- gene_name
    dp$Ensembl_ID <- ensembl_id
    dp$SampleID <- rownames(dp)

    gi_all_counts <- rbind(gi_all_counts, dp)
    create_gene_plot(dp, gene_name, output_path)
  }

  faceted_plot <- create_faceted_plot(gi_all_counts, output_path)
  return(list(faceted_plot = faceted_plot, gene_results_df = gene_results_list[[1]]))
}

# Function to filter, order, and save data
process_and_save_results <- function(df, padj_threshold, lfc_threshold = NULL, outpath, filename) {
  # Filter by padj
  filtered_df <- subset(df, padj < padj_threshold)

  # Filter by log2FoldChange if threshold is specified
  if (!is.null(lfc_threshold)) {
    filtered_df <- subset(filtered_df, log2FoldChange >= lfc_threshold | log2FoldChange <= -lfc_threshold)
  }

  # Order by padj
  ordered_df <- filtered_df[order(filtered_df$padj), ]

  # Save to CSV
  write.csv(ordered_df, file = file.path("output", outpath, filename), row.names = FALSE)
}
