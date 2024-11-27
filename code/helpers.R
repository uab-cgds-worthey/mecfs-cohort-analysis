library(tidyverse) # Available via CRAN
library(DESeq2) # Available via Bioconductor
library(RColorBrewer) # Available via CRAN
library(pheatmap) # Available via CRAN
library(genefilter) # Available via Bioconductor
library(limma) # Available via Bioconductor
library(gprofiler2) # Available via CRAN
library(biomaRt) # Available via Bioconductor
library(plotly) # Available via CRAN
library(ggpubr) # Available via CRAN
library(rmarkdown) # Available via CRAN
library(clusterProfiler) # Available via Bioconductor
library(org.Hs.eg.db) # Available via Bioconductor
library(ggrepel) # Available via CRAN
library(ReactomePA) # Available via Bioconductor
library(mygene) # Available via Bioconductor
library(DOSE) # Available via Bioconductor
library(enrichR) # Available via Bioconductor
library(STRINGdb) # Available via Bioconductor
library(EnhancedVolcano)
library(ComplexHeatmap)

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
create_heatmap <- function(mat, gene_info, vsd_coldata, output_path, ann_colors) {
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
process_and_save_filtered_results <- function(df, padj_threshold, lfc_threshold = NULL, outpath, filename) {
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


# Function to perform PCA and extract data
prepare_pca_data <- function(vsd_data, intgroup) {
  pca_data <- plotPCA(vsd_data, intgroup = intgroup, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  list(pca_data = pca_data, percent_var = percent_var)
}

# Function to generate PCA plot
generate_pca_plot <- function(pca_data, percent_var, colour = "Affected", fill = NULL, shape = "Batch") {
  ggplot(pca_data, aes(x = PC1, y = PC2)) +
    # Map colour, shape, and fill
    geom_point(aes_string(colour = colour, shape = shape, fill = fill), size = 4, stroke = 1.5) +
    # Define scales for shape and fill
    scale_shape_manual(values = c(21, 22, 24)) +  # Circle, square, triangle
    scale_fill_manual(values = c("white", "darkgray"), guide = guide_legend(title = fill)) +
    # Define axis labels
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    coord_fixed() +
    # Custom guides to ensure proper legends
    guides(
      shape = guide_legend(override.aes = list(fill = NA)),  # Batch shapes, no fill
      colour = guide_legend(title = colour),  # Affected legend
      fill = guide_legend(override.aes = list(shape = 21, size = 4))  # Fill legend with consistent shape
    ) +
    theme_minimal()
}

# Function to generate a heatmap of the top N variable genes
generate_top_variable_genes_heatmap <- function(vsd_data, gene_info, annotation_columns, annotation_colors, top_n = 50) {
  # Identify top N variable genes
  top_var_genes <- head(order(-rowVars(assay(vsd_data))), top_n)

  # Subset and center the data matrix
  mat <- assay(vsd_data)[top_var_genes, ]
  mat <- mat - rowMeans(mat)

  # Create annotation dataframe
  annotation_df <- as.data.frame(colData(vsd_data)[, annotation_columns])

  # Map Ensembl IDs to gene names
  ensembl_to_gene <- setNames(gene_info$gene_name, gene_info$Ensembl_ID)
  rownames(mat) <- ensembl_to_gene[rownames(mat)]  # Set gene names as row names

  # Generate the heatmap
  ComplexHeatmap::pheatmap(
    mat = mat,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    fontsize = 5,
    angle_col = c("0")
  )
}

# Function to generate a heatmap for significant genes
generate_significant_genes_heatmap <- function(vsd_data, res_data, gene_info, annotation_columns, annotation_colors, output_file) {
  # Extract Ensembl IDs for significant genes
  top_genes_ensembl <- rownames(res_data)

  # Subset and center the data matrix
  significant_genes_matrix <- assay(vsd_data)[top_genes_ensembl, ]
  significant_genes_matrix <- significant_genes_matrix - rowMeans(significant_genes_matrix)

  # Map Ensembl IDs to gene names
  ensembl_to_gene <- setNames(gene_info$gene_name, gene_info$Ensembl_ID)
  rownames(significant_genes_matrix) <- ensembl_to_gene[rownames(significant_genes_matrix)]  # Replace Ensembl IDs with gene names

  # Order rows alphabetically by gene name
  significant_genes_matrix <- significant_genes_matrix[order(rownames(significant_genes_matrix)), ]

  # Create annotation dataframe
  annotation_df <- colData(vsd_data) %>%
    as.data.frame() %>%
    dplyr::select(all_of(annotation_columns))

  # Generate the heatmap and save to file
  png(output_file, width = 12, height = 10, units = "in", res = 1200)
  draw(
    ComplexHeatmap::pheatmap(
      mat = significant_genes_matrix,
      color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(50),
      annotation_col = annotation_df,
      annotation_colors = annotation_colors,
      fontsize = 6,
      angle_col = c("0"),
      cellheight = 10,
      cellwidth = 24,
      border_color = "black",
      heatmap_legend_param = list(title = "Matrix")
    )
  )
  dev.off()
}

# Function to generate and save a volcano plot
generate_volcano_plot <- function(res_data, gene_labels, x_col, y_col, select_genes,
                                  xlab_text, ylab_text, p_cutoff = 0.05, fc_cutoff = 1.0,
                                  xlim_range = c(-25, 25), ylim_range = c(0, 7),
                                  output_file = "volcano_plot.png") {
  # Create the EnhancedVolcano plot
  volcano_plot <- EnhancedVolcano(res_data,
                                  lab = res_data[[gene_labels]],
                                  x = x_col,
                                  y = y_col,
                                  selectLab = select_genes,
                                  xlab = xlab_text,
                                  ylab = ylab_text,
                                  pCutoff = p_cutoff,
                                  FCcutoff = fc_cutoff,
                                  pointSize = 2.0,
                                  labSize = 4.0,
                                  labCol = 'black',
                                  labFace = 'bold',
                                  boxedLabels = TRUE,
                                  colAlpha = 4 / 5,
                                  legendPosition = 'bottom',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 1.0,
                                  colConnectors = 'black',
                                  legendLabels = c(
                                    " NS",
                                    bquote(~Log[2]~ 'fold change'),
                                    " padj",
                                    bquote(~Log[2]~ 'fold change and padj')
                                  ),
                                  title = 'Affected versus Unaffected',
                                  subtitle = "Differential Gene Expression Analysis in ME/CFS Patients")

  # Add custom scaling
  volcano_plot <- volcano_plot +
    ggplot2::coord_cartesian(xlim = xlim_range) +
    ggplot2::scale_x_continuous(breaks = seq(xlim_range[1], xlim_range[2], 5)) +
    ggplot2::scale_y_continuous(breaks = seq(ylim_range[1], ylim_range[2], 1))

  # Save the plot
  ggsave(plot = volcano_plot, filename = output_file, device = "png", width = 14, dpi = 1200, units = "in")
}

# Function to map data to STRING, add color coding, and plot the network
generate_string_network <- function(data, gene_col, lfc_col, species = 9606, score_threshold = 100,
                                    version = "11.5", network_type = "full", top_hits = 29,
                                    input_directory = "") {
  # Initialize the STRINGdb object
  string_db <- STRINGdb$new(
    version = version,
    species = species,
    score_threshold = score_threshold,
    network_type = network_type,
    input_directory = input_directory
  )

  # Map the input data to STRING
  mapped_data <- string_db$map(data, gene_col, removeUnmappedRows = TRUE)

  # Get the top hits
  hits <- mapped_data$STRING_id[1:top_hits]

  # Add differential expression color coding
  mapped_data_colored <- string_db$add_diff_exp_color(mapped_data, logFcColStr = lfc_col)

  # Post the payload for network visualization
  payload_id <- string_db$post_payload(mapped_data_colored$STRING_id, colors = mapped_data_colored$color)

  # Plot the STRING network with the payload
  string_db$plot_network(hits, payload_id = payload_id)
}

# Utility function to generate an MA plot
generate_ma_plot <- function(data, genenames, main_title = "MA Plot",
                             output_file = NULL, fdr = 0.05, fc = 1, size = 0.4,
                             top_genes = 30, theme = ggplot2::theme_minimal()) {

  # Create the MA plot
  plot <- ggmaplot(
    data = data,
    main = main_title,
    fdr = fdr,
    fc = fc,
    size = size,
    genenames = genenames,  # Dynamically passed gene names
    ggtheme = theme,
    legend = "top",
    top = top_genes,
    font.label = c("bold", 6),
    label.rectangle = TRUE,
    font.legend = "bold",
    font.main = "bold"
  )

  # Save plot if output_file is specified
  if (!is.null(output_file)) {
    ggsave(plot = plot, filename = output_file, device = "png", dpi = 300, width = 10, height = 8, units = "in")
  }

  return(plot)
}

