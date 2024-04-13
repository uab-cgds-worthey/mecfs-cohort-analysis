library(biomaRt)

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
