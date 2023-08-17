library(biomaRt)

check_order <- function(sample_metadata, counts) {
  # Check if data is ordered properly
  if (!all(rownames(sample_metadata) %in% colnames(counts)) ||
    !all(rownames(sample_metadata) == colnames(counts))) {
    # Order sample_metadata
    sample_metadata <- arrange(sample_metadata, RNA_Samples_id)

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
