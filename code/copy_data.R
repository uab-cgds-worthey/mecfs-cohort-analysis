# This file copies data from the necessary folders of the analysis
# and creates necessary data folders within this workflowr project.
# We will use the salmon.merged.gene_counts_length_scaled.tsv

library("fs")
analysis_dir <- "/data/project/worthey_lab/projects/MECFS_Ramsey/analysis/project_level_analysis/bulk_rnaseq"

counts_file <- "salmon.merged.gene_counts_length_scaled.tsv"

# Create subdirectories in workflowr data folder
dir_create(path = "data/star-salmon", recurse = TRUE)

# Copy files from analysis directory to workflowr project directory
# This line ONLY works on Cheaha
file_copy(path(analysis_dir, "star_salmon", counts_file), path("data", "star-salmon", counts_file))

# Create subdirectories in workflowr output folder
dir_create(path = "output/batch-correction-limma/plot-counts", recurse = TRUE)
