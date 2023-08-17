# Load the required library
library(tidyverse)

# Read the CSV file
data <- read.csv("data/MECFS_RNAseq_metadata_2023_08_16.csv")

# Define a mapping between numbers and letters
batch_mapping <- c("1" = "A", "2" = "B", "3" = "C")

# Update the "Batch" column using case_when
data <- data %>%
  mutate(Batch = case_when(
    Batch == "B1" ~ "A",
    Batch == "B2" ~ "B",
    Batch == "B3" ~ "C",
    TRUE ~ Batch
  ))



# Write the updated data back to the CSV file
write.csv(data, "data/MECFS_RNAseq_metadata_2023_08_17.csv", row.names = FALSE)
