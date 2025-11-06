# Load the required library
library(tidyverse)

# Read the CSV file
data <- read.csv("data/MECFS_RNAseq_metadata_2023_08_16.csv")

# Update the "Batch" column using case_when
data <- data %>%
  mutate(Batch = case_when(
    Batch == "A" ~ "1",
    Batch == "B" ~ "2",
    Batch == "C" ~ "3",
    TRUE ~ Batch
  ))



# Write the updated data back to the CSV file
write.csv(data, "data/MECFS_RNAseq_metadata_2023_08_18.csv", row.names = FALSE)
