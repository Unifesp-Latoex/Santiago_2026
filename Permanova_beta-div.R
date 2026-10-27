# load libraries
library(Polychrome)
library(RColorBrewer)
library(vegan)
library(ggplot2)
library(qiime2R)
library(dplyr)



# import metadata
metadata <- read_q2metadata("q2-metadata.txt")
head(metadata)


uwu_df <- read.csv(
  "uwunifrac-distance-matrix.csv", 
  header = TRUE, 
  row.names = 1, 
  check.names = FALSE
)
u_dist_matrix <- as.dist(as.matrix(uwu_df))

wu_df <- read.csv(
  "wunifrac-distance-matrix.csv", 
  header = TRUE, 
  row.names = 1, 
  check.names = FALSE
)
wu_dist_matrix <- as.dist(as.matrix(wu_df))

metadata <- metadata %>%
  mutate(
    bt = case_when(
      grepl("CS-CF", birth, ignore.case = TRUE) ~ "CS",
      grepl("VD|VD-CF", birth, ignore.case = TRUE) ~ "VD",
      TRUE ~ NA_character_ 
    )
  )



# Select only the relevant columns for metadata
metadata <- metadata %>% select(SampleID, time_point, birth, id)

# Ensure that metadata_cols is a character vector of column names
metadata_cols <- colnames(metadata)[2:length(metadata)]  # Assuming metadata columns start at index 4


common_samples <- intersect(rownames(as.matrix(u_dist_matrix)), metadata$SampleID)

u_dist_sub <- as.dist(as.matrix(u_dist_matrix)[common_samples, common_samples])

metadata_sub <- metadata %>% 
  filter(SampleID %in% common_samples) %>% 
  arrange(match(SampleID, common_samples))

adonis_interaction <- adonis2(
  u_dist_sub ~ birth * time_point, 
  data = metadata_sub, 
  strata = metadata$id,
  permutations = 999, 
  by = NULL
)
print(adonis_interaction)



