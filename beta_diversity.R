# Script to analyze 16S data - beta diversity 
# Santiago et al., 2026 


# Load data ----
library(tidyverse)
library(ggplot2)
library(qiime2R)  


# calculating all the PCs available based on qiime file
file_path <- "weighted_ordination.txt"

# Import metadata
metadata <- read_q2metadata("q2-metadata.txt")
head(metadata)
metadata <- metadata %>%
  mutate(birth = recode(birth,
                        "VDCF" = "VD-CF",
                        "CSCF" = "CS-CF"))


# Define time points and colors
time_points <- c("25", "53", "90")
group.colors <- c("#6C9A8B","#E8998D", "#7E737350") 


pcoa_data <- read.delim(file_path, header = FALSE, sep = "\t", skip = 7, fill = TRUE)
head(pcoa_data)
pcoa_data_clean <- pcoa_data[, colSums(is.na(pcoa_data)) < nrow(pcoa_data)]
colnames(pcoa_data_clean) <- c("Site", paste("PC", 1:(ncol(pcoa_data_clean)-1), sep=""))
head(pcoa_data_clean)


pcoa_data_clean <- pcoa_data[!(pcoa_data$V1 %in% c("Site", "Biplot", "Site constraints")),]
colnames(pcoa_data_clean) <- c("Site", paste("PC", 1:(ncol(pcoa_data_clean)-1), sep=""))
head(pcoa_data_clean)



# Calculating PCoA based on pcoa_data_clean
calculate_global_pcoa <- function(data) {
  full_distance_matrix <- dist(data[, c("PC1", "PC2", "PC3")]) 
  pcoa_result <- cmdscale(full_distance_matrix, eig = TRUE, k = 3)
  eigenvalues <- pcoa_result$eig
  total_variance <- sum(eigenvalues[eigenvalues > 0])
  percent_variance <- (eigenvalues / total_variance) * 100
  list(
    pcoa_result = pcoa_result,
    percent_variance = percent_variance
  )
}
global_pcoa <- calculate_global_pcoa(pcoa_data_clean)


metadata$PC1 <- global_pcoa$pcoa_result$points[, 1]
metadata$PC2 <- global_pcoa$pcoa_result$points[, 2]
metadata$PC3 <- global_pcoa$pcoa_result$points[, 3]
global_percent_variance <- global_pcoa$percent_variance



plot_pcoa_for_time_point <- function(data, time_point_filter = NULL) {
  library(dplyr)
  library(ggplot2)
  
  if (!is.null(time_point_filter)) {
    data <- data %>% filter(time_point == time_point_filter)
  }
  if (!all(c("PC1", "PC2", "birth") %in% colnames(data))) {
    stop("Columns PC1, PC2 or birth not found")
  }
  
  hulls <- data %>%
    group_by(birth) %>%
    filter(n() >= 3) %>%
    slice(chull(PC1, PC2)) %>%
    ungroup()
  
  # ----- plot -----
  plot <- ggplot(data, aes(x = PC1, y = PC2, color = birth)) +
    geom_polygon(
      data = hulls,
      aes(x = PC1, y = PC2, fill = birth, group = birth),
      inherit.aes = FALSE,
      alpha = 0.25,
      color = NA
    ) +
    geom_point(size = 6, shape = 16, alpha = 0.95) +
    geom_hline(yintercept = 0, color = "gray40", size = 0.5) +
    geom_vline(xintercept = 0, color = "gray40", size = 0.5) +
    labs(
      title = ifelse(is.null(time_point_filter),
                     "Weighted UniFrac for all TP",
                     paste("PND", time_point_filter)),
      x = paste0("PC1 (", round(global_percent_variance[1], 2), "%)"),
      y = paste0("PC2 (", round(global_percent_variance[2], 2), "%)")
    ) +
    scale_color_manual(values = group.colors, name = "Birth") +
    scale_fill_manual(values = group.colors, guide = "none") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 22),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      plot.title = element_text(size = 22, hjust = 0.5)
    )
  
  return(plot)
}



time_points <- unique(metadata$time_point)
for (time_point in time_points) {
  time_point_plot <- plot_pcoa_for_time_point(metadata, time_point_filter = time_point)
  print(time_point_plot)
  ggsave(paste0("PCoA_WU_PND90_", time_point, ".pdf"),
         plot = time_point_plot, width = 9, height = 4)
}
