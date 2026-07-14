# Script to analyze 16S data  - GLM
# Santiago et al., 2026

# --- Import data ----
library(qiime2R)
library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library(openxlsx)
library(ggbreak)
library(openxlsx)
library(writexl)
library(tidyr)
library(emmeans)
library(openxlsx)

# ASV table
dada_res <- read_qza("filtered-table.qza")
asv_table <- dada_res$data

# Metadata
metadata <- read_q2metadata("q2-metadata.txt")
metadata <- metadata %>% rename(time_point = 'time_point ')
metadata <- metadata %>%
  mutate(birth = recode(birth,
                        "VDCF" = "VD-CF",
                        "CSCF" = "CS-CF"))
metadata$birth <- factor(metadata$birth, levels = c("VD", "VD-CF", "CS-CF"))

# Taxonomy
taxonomy <- read_qza("taxonomy.qza")
taxonomy_wide <- parse_taxonomy(taxonomy$data)

# --- Prepare top genera ----
genus_sums <- summarize_taxa(asv_table, taxonomy_wide)$Genus
genus_sums_rel <- prop.table(as.matrix(genus_sums), margin = 2)

# Filter samples by time points
time_points_to_keep <- c("25", "53", "90")
metadata_filtered <- metadata %>% filter(time_point %in% time_points_to_keep)
genus_sums_rel <- genus_sums_rel[, colnames(genus_sums_rel) %in% metadata_filtered$SampleID]

# Determine top 10 genera per time point and combine to top 30
top_10_per_time <- lapply(time_points_to_keep, function(tp) {
  samples_tp <- metadata_filtered %>% filter(time_point == tp) %>% pull(SampleID)
  top_genus <- sort(rowSums(genus_sums_rel[, colnames(genus_sums_rel) %in% samples_tp]), decreasing = TRUE)[1:10]
  names(top_genus)
})
top_30_genus <- unique(unlist(top_10_per_time))

# Convert to long format
genus_df <- as.data.frame(genus_sums_rel[top_30_genus, ])
genus_df$Genus <- rownames(genus_df)
genus_sums_filtered_long <- genus_df %>%
  pivot_longer(-Genus, names_to = "SampleID", values_to = "Abundance") %>%
  left_join(metadata_filtered, by = "SampleID")

# Shortened genus names
extract_genus <- function(x) sapply(strsplit(x, "; "), function(y) tail(y, 1))
genus_name_map <- setNames(extract_genus(top_30_genus), top_30_genus)

# --- Top 30 plotting with statistics ----
top_30_genus_names <- top_30_genus
statistical_results_list <- list()

for (genus_name in top_30_genus_names) {
  
  genus_data <- genus_sums_filtered_long %>%
    filter(Genus == genus_name, !is.na(time_point))
  
  genus_data$birth <- factor(genus_data$birth, levels = c("VD", "VD-CF", "CS-CF"))
  
  genus_summary <- genus_data %>%
    group_by(time_point, birth) %>%
    summarise(
      mean_abundance = mean(Abundance, na.rm = TRUE),
      sem_abundance = sd(Abundance, na.rm = TRUE)/sqrt(n()),
      .groups = "drop"
    )
  
  genus_short_name <- genus_name_map[genus_name]
  
  p <- ggplot(genus_summary, aes(x = time_point, y = mean_abundance, fill = birth)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean_abundance - sem_abundance, ymax = mean_abundance + sem_abundance),
                  width = 0.2, position = position_dodge(width = 0.8)) +
    geom_jitter(data = genus_data,
                aes(x = time_point, y = Abundance, fill = birth),  
                color = "gray30",                                
                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                size = 2, alpha = 0.8, inherit.aes = FALSE) +
    scale_fill_manual(values = c("VD"="#7E737350", "VD-CF"="#E8998D", "CS-CF"="#6C9A8B")) +
    labs(title = genus_short_name, 
         y = "Relative Abundance", 
         x = "PND") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 23),   
      axis.text.y = element_text(size = 23),                           
      axis.title.x = element_text(size = 26),      
      axis.title.y = element_text(size = 26),           
      axis.line.x = element_line(size = 1.5),
      axis.line.y = element_line(size = 1.5),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 20),                         
      legend.title = element_text(size = 22)         
    )
  
  ggsave(paste0(genus_short_name,"_abundance.pdf"),
         plot=p, height=6, width=8)
}



# glm 

statistical_results_list <- list()

for (i in seq_along(top_30_genus)) {
  genus_name <- top_30_genus[i]
  genus_short_name <- genus_name_map[genus_name]
  
  genus_data <- genus_sums_filtered_long %>%
    filter(Genus == genus_name, !is.na(time_point))
  
  if (nrow(genus_data) == 0) next
  
  genus_data$birth <- factor(genus_data$birth, levels = c("VD", "VD-CF", "CS-CF"))
  
  
  glm_results <- glm(Abundance ~ birth * time_point, data = genus_data, family = gaussian())
  
  
  emm_results <- emmeans(glm_results, pairwise ~ birth | time_point, infer = c(TRUE, TRUE))
  emm_contrasts <- as.data.frame(emm_results$contrasts)
  
  
  if(!all(c("lower.CL","upper.CL") %in% colnames(emm_contrasts))) {
    emm_contrasts$lower.CL <- NA
    emm_contrasts$upper.CL <- NA
  }
  
  
  ancom_like <- emm_contrasts %>%
    select(time_point, contrast, estimate, SE, df, lower.CL, upper.CL, t.ratio, p.value) %>%
    rename(logFC = estimate,
           SE = SE,
           df = df,
           lower.CL = lower.CL,
           upper.CL = upper.CL,
           t = t.ratio,
           p = p.value) %>%
    mutate(
      p_adj = p.adjust(p, method = "fdr"),
      genus = genus_short_name
    ) %>%
    select(genus, time_point, contrast, logFC, SE, df, lower.CL, upper.CL, t, p, p_adj)
  
  statistical_results_list[[genus_short_name]] <- ancom_like
}


stat_results_df <- bind_rows(statistical_results_list)

write.xlsx(stat_results_df,
           file = "top30_genus.xlsx",
           rowNames = FALSE)


