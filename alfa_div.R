#Script to analyze 16S data  - alfa diversity
# Santiago et al., 2026

# Load packages
library(ggplot2)
library(dplyr)
library(vegan)
library(openxlsx)
library(qiime2R)
library(biomehorizon)
library(tidyverse)
library(ggpubr)
library(magrittr)
library(patchwork)
library(cowplot)
library(lme4)
library(emmeans)
library(pairwiseAdonis)
library(knitr)

# Load data ----

# import metadata
metadata <- read_q2metadata("q2-metadata.txt")
head(metadata)

# import indexes
shannon<-read_qza("shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID")
gplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))
metadata <-
  metadata %>%
  left_join(shannon)

obs_species<-read_qza("observed_features_vector.qza")
obs_species <- obs_species$data %>% rownames_to_column("SampleID")
gplots::venn(list(metadata=metadata$SampleID, obs_species=obs_species$SampleID))
metadata <-
  metadata %>%
  left_join(obs_species)

evenness<-read_qza("evenness_vector.qza")
evenness <- evenness$data %>% rownames_to_column("SampleID")
gplots::venn(list(metadata=metadata$SampleID, evenness=evenness$SampleID))
metadata <-
  metadata %>%
  left_join(evenness)

faithpd<-read_qza("faith_pd_vector.qza")
faithpd <- faithpd$data %>% rownames_to_column("SampleID")
faithpd <- faithpd %>%
  rename(index = SampleID, SampleID = V1, faith = V2)
gplots::venn(list(metadata=metadata$SampleID, faithpd=faithpd$SampleID))
metadata <-
  metadata %>%
  left_join(faithpd)


# Set details
metadata <- metadata %>%
  mutate(birth = recode(birth,
                        "VDCF" = "VD-CF",
                        "CSCF" = "CS-CF"))

group.colors <- c('CS-CF' = "#6C9A8B", 'VD-CF' = "#E8998D", 'VD' ="#7E737350")
metadata$time_point <- as.factor(metadata$time_point)
metadata$birth <- as.factor(metadata$birth)
metadata$birth <- factor(metadata$birth, levels = c("VD", "VD-CF", "CS-CF"))
metadata$time_point <- factor(metadata$time_point, levels = c("25", "53", "90"))


metrics <- list(
  shannon = list(col = "shannon_entropy", label = "Shannon Diversity"),
  faith   = list(col = "faith", label = "Faith PD"),
  pielou  = list(col = "pielou_evenness", label = "Pielou Evenness"),
  observed = list(col = "observed_features", label = "Observed Features")
)


all_stats <- data.frame()

for (metrics_name in names(metrics)) {
  cat("Metrics:", metrics_name, "\n")
  
  metric_col <- metrics[[metrics_name]]$col
  metric_label <- metrics[[metrics_name]]$label
  if (!metric_col %in% colnames(metadata)) {
    stop(paste("Column", metric_col, "not found in metadata"))
  }

  df_metric <- metadata %>%
    dplyr::select(SampleID, birth, time_point, !!sym(metric_col)) %>%
    drop_na()
  dist_matrix <- dist(df_metric[[metric_col]])
  dispersion_result <- betadisper(dist_matrix, group = df_metric$birth)
  dispersion_test <- permutest(dispersion_result)
  
  # PERMANOVA
  permanova_result <- adonis2(dist_matrix ~ birth * time_point,
                              data = df_metric, permutations = 999)
  
  # Pairwise PERMANOVA
  pairwise_results <- data.frame()
  for (tp in unique(df_metric$time_point)) {
    df_tp <- df_metric %>% filter(time_point == tp)
    dist_tp <- dist(df_tp[[metric_col]])
    pw <- pairwise.adonis(dist_tp, factors = df_tp$birth)
    pw$time_point <- tp
    pairwise_results <- rbind(pairwise_results, pw)
  }
  
 
  sig_pairs <- pairwise_results %>%
    filter(p.adjusted < 0.05) %>%
    mutate(
      group1 = gsub(" vs .*", "", pairs),
      group2 = gsub(".* vs ", "", pairs),
      y.position = max(df_metric[[metric_col]], na.rm = TRUE) * 1.05,  
      label = case_when(
        p.adjusted < 0.001 ~ "***",
        p.adjusted < 0.01  ~ "**",
        p.adjusted < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
  
  
  stat_df <- data.frame(
    Metric = metric_label,
    PERMANOVA_p = permanova_result$`Pr(>F)`[1],
    pairwise_results
  )
  all_stats <- rbind(all_stats, stat_df)
  
 
  p <- ggplot(df_metric, aes(x = time_point, y = !!sym(metric_col), fill = birth)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA,
                 color = "black", linewidth = 0.8,
                 position = position_dodge(width = 0.8)) +
    geom_jitter(color = "gray30",
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                alpha = 0.8, size = 4) +
    scale_fill_manual(values = group.colors) +
    labs(x = "Postnatal day", y = metric_label, title = "") +
    theme_classic() +
    theme(
      legend.position = "bottom",
      text = element_text(size = 38),
      axis.text = element_text(size = 40),
      axis.title = element_text(size = 44),
      legend.text = element_text(size = 34),
      legend.title = element_text(size = 42),
      axis.line.x = element_line(size = 1.5),
      axis.line.y = element_line(size = 1.5)
    ) +
    annotate("text", x = 1.5, y = min(df_metric[[metric_col]], na.rm = TRUE) * 0.95,
             label = paste("PERMANOVA p =", signif(permanova_result$`Pr(>F)`[1], 3)),
             hjust = 0.5, size = 5, color = "black") +
    stat_pvalue_manual(sig_pairs, label = "label",
                       xmin = "group1", xmax = "group2",
                       y.position = "y.position", tip.length = 0.02)

  ggsave(paste0("plot_", metrics_name, ".pdf"), plot = p, height = 6, width = 8)
}

write.xlsx(all_stats, "stats_all_metrics.xlsx", rowNames = FALSE)

