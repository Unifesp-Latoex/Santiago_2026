#Script to analyze 16S data - relative abundances 
# Santiago et al., 2026

# load libraries ----
library(biomehorizon)
library(tidyverse)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(magrittr)
library(patchwork)
library(cowplot)
library(tibble)
library(ggembl)


# import data ----
# import ASVs
dada_res <- read_qza("filtered-table.qza")
names(dada_res)
dada_res$data[1:5,1:5]

# import metadata
metadata <- read_q2metadata("/q2-metadata.txt")
head(metadata)

# import taxonomy
taxonomy<-read_qza("taxonomy.qza")
head(taxonomy$data)
taxonomy_wide<-parse_taxonomy(taxonomy$data)
head(taxonomy_wide)

renamed_dada <- dada_res$data
genus_sums <- summarize_taxa(renamed_dada, taxonomy_wide)$Genus



# Normalize the abundances
genus_sums <- as.data.frame(genus_sums) %>%
  rownames_to_column(var = "ASVs")
abund_long <- genus_sums %>%
  pivot_longer(-ASVs, names_to = "SampleID", values_to = "Abundance")
abund_long_normalized <- abund_long %>%
  left_join(metadata %>% dplyr::select(SampleID, `time_point `, id, birth), by = "SampleID") %>%
  group_by(`time_point `, SampleID) %>%
  mutate(RelAbundance = Abundance / sum(Abundance)) %>%
  ungroup()
abund_long_normalized <- as.data.frame(abund_long_normalized)
metadata <- as.data.frame(metadata)

abund_long_normalized <- abund_long_normalized %>%
  separate(ASVs, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";", remove = FALSE, fill = "right", extra = "merge") %>%
  mutate(Genus = trimws(Genus))  # Remover espaços em branco

top_genera_df <- abund_long_normalized %>%
  dplyr::group_by(time_point, Genus) %>%
  dplyr::summarise(MeanAbundance = mean(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(time_point) %>%
  dplyr::slice_max(order_by = MeanAbundance, n = 10) %>%
  dplyr::ungroup()


abund_classified <- abund_long_normalized %>%
  dplyr::left_join(
    top_genera_df %>%
      dplyr::select(time_point, Genus) %>%
      dplyr::mutate(InTop = TRUE),
    by = c("time_point", "Genus")
  ) %>%
  dplyr::mutate(Genus = ifelse(is.na(InTop), "Other", Genus)) %>%
  dplyr::select(-InTop)


abund_classified <- abund_classified %>%
  group_by(id, time_point) %>%
  mutate(
    Genus = fct_reorder2(Genus, RelAbundance, interaction(id, time_point), .fun = sum)
  ) %>%
  ungroup()

abund_classified_plot <- abund_classified %>%
  group_by(id, time_point, birth, Genus) %>%  
  summarise(RelAbundance = sum(RelAbundance, na.rm = TRUE), .groups = "drop")

custom_colors <- c(
  "Lactobacillus" = "#6C9A8B", 
  "Bacteroides_H" = "#A0CBAE", 
  "Duncaniella" = "#B7D3C5", 
  "Limosilactobacillus" = "#E8A8A0", 
  "CAG-485" = "#F3C3B0",
  "Romboutsia_B" = "#F9D6C8",
  "Turicibacter" = "#867875",
  "Prevotella" = "#E2B8A8", 
  "d__Bacteria; Firmicutes_A; Clostridia_258483; Lachnospirales; Lachnospiraceae; NA" = "#B7A29D",
  "Phocaeicola_A_858004" = "#F5D0A9",
  "Ligilactobacillus" = "#DCE3D7", 
  "Copromonas" = "#9FB3A9", 
  "SFMI01" = "#BFD3C4",
  "UBA7173" = "#FFB34790",      
  "d__Bacteria; Bacteroidota; Bacteroidia; Bacteroidales; Muribaculaceae; NA" = "#FF9F9F", 
  "Alloprevotella" = "#F5D0A990",  
  "Paramuribaculum" = "#D8E4E0", 
  "Other" = "#A6948D" 
)




timepoint_labeller <- function(values) {
  paste0("PND ", values)
}

for (grupo in unique(abund_classified_plot$birth)) {
  
  plot_group <- abund_classified_plot %>%
    filter(birth == grupo) %>%
    ggplot(aes(x = id, y = RelAbundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = custom_colors) +  # sua paleta
    facet_wrap(~ time_point, 
               scales = "free_x", 
               labeller = labeller(time_point = timepoint_labeller),
               strip.position = "bottom") +
    theme_publication() +  # tema ggembl
    theme(
      axis.text.x = element_blank(),        # remove labels do eixo x
      axis.ticks.x = element_blank(),       # remove ticks do eixo x
      axis.text.y = element_text(size = 18),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 18),
      strip.text = element_text(size = 18),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 15),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5)  # centralizado
    ) +
    labs(
      title = paste(grupo),
      y = "Relative abundance",
      fill = "Genus"
    )
  
  print(plot_group)
  
  # salvar pdf
  ggsave(
    paste0("paper_rel_", grupo, ".pdf"),
    plot = plot_group,
    height = 8, width = 12, device = "pdf"
  )
}




# to know how many genera are shared
calcular_venn_pnd <- function(df_pnd, key) {
  l_vd    <- df_pnd %>% filter(birth == "VD", RelAbundance > 0) %>% pull(Genus) %>% unique()
  l_vd_cf <- df_pnd %>% filter(birth == "VD-CF", RelAbundance > 0) %>% pull(Genus) %>% unique()
  l_cs_cf <- df_pnd %>% filter(birth == "CS-CF", RelAbundance > 0) %>% pull(Genus) %>% unique()
  
  total_pnd <- length(unique(c(l_vd, l_vd_cf, l_cs_cf)))
  compartilhados_todos <- length(intersect(l_vd, intersect(l_vd_cf, l_cs_cf)))
  return(data.frame(
    Total_Detectado = total_pnd,
    Compartilhado_3_Grupos = compartilhados_todos
  ))
}


resumo_venn_final <- abund_long_normalized %>%
  group_by(time_point) %>%
  group_modify(~ calcular_venn_pnd(.x)) %>%
  ungroup()
print(resumo_venn_final)



# to extract the genus
extrair_listas_generos <- function(df_pnd, key) {
  
  l_vd    <- df_pnd %>% filter(birth == "VD", RelAbundance > 0) %>% pull(Genus) %>% unique()
  l_vd_cf <- df_pnd %>% filter(birth == "VD-CF", RelAbundance > 0) %>% pull(Genus) %>% unique()
  l_cs_cf <- df_pnd %>% filter(birth == "CS-CF", RelAbundance > 0) %>% pull(Genus) %>% unique()

  todos_3   <- intersect(l_vd, intersect(l_vd_cf, l_cs_cf))
  so_vd     <- setdiff(l_vd, union(l_vd_cf, l_cs_cf))
  so_vd_cf  <- setdiff(l_vd_cf, union(l_vd, l_cs_cf))
  so_cs_cf  <- setdiff(l_cs_cf, union(l_vd, l_vd_cf))
  
  res <- bind_rows(
    data.frame(Genus = todos_3,  Category = "Shared_by_all_3"),
    data.frame(Genus = so_vd,    Category = "Exclusive_to_VD"),
    data.frame(Genus = so_vd_cf, Category = "Exclusive_to_VD-CF"),
    data.frame(Genus = so_cs_cf, Category = "Exclusive_to_CS-CF")
  )
  
  return(res)
}

# generate complete list
lista_generos_suplementar <- abund_long_normalized %>%
  group_by(time_point) %>%
  group_modify(~ extrair_listas_generos(.x)) %>%
  ungroup() %>%
  rename(PND = time_point)


wb <- createWorkbook()
addWorksheet(wb, "Summary_Numbers")
writeData(wb, "Summary_Numbers", tabela_suplementar_venn) 

addWorksheet(wb, "Full_Genus_Lists")
writeData(wb, "Full_Genus_Lists", lista_generos_suplementar)

saveWorkbook(wb, 
             file = "Supplementary_Table_Complete_Biodiversity.xlsx", 
             overwrite = TRUE)

