# >>>>>>>>>>>>>>>> FIGURE 3A <<<<<<<<<<<<<<<<<<<
library(readxl)
CountsKD <-  read_excel("/media/4tb2/Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig3a-MicroorgPerLevel.xlsx")
library(tidyverse)
CountsKD_l <- CountsKD %>%
  pivot_longer(cols = c(Phylum, Genus, Species),
               names_to = "Taxonomic_Levels",
               values_to = "Counts")

CountsKD_l$Counts <- as.numeric(CountsKD_l$Counts)
order_levels <- colnames(CountsKD)[2:10]
CountsKD_l$Taxonomic_Levels <- factor(CountsKD_l$Taxonomic_Levels, levels = order_levels)
unique(CountsKD_l$Methodology)
CountsKD_l$Methodology <- factor(CountsKD_l$Methodology, levels = c("-D", "DdhD", "BoD", "RsD", "bwaD",
                                                                              "-K", "BoK", "RsK", "bwaK"))
CountsKD_l <- CountsKD_l %>%
  mutate(GroupColor = ifelse(grepl("D", Methodology), "DRAGEN", "Kraken"))

#Color:
colores_metodologia <- c("Kraken" = "#FF7F00",
                         "DRAGEN" = "#9370DB")

ggplot(CountsKD_l, aes(x = Methodology, y = Counts, fill = GroupColor)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +  # Boxplots agrupados
  geom_line(aes(group = interaction(ID, Taxonomic_Levels)),
            colour = "black", linetype = "11", size = 0.3) +  # Líneas finas conectando puntos

  theme_minimal() +
  labs(x = "Taxonomic Level", y = "Number of Microorganisms",
       title = "") +

  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 15)),  # Separar título eje x de los valores
    axis.title.y = element_text(margin = margin(r = 15)),  # Separar título eje y de los valores
    legend.position = "bottom",  # **Ubica la leyenda abajo**
    legend.box = "horizontal"  # **Hace que los elementos de la leyenda se distribuyan horizontalmente**
  ) +

  scale_fill_manual(values = colores_metodologia, name = "Taxonomic Classifier") +
  facet_wrap(~ Taxonomic_Levels, scales = "free_y")  # Crear gráficos por cada nivel taxonómico con escalas independientes


# NEW FIGURE DUE TO REVIEWERS' COMMENTS -----------------------------
library(dplyr)

ddhd_vals <- CountsKD_l %>%
  filter(Methodology == "DdhD") %>%
  select(ID, Taxonomic_Levels, DdhD_Counts = Counts)

diff_df <- CountsKD_l %>%
  filter(Methodology != "DdhD") %>%
  left_join(ddhd_vals, by = c("ID", "Taxonomic_Levels")) %>%
  mutate(Difference = 100 * (Counts - DdhD_Counts) / DdhD_Counts)


ggplot(diff_df, aes(x = Methodology, y = Difference, fill = GroupColor)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  
  #geom_line(aes(group = interaction(ID, Taxonomic_Levels)),
  #          colour = "black", linetype = "11", size = 0.3) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  
  theme_minimal() +
  labs(x = "Taxonomic Classifier",
       #y = "Difference  Nº Microorganism \n(Method - DdhD)",
       y = "Relative Difference (%) \n(Method - DdhD)",
       title = "") +
  
  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  
  scale_fill_manual(values = colores_metodologia, name = "") +
  facet_wrap(~ Taxonomic_Levels, scales = "free_y")

library(dplyr)


wilcox_results <- diff_df %>%
  group_by(Methodology, Taxonomic_Levels) %>%
  summarise(
    p_value = wilcox.test(Counts, DdhD_Counts, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_label = paste0("p = ", signif(p_value, 3)))




#------------------------------------------



#Patterns:
library(ggpattern)
patterns <- c("DRAGEN" = "stripe", "KRAKEN" = "none")

ggplot(CountsKD_l, aes(x = Methodology, y = Counts, pattern = GroupColor)) +
  geom_boxplot_pattern(
    aes(group = interaction(Methodology, GroupColor)),
    pattern_fill = "black",
    pattern_density = 0.4,
    pattern_spacing = 0.05,
    pattern_angle = 45,
    fill = "grey70",
    color = "black",
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  geom_line(
    aes(group = interaction(ID, Taxonomic_Levels)),
    colour = "black", linetype = "11", size = 0.2
  ) +

  theme_minimal() +
  labs(x = "Taxonomic Level", y = "Number of Microorganisms",
       title = "") +

  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +

  scale_pattern_manual(values = patterns, name = "Taxonomic Classifier") +
  facet_wrap(~ Taxonomic_Levels, scales = "free_y")


#>>>>>>>>>>>>>>>>>>>> FIGURE 3B <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

load("/media/4tb2/Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig3-data_list.RData")
library(UpSetR)
library(ComplexHeatmap) 
data_matrix <- fromList(data_list)

#>>>>>>>>>>>>>>>> FIGURE 3B <<<<<<<<<<<<<<<<<<

bar_colors <- rep("grey15", 25)
bar_colors[c(2,3,9)] <- c( "#FF7F00", "#9370DB", "#56B4E9")  
upset_plot <- upset(
  data_matrix,
  sets = names(data_list),
  order.by = "freq",
  mainbar.y.label = "Number of Species",
  #mainbar.y.label = "Number of Genus",
  sets.x.label = "Nº Detected Species",
  text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.3),  
  point.size = 2,  
  line.size = 0,  
  keep.order = TRUE,  
  sets.bar.color = "grey50",  
  matrix.color = "black",  
  main.bar.color = bar_colors,
  shade.color = "gray90", 
  number.angles = 0
)

upset_plot


##########################################################

#Process one sample with reporting zero counts -----------------------------
args <- c(
  sprintf("--db %s", db_kraken2),
  sprintf("--confidence %s", confidence),
  "--threads 15",
  "--use-names",
  "--report-zero-counts",
  sprintf("--output %s/Resultados_KRAKEN/output_%s.txt", patient_dir, de_host_file),
  sprintf("--report %s/Resultados_KRAKEN/report_%s.sequences_CON_CEROS", patient_dir, de_host_file),
  sprintf("--paired %s %s", fileR1, fileR2)
)

#system2(path_kraken2, args = args)

#Here we can see all the taxa included in the Kraken's  Data Base 

species <- unique(unlist(data_list))
df_presence <- as.data.frame(data_matrix)
df_presence$Specie <- species  
df_presence$Intersec <- apply(df_presence[,-ncol(df_presence)], 1, paste, collapse = "-")

# Agrupar géneros por intersección
intersec <- split(df_presence$Specie, df_presence$Intersec)
colnames(df_presence)
e_enD_noK <- df_presence$Specie[df_presence$Intersec == "0-0-0-0-1-1-1-1-1"]

All_Taxa_Kraken <- read_excel("/media/4tb2/Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig3-All_Taxa_Kraken.xlsx")
all(e_enD_noK %in% All_Taxa_Kraken$V6)
any(e_enD_noK %in% All_Taxa_Kraken$V6)

length(e_enD_noK)
no_DB_Kraken <- e_enD_noK[ -which( e_enD_noK %in% All_Taxa_Kraken$V6)]
length(no_DB_Kraken)/length(e_enD_noK) *100
# 90%




