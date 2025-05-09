# >>>>>>>>>>>>>>>> FIGURE 3A <<<<<<<<<<<<<<<<<<<
CountsKD <-  read_excel("Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig3a-MicroorgPerLevel.xlsx")

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

#  p-values - Wilcoxon by pares
calculate_pval <- function(data) {
  data <- data %>%
    arrange(ID)

  metodologies <- unique(data$Methodology)
  combinations <- combn(metodologies, 2, simplify = FALSE)

  par <- combinations[[3]]
  par <- combinations[[2]]
  resultados <- lapply(combinations, function(par) {

    subset_data <- data %>% filter(Methodology %in% par)

    spread_data <- subset_data %>%
      group_by(ID, Metodologia) %>%
      summarise(Conteos = first(Counts)) %>%
      ungroup() %>%
      pivot_wider(names_from = Methodology, values_from = Counts)

    spread_data$dif = spread_data[[par[1]]] - spread_data[[par[2]]]

    if(all(spread_data$dif == 0)) {
      return(message(sprintf("The differences in %s are all equal to 0", par)))
    }

    test <- wilcox.test(
      Counts ~ Methodology,
      data = subset_data,
      paired = TRUE
    )

    tibble(
      Taxonomic_Levels = unique(data$Taxonomic_Levels),
      Metodology1 = par[1],
      Metodology2 = par[2],
      p_value = test$p.value
    )
  })

  bind_rows(resultados)
}

# Aplicar la función a cada nivel taxonómico
pval_pares <- CountsKD_l %>%
  group_split(Taxonomic_Levels) %>%
  lapply(calculate_pval) %>%
  bind_rows()

print(pval_pares)

write.xlsx(pval_pares, file = "~/Daniela/Biota/PipelineBiota/paraPaper/Pvalores_Metodologias_NumberMicroog.xlsx")


#>>>>>>>>>>>>>>>>>>>> FIGURE 3B <<<<<<<<<<<<<<<<<<<<

library(UpSetR)
library(ComplexHeatmap)  # Necesario para fromList()
data_matrix <- fromList(data_list)

#>>>>>>>>>>>>>>>> FIGURE 3B <<<<<<<<<<<<<<<<<<
# Definir los colores: negro para la mayoría, colores solo en las primeras 4 barras
bar_colors <- rep("black", 25)
bar_colors[c(2,3,9)] <- "grey40"  # Colores para destacar

bar_colors <- rep("grey15", 25)
bar_colors[c(2,3,9)] <- c("#9370DB", "#FF7F00", "#56B4E9")  # Colores para destacar

upset_plot <- upset(
  data_matrix,
  sets = names(data_list),
  order.by = "freq",
  mainbar.y.label = "Number of Species",
  #mainbar.y.label = "Number of Genus",
  sets.x.label = "Nº Detected Species",
  text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.3),  # Aumenta la legibilidad del texto
  point.size = 2,  # Puntos más grandes en la matriz
  line.size = 0,  # Líneas más gruesas para mayor visibilidad
  keep.order = TRUE,  # Mantiene el orden original de los conjuntos
  sets.bar.color = "grey50",  # Todas las barras de conjuntos en negro
  matrix.color = "black",  # Matriz en negro
  main.bar.color = bar_colors,
  shade.color = "gray90",  # Fondo de sombreado más claro para contraste
  number.angles = 0,  # Mejor legibilidad de números
  #mainbar.y.max = max(rowSums(data_matrix)) * 0.8,  # Reduce la altura de la barra principal
  #mb.ratio = c(0.6, 0.4)  # Ajusta la proporción entre matriz y barra principal
)

upset_plot



