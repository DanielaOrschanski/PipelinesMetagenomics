library(ggpattern)

# >>>>>>>>>>>>>>>> Figure 5 A <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
library(PipelineBiota)
library(readxl)
library(dplyr)
df_proteobacteria <- read_excel("Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig5-DFProteobacteria.xlsx")
df_proteobacteria$Rango <- factor(df_proteobacteria$`Rango etario`, levels = c("18-35", "35-55", ">55") )

kruskal_results <- df_proteobacteria %>%
  group_by(Methodology) %>%
  summarise(p_value = kruskal.test(RelFreq ~ Rango)$p.value) %>%
  mutate(p_text = paste0("p = ", signif(p_value, 2)))

df_proteobacteria <- df_proteobacteria %>%
  mutate(`Taxonomic Clasiffier` = ifelse(grepl("D", Methodology), "DRAGEN",
                                         ifelse(grepl("K", Methodology), "Kraken", "Other")))

unique(df_proteobacteria$Methodology)

df_proteobacteria$`Taxonomic Clasiffier` <- as.factor(df_proteobacteria$`Taxonomic Clasiffier`)
df_proteobacteria$Methodology <- factor(df_proteobacteria$Methodology,
                                        levels = c("-D", "BoD", "RsD", "bwaD","DdhD","-K", "BoK", "RsK", "bwaK"),
                                        ordered = TRUE)
levels(df_proteobacteria$Methodology)
colores <- c("DRAGEN" = "#8470FF" , "Kraken" = "#FFA500")

library(ggplot2)
ggplot(df_proteobacteria, aes(x = Rango, y = RelFreq, pattern = `Taxonomic Clasiffier`, fill = `Taxonomic Clasiffier`)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_jitter(aes(color = `Taxonomic Clasiffier`), width = 0.2, alpha = 0.5, size = 1) +
  geom_boxplot(width = 0.6, outlier.shape = NA, color = "black") +
  scale_fill_manual(values = colores) +
  scale_color_manual(values = colores) +
  theme_minimal() +
  labs(x = "Age Range", y = "Proteobacteria's Frequency", title = "") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
    #legend.box = "horizontal"
  ) +
  facet_wrap(~ Methodology, scales = "free", ncol = 5, strip.position = "top") +
  geom_text(data = kruskal_results,
            aes(x = 2, y = max(df_proteobacteria$RelFreq, na.rm = TRUE) + 0.05, label = p_text),
            inherit.aes = FALSE, size = 3)


# >>>>>>>>>>>> FIGURE S1 and 5b <<<<<<<<<<<<

#Fig 5b ---------------------------------------------------------------
df_sign_sex <- read_excel("Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig5-Genera_DfSign_Sex.xlsx")
data_list <- list(
  BoK = unique(df_sign_sex$Via[df_sign_sex$Wilcoxon_KBo < 0.05]),
  RsK = unique(df_sign_sex$Via[df_sign_sex$Wilcoxon_KRs < 0.05]),
  bwaK = unique(df_sign_sex$Via[df_sign_sex$Wilcoxon_KBWA < 0.05]),
  '-K' = unique(df_sign_sex$Via[df_sign_sex$Wilcoxon_K < 0.05]),

  BoD = unique(df_sign_sex$Via[df_sign_sex$Wilcoxon_DBo < 0.05]),
  RsD = unique(df_sign_sex$Via[df_sign_sex$Wilcoxon_DRs < 0.05]),
  bwaD = unique(df_sign_sex$Via[df_sign_sex$Wilcoxon_DBWA < 0.05]),
  '-D' = unique(df_sign_sex$Via[df_sign_sex$Wilcoxon_D < 0.05]),
  DdhD = unique(df_sign_sex$Via[df_sign_sex$Wilcoxon_DdhD< 0.05])
)

bar_colors <- rep("black", 45)
bar_colors[c(9, 10, 16)] <- "#56B4E9"
bar_colors[c(1,3)] <- c("#8470FF","#FFA500")
data_matrix <- fromList(data_list)

library(UpSetR)
upset(
  data_matrix,
  sets = names(data_list),
  order.by = "freq",
  mainbar.y.label = "Nº Genera signif dif between sexes",
  sets.x.label = "Nº Signif Genera Detected",
  text.scale = c(1.5, 1.3, 1.2, 1, 1.3, 1.3),  # Aumenta la legibilidad del texto
  point.size = 2,  # Puntos más grandes en la matriz
  line.size = 0,  # Líneas más gruesas para mayor visibilidad
  keep.order = TRUE,  # Mantiene el orden original de los conjuntos
  sets.bar.color = "grey50",  # Todas las barras de conjuntos en negro
  matrix.color = "black",  # Matriz en negro
  main.bar.color = bar_colors,
  shade.color = "gray90",  # Fondo de sombreado más claro para contraste
  number.angles = 0,  # Mejor legibilidad de números
  mb.ratio = c(0.7, 0.3)
)



# Fig S1: -------------------------------------------------------------------
df_sign_age <- read_excel("Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig5-Genera_DfSign_Age.xlsx")

data_list <- list(
  BoK = unique(df_sign_age$Via[df_sign_age$Wilcoxon_KBo < 0.05]),
  RsK = unique(df_sign_age$Via[df_sign_age$Wilcoxon_KRs < 0.05]),
  bwaK = unique(df_sign_age$Via[df_sign_age$Wilcoxon_KBWA < 0.05]),
  '-K' = unique(df_sign_age$Via[df_sign_age$Wilcoxon_K < 0.05]),

  BoD = unique(df_sign_age$Via[df_sign_age$Wilcoxon_DBo < 0.05]),
  RsD = unique(df_sign_age$Via[df_sign_age$Wilcoxon_DRs < 0.05]),
  bwaD = unique(df_sign_age$Via[df_sign_age$Wilcoxon_DBWA < 0.05]),
  '-D' = unique(df_sign_age$Via[df_sign_age$Wilcoxon_D < 0.05]),
  DdhD = unique(df_sign_age$Via[df_sign_age$Wilcoxon_DdhD< 0.05])
)

library(UpSetR)
library(ComplexHeatmap)
data_matrix <- fromList(data_list)

bar_colors <- rep("black", 47)
bar_colors[c(5,11,14)] <- "#56B4E9"
bar_colors[c(1,2)] <- c("#8470FF","#FFA500")

upset(
  data_matrix,
  sets = names(data_list),
  order.by = "freq",
  mainbar.y.label = "Nº Genera signif dif between age ranges",
  sets.x.label = "Nº Signif Genera Detected",
  text.scale = c(1.5, 1.3, 1.2, 1, 1.3, 1.3),
  point.size = 2,
  line.size = 0,
  keep.order = TRUE,  # Mantiene el orden original de los conjuntos
  sets.bar.color = "grey50",  # Todas las barras de conjuntos en negro
  matrix.color = "black",  # Matriz en negro
  main.bar.color = bar_colors,
  shade.color = "gray90",  # Fondo de sombreado más claro para contraste
  number.angles = 0,  # Mejor legibilidad de números
  mb.ratio = c(0.7, 0.3)
)


# >>>>>>>>>>>>>>>>>>> FIGURE 5C <<<<<<<<<<<<<<<<<<<<<<<
gen_sign_sex <- c("Propionibacterium", "Elizabethkingia", "Rhodopseudomonas", "Dermabacter", "Enterococcus",
                   "Aeromonas", "Flavobacterium", "Gemella")
df_gen_sign_sex <- df_sign_sex[df_sign_sex$Via %in% gen_sign_sex,]
df_gen_sign_sex <- df_gen_sign_sex[,-1]

df_long <- df_gen_sign_sex %>%
  pivot_longer(cols = -Via,
               names_to = c(".value", "Method"),
               names_sep = "_")


df_long$pValor <- ifelse(df_long$Wilcoxon<0.05, -log10(df_long$Wilcoxon), NA)

df_long$Method[df_long$Method == "KBo"] <- "BoK"
df_long$Method[df_long$Method == "KBWA"] <- "bwaK"
df_long$Method[df_long$Method == "KRs"] <- "RsK"
df_long$Method[df_long$Method == "K"] <- "-K"
df_long$Method[df_long$Method == "DBo"] <- "BoD"
df_long$Method[df_long$Method == "DBWA"] <- "bwaD"
df_long$Method[df_long$Method == "DRs"] <- "RsD"
df_long$Method[df_long$Method == "D"] <- "-D"


ggplot(df_long, aes(x = Method,
                    y = Via,
                    size = pValor,
                    fill = GrupoDominante
)) +
  geom_point(color = "black", shape = 21, alpha = 0.8, stroke = 0.5) +
  scale_size_continuous(
    name = "p-value",
    breaks = c(-log10(0.049), -log10(0.01), -log10(0.0015)),
    labels = c("0.05", "0.01", "0.001"),
    range = c(2, 12)
  ) +
  scale_y_discrete(limits = rev(c("Propionibacterium",
                                  "Dermabacter","Rhodopseudomonas",
                                  "Enterococcus", "Elizabethkingia",
                                  "Aeromonas", "Flavobacterium", "Gemella"))) +
  scale_x_discrete(
    labels = c("-D", "DdhD", "BoD", "RsD", "bwaD", " ", "-K", "BoK", "RsK", "bwaK"),  # “ ” es un espacio unicode más grande
    limits = c("-D", "DdhD", "BoD", "RsD", "bwaD", " ", "-K", "BoK", "RsK", "bwaK")
  ) +

  scale_fill_manual(name = "Dominant Group",
                    values = c("Femenino" = "#27408B", "Masculino" = "#CD4F39"),
                    labels = c("Femenino" = "Female", "Masculino" = "Male")) +
  theme_minimal() +
  labs(x = "Methodology", y = "Genus", title = "") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  ) + guides(
    fill = guide_legend(
      override.aes = list(size = 5)
    )
  )


