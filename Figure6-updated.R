# >>>>>>>>>>>>>>>> FIGURE 6 A <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

NumberPathways <- read_excel("/media/4tb2/Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig6-NumberPathways.xlsx")

colors_de_host <- c("BoH"= "#E69F00",
                    "bwaH" = "#A52A2A",
                    "RsH" = "#009E73",
                    "-H" = "#9467BD" )

ggplot(NumberPathways, aes(x = de_host, y = Pathway, fill = de_host)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +  # Violin plot con transparencia
  geom_jitter(aes(color = de_host), width = 0.2, size = 1.2, alpha = 0.7) +  # Jitter con mismo color que boxplot
  geom_boxplot(width = 0.6, outlier.shape = NA, color = "black", alpha = 0.8) +  # Boxplot sin outliers visibles
  geom_line(aes(group = ID), colour = "black", linetype = "11", alpha = 0.4) +  # Conectar mismo ID

  theme_minimal(base_size = 15) +
  labs(x = "",
       y = "Number of Identified Pathways",
       title = "") +  # Título en negrita
  theme(
    plot.title = element_text(size = 14, face = "bold"),  # Título en negrita
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",       # Mueve la leyenda abajo
    legend.direction = "vertical"   # Ordena la leyenda en una fila
  ) +
  scale_fill_manual(values = colors_de_host, name = "De-Host") +
  scale_color_manual(values = colors_de_host, guide = "none")  # Usa los mismos colores para jitter y oculta la leyenda extra

groups <- unique(NumberPathways$de_host)
pares <- combn(groups, 2, simplify = FALSE)

results <- lapply(pares, function(par) {
  data_par <- NumberPathways[NumberPathways$de_host %in% par, ]

  # Same order of ID
  data_par1 <- data_par[data_par$de_host == par[1], ]
  data_par2 <- data_par[data_par$de_host == par[2], ]

  data_par1 <- data_par1[order(data_par1$ID), ]
  data_par2 <- data_par2[order(data_par2$ID), ]

  all(data_par1$ID == data_par2$ID)
  if(all((data_par1$Pathway - data_par2$Pathway) == 0)) {
    p_v <- "Non-significat"
  } else {
    # Paired wilcoxon test
    test <- wilcox.test(data_par1$Pathway, data_par2$Pathway, paired = TRUE)
    p_v <- test$p.value
  }

  data.frame(Group1 = par[1],
             Group2 = par[2],
             p.value = p_v )
})


results_df <- do.call(rbind, results)
results_df$p_adj <- p.adjust(results_df$p.value, method = "BH")
print(results_df)

# Reviewed Fig 6a : everyone against BoH ---------------------
library(dplyr)

# Crear una tabla con solo los valores de BoH
BoH_values <- NumberPathways %>%
  filter(de_host == "BoH") %>%
  select(ID, Pathway_BoH = Pathway)

# Unir con los demás métodos
NumberPathways_diff <- NumberPathways %>%
  filter(de_host != "BoH") %>%
  left_join(BoH_values, by = "ID") %>%
  mutate(Diff_to_BoH = Pathway - Pathway_BoH)

ggplot(NumberPathways_diff, aes(x = de_host, y = Diff_to_BoH, fill = de_host)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_jitter(aes(color = de_host), width = 0.2, size = 1.2, alpha = 0.7) +
  geom_boxplot(width = 0.6, outlier.shape = NA, color = "black", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +  # Línea base en 0
  
  theme_minimal(base_size = 15) +
  labs(x = "",
       y = "Difference in Number of Pathways (vs BoH)",
       title = "") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    legend.direction = "vertical"
  ) +
  scale_fill_manual(values = colors_de_host, name = "De-Host") +
  scale_color_manual(values = colors_de_host, guide = "none")



# >>>>>>>>>>>>>>>> FIGURE 6 B <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

load("~/Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig6b-DataList-Pathways.RData" )

ggvenn(data_list,
       fill_color = c(  "#9467BD", "#E69F00","#A52A2A", "#009E73" ),
       stroke_size = 0,
       set_name_size = 4,
       text_size = 3)

vias_data <- tibble(
  Vias = unique(unlist(data_list)),
  Bo = as.integer(Vias %in% data_list$BoH),
  BWA = as.integer(Vias %in% data_list$bwaH),
  Rs = as.integer(Vias %in% data_list$RsH),
  SinDH = as.integer(Vias %in% data_list$`-H`)
)

vias_summary <- vias_data %>%
  group_by(Bo, BWA, Rs, SinDH) %>%
  summarise(Vias_List = list(Vias), Count = n(), .groups = "drop")


#>>>>>>>>>>>>>>>>> Figure 6c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
df_long <-  read_excel(".../Fig7a-Pathways-Sex.xlsx")

ggplot(df_long, aes(x = Method,
                    y = Via,
                    size = SizeP,  # Controla el tamaño con -log10(Wilcoxon)
                    fill = GrupoDominante,
                    shape = GrupoDominante)) +  # Forma según GrupoDominante
  geom_point(color = "black", alpha = 0.8) +
  scale_size_continuous(
    name = "p-value",
    breaks = c(-log10(0.043), -log10(0.01), -log10(0.0035)),  # Puntos clave de significancia
    labels = c("0.05", "0.01", "0.001"),  # Mostrar valores originales de Wilcoxon
    range = c(2, 12)  # Controla el tamaño de los puntos
  ) +
  scale_y_discrete(limits = rev(c("PWY-6922", "ARGININE-SYN4-PWY","SER-GLYSYN-PWY",
                                  "RIBOSYN2-PWY", "P4-PWY",
                                  "THRESYN-PWY","PWY-7761" )
  )) +
  scale_x_discrete(limits = (c("-H", "BoH", "RsH" ,"bwaH"))) +
  scale_fill_manual(values = c("Femenino" = "grey", "Masculino" = "#CD4F39"),
                    guide = "none",
                    name = "Dominant Group",
                    labels = c("Femenino" = "Female", "Masculino" = "Male")) +
  scale_shape_manual(name = "Sex",
                     values = c("Femenino" = 21, "Masculino" = 21),
                     labels = c("Femenino" = "Female", "Masculino" = "Male"),
                     guide = "none") +
  theme_minimal() +
  labs(x = "Methodology", y = "Pathways",
       fill = "Dominant Group",
       title = "") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size=10)
    #,
    #legend.position = "bottom",
    #legend.box = "vertical"
  )+ guides(
    size = guide_legend(
      override.aes = list(
        shape = 21,         # círculo
        fill = "white",     # blanco adentro
        color = "black"     # borde negro
      )
    ),
    fill = guide_legend( override.aes = list(
      shape = 21,
      size = 5
    )
    )
  )

# Figure 6d <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
df_long <- read_excel(".../Fig7b-Pathways-Age.xlsx")

#TABLE S3 -------------------

df_combined <- bind_rows(lapply(names(lista_vias), function(nombre) {
  lista_vias[[nombre]] %>%
    mutate(Metodologia = nombre)
}), .id = "ID")
df_combined <- df_combined[,-1]

#sex:
vias_select <- df_sign_sexo$Via

#age:
vias_select <- df_sign_rango$Via

df_filtered <- df_combined[df_combined$Pathway %in% vias_select,]
colnames(df_filtered)[1] <- "Pathway"

df_long <- df_filtered %>%
  pivot_longer(cols = -c(Pathway, Clase, Description, Metodologia),
               names_to = "ID",
               values_to = "Abundancia")
library(tidyverse)

MetadataB <- as.data.frame(read_excel("~/Daniela/Biota/PipelineBiota/data/Metadata-soloColumnasUsables.xlsx"))
df_merged <- df_long %>%
  left_join(MetadataB[, c("ID", "Sexo", "Rango etario")], by = "ID")  # Asegúrate de que MetadataB tenga una columna "Muestra"

str(df_merged)
df_merged$Sexo <- as.factor(df_merged$Sexo)
df_merged$`Rango etario` <- factor( df_merged$`Rango etario` , levels = c("18-35", "35-55", ">55"))
df_merged$Metodologia <- as.factor(df_merged$Metodologia)

p_values_df <- df_merged %>%
  group_by(Pathway, Metodologia) %>%
  summarise(p_value = wilcox.test(Abundancia ~ Sexo)$p.value) %>%
  mutate(p_label = sprintf("p = %.3f", p_value))  # Formatear el p-valor
p_values_df$simbolo <- ifelse(p_values_df$p_value < 0.05, ifelse(p_values_df$p_value < 0.01,"**", "*"), "-")


p_values_wide <- p_values_df[, -which(colnames(p_values_df) %in% c("p_label", "simbolo", "max_value", "y_position"))] %>%
  pivot_wider(names_from = Metodologia,
              values_from = p_value,
              names_sort = TRUE)
CPM_Vias_Clases_Anotacion_DifSign_82p <- read_excel("~/Daniela/Biota/CPM_Vias_Clases_Anotacion_DifSign_82p.xlsx")
colnames(CPM_Vias_Clases_Anotacion_DifSign_82p)[1] <- "Pathway"
p_values_wide <- merge(p_values_wide, CPM_Vias_Clases_Anotacion_DifSign_82p[, c("Pathway", "Clase", "Description")], by = "Pathway")
p_values_wide[, c(2:5)] <- round(p_values_wide[,c(2:5)], 3)


