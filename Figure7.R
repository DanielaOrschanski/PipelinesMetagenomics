#>>>>>>>>>>>>>>>>> Figure 7a <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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

# Figure 7b <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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

