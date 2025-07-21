# >>>>>>>>>>>>>>>> FIGURE 6 A <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

NumberPathways <- read_excel(".../Fig6-NumberPathways.xlsx")

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


# >>>>>>>>>>>>>>>> FIGURE 6 B <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

load(".../Fig6b-DataList-Pathways.RData" )

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




