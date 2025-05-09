library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(purrr)

df_wide <- read_excel("Daniela/Biota/PipelineBiota/paraPaper/Scripts Reproducir Paper/Fig4-RA6genera.xlsx")

# Pares of methods to compare
#Fig4a:
comparisons <- list(
  c("BoK", "BoD"),
  c("bwaK", "bwaD"),
  c("RsK", "RsD"),
  c("-K", "-D")
)

#Fig4b:
comparisons <- list(
  c("BoK", "-K"),
  c("BoK", "RsK"),
  c("BoK", "bwaK")
)

#Fig4c:
comparisons <- list(
  c("DdhD", "-D"),
  c("DdhD", "BoD"),
  c("DdhD", "RsD"),
  c("DdhD", "bwaD")
)


df_long <- df_wide %>%
  pivot_longer(cols = -Genus,
               names_to = c("Sample", "Methodology"),
               names_sep = "_",
               values_to = "Value")

unique(df_long$Methodology)

# Generate data for Bland-Altman plots:
df_bland <- bind_rows(lapply(comparisons, function(par) {
  df_long %>%
    filter(Methodology %in% par) %>%
    pivot_wider(names_from = Methodology, values_from = Value) %>%
    mutate(
      Comparison = paste(par[1], "vs", par[2]),
      Mean = (as.numeric(.data[[par[1]]]) + as.numeric(.data[[par[2]]])) / 2,
      Difference = as.numeric(.data[[par[1]]]) - as.numeric(.data[[par[2]]]),
      DifferencePerc = (as.numeric(.data[[par[1]]]) - as.numeric(.data[[par[2]]])) / as.numeric(.data[[par[1]]]) * 100
    )
}))

# Calculate limits of Bland-Altman for each comparison:
df_lims <- df_bland %>%
  group_by(Genus, Comparison) %>%
  summarise(
    Diff_mean = mean(Difference, na.rm = TRUE),
    Diff_SD = sd(Difference, na.rm = TRUE),
    Upper_limit = Diff_mean + 1.96 * Diff_SD,
    Lower_limit = Diff_mean - 1.96 * Diff_SD
  )

df_final <- df_bland %>%
  left_join(df_lims, by =c("Genus", "Comparison"))

df_final <- df_final %>%
  mutate(Comparison = as.factor(Comparison))


library(ggplot2)
library(gridExtra)
library(dplyr)
library(purrr)

#Plots:
plots <- df_final %>%
  split(.$Genus) %>%
  map(~ ggplot(.x, aes(x = Mean, y = Difference, color = Comparison, shape = Comparison)) +
        geom_point(size = 1.5, alpha = 0.8) +
        geom_hline(aes(yintercept = Diff_mean, color = Comparison), linetype = "dashed", linewidth = 0.4) +
        geom_hline(aes(yintercept = Upper_limit, color = Comparison), linetype = "dotted", linewidth = 0.4) +
        geom_hline(aes(yintercept = Lower_limit, color = Comparison), linetype = "dotted", linewidth = 0.4) +
        labs(title = unique(.x$Genus)) +
        theme_classic() +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      +
        scale_color_manual(values = c("black", "blue", "red", "green")) +
        scale_shape_manual(values = c(16, 17, 15, 18))
  )

# Extracte legends:
legend_plot <- ggplot(df_final, aes(x = Mean, y = Difference, color = Comparison, shape = Comparison)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("black", "blue", "red", "green")) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  theme_minimal() +
  theme(legend.position = "right")

legend <- cowplot::get_legend(legend_plot)

#Adjust numbers to align graphs:
plots <- lapply(plots, function(p) {
  p + scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
})

Fig4a <- plots
legend_4a <- legend

Fig4b <- plots
legend_4b <- legend

Fig4c <- plots
legend_4c <- legend

grid.arrange(
  grobs = Fig4a,
  ncol = 2,
  bottom = textGrob("Average Relative Abundance", gp = gpar(fontsize = 12)),
  #Figure 4a:
  left = textGrob("Difference: Kraken - DRAGEN", gp = gpar(fontsize = 12), rot = 90),
  right = legend
)

grid.arrange(
  grobs = Fig4b,
  ncol = 2,
  bottom = textGrob("Average Relative Abundance", gp = gpar(fontsize = 12)),
  left = textGrob("Difference: BoK vs Other Kraken Methods", gp = gpar(fontsize = 12), rot = 90),
  right = legend
)


grid.arrange(
  grobs = Fig4c,
  ncol = 2,
  bottom = textGrob("Average Relative Abundance", gp = gpar(fontsize = 12)),
  #Figure 4c:
  left = textGrob("Difference: DdhD vs Other DRAGEN Methods", gp = gpar(fontsize = 12), rot = 90),
  right = legend
)

grid.arrange(
  grobs = c(Fig4a, Fig4b, Fig4c),
  ncol = 3,
  bottom = textGrob("Average Relative Abundance", gp = gpar(fontsize = 12)),

  left = c(
    textGrob("Difference: Kraken - DRAGEN", gp = gpar(fontsize = 12), rot = 90),
    textGrob("Difference: BoK vs Other Kraken Methods", gp = gpar(fontsize = 12), rot = 90),
    textGrob("Difference: DdhD vs Other DRAGEN Methods", gp = gpar(fontsize = 12), rot = 90)
           ),
  right = legend
)

library(gridExtra)
library(grid)

# Primero apilamos los gráficos verticalmente dentro de cada grupo
stacked_Fig4a <- arrangeGrob(grobs = Fig4a, ncol = 1)
stacked_Fig4b <- arrangeGrob(grobs = Fig4b, ncol = 1)
stacked_Fig4c <- arrangeGrob(grobs = Fig4c, ncol = 1)

# Luego unimos el eje y con cada stack
col1 <- arrangeGrob(
  textGrob("Difference: Kraken - DRAGEN", gp = gpar(fontsize = 12), rot = 90),
  stacked_Fig4a,
 # legend_4a,
  ncol = 2,
  widths = c(1, 5)
)

col2 <- arrangeGrob(
  textGrob("Difference: BoK vs Other Kraken Methods", gp = gpar(fontsize = 12), rot = 90),
  stacked_Fig4b,
  #legend_4b,
  ncol = 2,
  widths = c(1, 5)
)

col3 <- arrangeGrob(
  textGrob("Difference: DdhD vs Other DRAGEN Methods", gp = gpar(fontsize = 12), rot = 90),
  stacked_Fig4c,
  #legend_4c,
  ncol = 2,
  widths = c(1, 5)
)

# Finalmente los unimos todos
final_plot <- arrangeGrob(
  col1, col2, col3,
  ncol = 3,
  bottom = textGrob("Average Relative Abundance", gp = gpar(fontsize = 12))
  #right = c(legend_4a, legend_4b, legend_4c)
)

# Mostrar
grid.newpage()
grid.draw(final_plot)




#Statistics for Figure 4a: ---------------------------
#Overestimated by kraken:
df_final %>%
  filter(Genus == "Cutibacterium", DifferencePerc > 0) %>%
  group_by(Comparison) %>%
  summarise(n_samples_over = n_distinct(Sample))
df_final %>%
  filter(Genus == "Cutibacterium") %>%
  group_by(Comparison) %>%
  summarise(
    n_samples = n_distinct(Sample),
    avg_diff_perc = mean(DifferencePerc, na.rm = TRUE),
    .groups = "drop"
  )

#Underestimated by kraken:
df_final %>%
  filter(Genus == "Acinetobacter", DifferencePerc < 0) %>%
  group_by(Comparison) %>%
  summarise(n_samples_over = n_distinct(Sample))
df_final %>%
  filter(Genus == "Acinetobacter") %>%
  group_by(Comparison) %>%
  summarise(
    n_samples = n_distinct(Sample),
    avg_diff_perc = mean(DifferencePerc, na.rm = TRUE),
    .groups = "drop"
  )

#-----------------------------------------------------

# p-values with Wilcoxon
calculate_p <- function(df, method1, method2) {
  df %>%
    group_by(Genus) %>%
    summarise(
      p_value = ifelse(
        sum(!is.na(get(method1))) > 1 & sum(!is.na(get(method2))) > 1, # Asegura que haya más de un valor
        wilcox.test(get(method1), get(method2), paired = TRUE)$p.value,
        NA
      )
    ) %>%
    mutate(Comparacion = paste(method1, "vs", method2))
}

res <- lapply(comparisons, function(par) calculate_p(df_final, par[1], par[2]))
df_pvalues <- bind_rows(res)
