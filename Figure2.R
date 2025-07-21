#>>>>>>>>>>>>>>>>>>>>>>>>>> FIGURE 2A <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

library(ggplot2)
MappedReads_AllMethods <- read_excel(".../Fig2a-MappedReads.xlsx")
MappedReads_AllMethods$Method <- factor(MappedReads_AllMethods$Method, levels = c("Bowtie","BWA", "RSubread", "DRAGEN"))

ggplot(MappedReads_AllMethods, aes(x = Method, y = MappedReads, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_jitter(aes(color = Method), width = 0.2, size = 1.2, alpha = 0.7) +
  geom_boxplot(width = 0.6, outlier.shape = NA, color = "black", alpha = 0.8) +
  geom_line(aes(group=Sample), colour="black", linetype="11", alpha = 0.3) +
  theme_minimal() +
  labs(x = "", y = "% Mapped Reads", title = "") +
  scale_y_continuous(breaks = seq(0, 100, by =20))+
  scale_fill_manual(values = c(
    "Bowtie" = "#E69F00",
    "BWA" = "#56B4E9",
    "RSubread" = "#009E73",
    "DRAGEN" = "#9467BD"
  )) +
  scale_color_manual(values = c(
    "Bowtie" = "#E69F00",
    "BWA" = "#56B4E9",
    "RSubread" = "#009E73",
    "DRAGEN" = "#9467BD"
  )) +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))


#P-values:

library(dplyr)
library(tidyr)
library(purrr)

df_wide <- MappedReads_AllMethods %>%
  pivot_wider(names_from = Method, values_from = MappedReads)

df_wide$Dif_BWA_Rs <- ((df_wide$BWA - df_wide$RSubread ) / df_wide$BWA)*100

methods <- setdiff(colnames(df_wide), "Sample")
method_pairs <- combn(methods, 2, simplify = FALSE)

wilcox_results <- map(method_pairs, function(pair) {
  test <- wilcox.test(df_wide[[pair[1]]], df_wide[[pair[2]]], paired = FALSE)
  tibble(
    Method1 = pair[1],
    Method2 = pair[2],
    p_value = test$p.value
  )
}) %>% bind_rows()

wilcox_results <- wilcox_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

#>>>>>>>>>>>>>>>>>>>>>>>>>> FIGURE 2B <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GC_content_AllMethods <- read_excel(".../Fig2b-GC.xlsx")

ggplot(GC_content_AllMethods, aes(x = `GC Content`, y = Count, group = interaction(ID, Method), color = Method)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ Method, ncol = 4) +
  labs(
    title = "",
    x = "GC content (%)",
    y = "Number of sequences"
  ) +
  scale_color_manual(values = c(
    "Bowtie" = "#E69F00",
    "BWA" = "#56B4E9",
    "RSubread" = "#009E73",
    "No-dehost" = "#EE6363"
  )) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(breaks = seq(0, 200000, by = 25000)) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    strip.text = element_text( size = 15)
  )


GC_40 <- GC_content_AllMethods %>%
  filter(`GC Content` >= 30, `GC Content` <= 50)

GC_60 <- GC_content_AllMethods %>%
  filter(`GC Content` >= 50, `GC Content` <= 70)

max40 <- GC_40 %>%
  group_by(ID, Method) %>%
  summarise(MaxCount = max(Count), .groups = "drop") %>%
  mutate(GC_peak = "GC ~40%")

max60 <- GC_60 %>%
  group_by(ID, Method) %>%
  summarise(MaxCount = max(Count), .groups = "drop") %>%
  mutate(GC_peak = "GC ~60%")

maxs <- bind_rows(max40, max60)

maxs$Method <- factor(maxs$Method, levels = c("No-dehost", "Bowtie", "BWA", "RSubread"))
maxs$GC_peak <- factor(maxs$GC_peak, levels = c("GC ~40%", "GC ~60%"))

maxs <- as.data.frame(maxs)
maxs$ID <- as.factor(maxs$ID)

wilcox_results_peak <- maxs %>%
  filter(Method %in% c("BWA", "RSubread", "Bowtie")) %>%
  group_by(GC_peak) %>%
  nest() %>%
  mutate(
    data_wide = map(data, ~ pivot_wider(.x, names_from = Method, values_from = MaxCount)),

    # Wilcoxon BWA vs RSubread
    test_BWA_Rs = map(data_wide, ~ wilcox.test(.x$BWA, .x$RSubread, paired = TRUE)),
    p_BWA_Rs = map_dbl(test_BWA_Rs, ~ .x$p.value),

    # Wilcoxon BWA vs Bowtie
    test_BWA_Bo = map(data_wide, ~ wilcox.test(.x$BWA, .x$Bowtie, paired = TRUE)),
    p_BWA_Bo = map_dbl(test_BWA_Bo, ~ .x$p.value),

    # Wilcoxon RSubread vs Bowtie
    test_Rs_Bo = map(data_wide, ~ wilcox.test(.x$RSubread, .x$Bowtie, paired = TRUE)),
    p_Rs_Bo = map_dbl(test_Rs_Bo, ~ .x$p.value),

  )

wide <- wilcox_results_peak$data_wide[[1]]
wide_60 <- wilcox_results_peak$data_wide[[which(wilcox_results_peak$GC_peak == "GC ~60%")]]
wide_40 <- wilcox_results_peak$data_wide[[which(wilcox_results_peak$GC_peak == "GC ~40%")]]

wilcox.test(wide_60$BWA, wide_60$RSubread, paired = TRUE)
wilcox.test(wide_60$Bowtie, wide_60$RSubread, paired = TRUE)

wide_60$Bo_BWA <- (wide_60$Bowtie - wide_60$BWA) / wide_60$Bowtie *100
wide_40$Bo_BWA <- (wide_40$Bowtie - wide_40$BWA) / wide_40$Bowtie *100

#>>>>>>>>>>>>>>>>>>>>>>>>>> FIGURE 2C <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

HomoSapiens_AllMethods <-read_excel(".../Fig2c-HomoSapiens.xlsx")
HomoSapiens_AllMethods$Method <- factor(HomoSapiens_AllMethods$Method, levels = c("-K","-D", "DdhD",  "BoK", "BoD",
                                                    "bwaK", "bwaD", "RsK", "RsD"))

ggplot(HomoSapiens_AllMethods, aes(x = Method, y = Abundance, fill = Method)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = NA) +  # Violin con transparencia
  geom_jitter(aes(color = Method), width = 0.2, size = 1.2, alpha = 0.4) +  # Jitter sutil y transparente
  geom_boxplot(width = 0.6, outlier.shape = NA, color = "black", alpha = 0.8) +  # Boxplot clásico
  theme_minimal() +
  labs(x = "", y = "Relative Abundance (log10)", title = "") +
  scale_y_log10() +
  scale_fill_manual(values = c(
    "BoK" = "#E69F00",  # Naranja suave
    "BoD" = "#F0C987",

    "bwaK" = "#56B4E9", # Azul suave
    "bwaD" = "#A0D8F0",

    "RsK" = "#009E73",  # Verde suave
    "RsD" = "#8FCB88",

    "-K" = "#666666",  # Gris oscuro
    "DdhD" = "#9467BD", # Violeta suave
    "-D" = "#999999"
  )) +
  scale_color_manual(values = c(
    "BoK" = "#E69F00", "BoD" = "#F0C987",
    "bwaK" = "#56B4E9", "bwaD" = "#A0D8F0",
    "RsK" = "#009E73", "RsD" = "#8FCB88",
    "-K" = "#666666", "DdhD" = "#9467BD", "-D" = "#999999"
  )) +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))  # Rotar etiquetas suavemente


#P-values:

library(dplyr)
library(tidyr)
library(purrr)

method_pairs <- combn(unique(HomoSapiens_AllMethods$Method), 2, simplify = FALSE)

get_pvalue <- function(df, method1, method2) {
  df_pair <- df %>%
    filter(Method %in% c(method1, method2)) %>%
    pivot_wider(names_from = Method, values_from = Abundance)

  test <- wilcox.test(df_pair[[method1]], df_pair[[method2]], paired = TRUE)

  diff_pct <- mean(((df_pair[[method1]] - df_pair[[method2]]) / df_pair[[method1]]) * 100)

  tibble(
    Method1 = method1,
    Method2 = method2,
    p_value = test$p.value,
    diff_pct = diff_pct
  )
}

wilcoxon_results <- map_dfr(method_pairs, ~ get_pvalue(HomoSapiens_AllMethods, .x[1], .x[2])) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

