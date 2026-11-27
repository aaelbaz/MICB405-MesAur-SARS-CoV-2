library(tidyverse)

# Automatically switch working directory to this script's folder if needed
if (!dir.exists("data")) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_path <- if (length(file_arg) > 0) gsub("~\\+~", " ", sub("^--file=", "", file_arg[1])) else NULL
  if (is.null(script_path)) {
    for (f in rev(sys.frames())) {
      if (!is.null(f$ofile)) { script_path <- f$ofile; break }
    }
  }
  if (is.null(script_path) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    script_path <- rstudioapi::getActiveDocumentContext()$path
  }
  if (!is.null(script_path) && nzchar(script_path)) setwd(dirname(script_path))
}

dir.create("figures", showWarnings = FALSE)

# 1. 4 dpi: Mock vs. Cov 
dat_4dpi <- read_csv("data/filtered and sorted/4dpi_top5_up_down.csv") %>%
  mutate(
    Direction = if_else(log2FoldChange > 0, "UP", "DOWN"),
    sig_star  = case_when(padj < 1e-10 ~ "***", padj < 1e-4 ~ "**", TRUE ~ "*"),
    bar_label = sprintf("%+.2f %s", log2FoldChange, sig_star)
  ) %>%
  arrange(log2FoldChange)

order_4dpi <- dat_4dpi %>% dplyr::select(Gene) %>% pull()
xmax_4dpi  <- max(abs(dat_4dpi$log2FoldChange) + dat_4dpi$lfcSE) * 1.38

p_4dpi <- dat_4dpi %>%
  ggplot(aes(x = Gene, y = log2FoldChange, fill = Direction)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "grey50") +
  geom_col(color = "black", width = 0.72) +
  geom_errorbar(
    aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE),
    width = 0.28, linewidth = 0.5
  ) +
  geom_text(
    aes(
      y     = if_else(log2FoldChange > 0, log2FoldChange + lfcSE, log2FoldChange - lfcSE),
      label = bar_label,
      hjust = if_else(log2FoldChange > 0, -0.15, 1.15)
    ),
    size = 3.2, fontface = "bold"
  ) +
  scale_x_discrete(limits = order_4dpi) +
  coord_flip(ylim = c(-xmax_4dpi, xmax_4dpi)) +
  scale_fill_manual(values = c("UP" = "#D95F02", "DOWN" = "#1F78B4")) +
  labs(title = "4 dpi: Mock vs. Cov",
       x = "Gene", y = "log2(Fold Change)", fill = "Regulation") +
  theme_light()

ggsave("figures/4dpi_top5_barchart.png", plot = p_4dpi, width = 7, height = 4.5, dpi = 300)


# 2. 80 dpi: Mock vs. Cov 
dat_80dpi <- read_csv("data/filtered and sorted/80dpi_top5_up_down.csv") %>%
  mutate(
    Direction = if_else(log2FoldChange > 0, "UP", "DOWN"),
    sig_star  = case_when(padj < 1e-10 ~ "***", padj < 1e-4 ~ "**", TRUE ~ "*"),
    bar_label = sprintf("%+.2f %s", log2FoldChange, sig_star)
  ) %>%
  arrange(log2FoldChange)

order_80dpi <- dat_80dpi %>% dplyr::select(Gene) %>% pull()
xmax_80dpi  <- max(abs(dat_80dpi$log2FoldChange) + dat_80dpi$lfcSE) * 1.38

p_80dpi <- dat_80dpi %>%
  ggplot(aes(x = Gene, y = log2FoldChange, fill = Direction)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "grey50") +
  geom_col(color = "black", width = 0.72) +
  geom_errorbar(
    aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE),
    width = 0.28, linewidth = 0.5
  ) +
  geom_text(
    aes(
      y     = if_else(log2FoldChange > 0, log2FoldChange + lfcSE, log2FoldChange - lfcSE),
      label = bar_label,
      hjust = if_else(log2FoldChange > 0, -0.15, 1.15)
    ),
    size = 3.2, fontface = "bold"
  ) +
  scale_x_discrete(limits = order_80dpi) +
  coord_flip(ylim = c(-xmax_80dpi, xmax_80dpi)) +
  scale_fill_manual(values = c("UP" = "#D95F02", "DOWN" = "#1F78B4")) +
  labs(title = "80 dpi: Mock vs. Cov",
       x = "Gene", y = "log2(Fold Change)", fill = "Regulation") +
  theme_light()

ggsave("figures/80dpi_top5_barchart.png", plot = p_80dpi, width = 7, height = 4.5, dpi = 300)


# 3. 80dpi vs 4dpi cov (Top 5 UP/DOWN)
dat_80vs4 <- read_csv("data/filtered and sorted/80dpi_vs_4dpi_cov_top5_up_down.csv") %>%
  mutate(
    Direction = if_else(log2FoldChange > 0, "UP", "DOWN"),
    sig_star  = case_when(padj < 1e-10 ~ "***", padj < 1e-4 ~ "**", TRUE ~ "*"),
    bar_label = sprintf("%+.2f %s", log2FoldChange, sig_star)
  ) %>%
  arrange(log2FoldChange)

order_80vs4 <- dat_80vs4 %>% dplyr::select(Gene) %>% pull()
xmax_80vs4  <- max(abs(dat_80vs4$log2FoldChange) + dat_80vs4$lfcSE) * 1.38

p_80vs4 <- dat_80vs4 %>%
  ggplot(aes(x = Gene, y = log2FoldChange, fill = Direction)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "grey50") +
  geom_col(color = "black", width = 0.72) +
  geom_errorbar(
    aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE),
    width = 0.28, linewidth = 0.5
  ) +
  geom_text(
    aes(
      y     = if_else(log2FoldChange > 0, log2FoldChange + lfcSE, log2FoldChange - lfcSE),
      label = bar_label,
      hjust = if_else(log2FoldChange > 0, -0.15, 1.15)
    ),
    size = 3.2, fontface = "bold"
  ) +
  scale_x_discrete(limits = order_80vs4) +
  coord_flip(ylim = c(-xmax_80vs4, xmax_80vs4)) +
  scale_fill_manual(values = c("UP" = "#D95F02", "DOWN" = "#1F78B4")) +
  labs(title = "80dpi vs 4dpi cov",
       x = "Gene", y = "log2(Fold Change)", fill = "Regulation") +
  theme_light()

ggsave("figures/80dpi_vs_4dpi_cov_top5_barchart.png", plot = p_80vs4, width = 7, height = 4.5, dpi = 300)
