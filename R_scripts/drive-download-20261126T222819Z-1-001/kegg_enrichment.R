library(tidyverse)
library(clusterProfiler)

# 1. Load genes from KEGG for Mesocricetus auratus
maua_raw <- read_tsv(
  "https://rest.kegg.jp/list/maua",
  col_names = c("kegg_id", "type", "location", "desc"),
  show_col_types = FALSE
)

# 2. Get KEGG ID 
kegg_map <- maua_raw %>%
  mutate(
    entrez_id = str_remove(kegg_id, "^maua:"),
    sym_part  = if_else(
      str_detect(desc, ";"),
      str_extract(desc, "^[^;]+"),
      paste0("LOC", entrez_id)
    )
  ) %>%
  separate_longer_delim(sym_part, delim = ", ") %>%
  mutate(Gene = str_trim(sym_part)) %>%
  dplyr::select(Gene, entrez_id) %>%
  distinct(Gene, .keep_all = TRUE)


# ==============================================================================
# 3. 4 dpi: Mock vs. Cov
# ==============================================================================
dat_4dpi <- read_csv("data/4dpi_results_filtered.csv", show_col_types = FALSE) %>%
  rename(Gene = 1) %>%
  left_join(kegg_map, by = "Gene") %>%
  mutate(
    entrez_id = if_else(
      is.na(entrez_id) & str_detect(Gene, "^LOC[0-9]+$"),
      str_remove(Gene, "^LOC"),
      entrez_id
    )
  ) %>%
  filter(!is.na(entrez_id))

up_4dpi <- dat_4dpi %>%
  filter(log2FoldChange >= 1, padj <= 0.05) %>%
  pull(entrez_id)

dn_4dpi <- dat_4dpi %>%
  filter(log2FoldChange <= -1, padj <= 0.05) %>%
  pull(entrez_id)

kegg_4dpi_up <- enrichKEGG(
  gene         = up_4dpi,
  organism     = "maua",
  keyType      = "ncbi-geneid",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

kegg_4dpi_dn <- enrichKEGG(
  gene         = dn_4dpi,
  organism     = "maua",
  keyType      = "ncbi-geneid",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

write_csv(as.data.frame(kegg_4dpi_up), "data/kegg_4dpi_up.csv")
write_csv(as.data.frame(kegg_4dpi_dn), "data/kegg_4dpi_dn.csv")

# 4 dpi UP plot
kegg_4up_df <- as.data.frame(kegg_4dpi_up) %>%
  filter(category != "Human Diseases") %>%
  head(15) %>%
  mutate(
    Direction     = "UP",
    GeneRatio_num = Count / DESize,
    Label         = str_wrap(Description, 45)
  )

p_kegg_4up <- kegg_4up_df %>%
  mutate(Label = fct_reorder(Label, GeneRatio_num)) %>%
  ggplot(aes(x = GeneRatio_num, y = Label, color = p.adjust)) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(3, 8)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    title = "KEGG Enrichment: 4 dpi Upregulated (Cov vs. Mock)",
    x     = "Gene Ratio",
    y     = NULL,
    color = "p.adjust",
    size  = "Gene Count"
  ) +
  theme_light() +
  theme(
    axis.text.y  = element_text(color = "black", size = 9.5),
    panel.border = element_rect(color = "black")
  )

ggsave("figures/kegg_4dpi_up.png", plot = p_kegg_4up, width = 9, height = 6.3, dpi = 300)

# 4 dpi DOWN plot
kegg_4dn_df <- as.data.frame(kegg_4dpi_dn) %>%
  filter(category != "Human Diseases") %>%
  head(15) %>%
  mutate(
    Direction     = "DOWN",
    GeneRatio_num = Count / DESize,
    Label         = str_wrap(Description, 45)
  )

p_kegg_4dn <- kegg_4dn_df %>%
  mutate(Label = fct_reorder(Label, GeneRatio_num)) %>%
  ggplot(aes(x = GeneRatio_num, y = Label, color = p.adjust)) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(3, 8), breaks = c(6, 8, 10, 12)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    title = "KEGG Enrichment: 4 dpi Downregulated (Cov vs. Mock)",
    x     = "Gene Ratio",
    y     = NULL,
    color = "p.adjust",
    size  = "Gene Count"
  ) +
  theme_light() +
  theme(
    axis.text.y  = element_text(color = "black", size = 9.5),
    panel.border = element_rect(color = "black")
  )

ggsave("figures/kegg_4dpi_dn.png", plot = p_kegg_4dn, width = 9, height = 3.4, dpi = 300)

# 4 dpi combined UP & DOWN plot
kegg_4comb_df <- bind_rows(kegg_4dn_df, kegg_4up_df) %>%
  mutate(
    Direction = factor(Direction, levels = c("DOWN", "UP")),
    Label     = fct_reorder(Label, GeneRatio_num, .fun = max)
  )

p_kegg_4comb <- kegg_4comb_df %>%
  ggplot(aes(x = GeneRatio_num, y = Label, color = p.adjust)) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(3, 8)) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_grid(. ~ Direction) +
  labs(
    title = "KEGG Enrichment: 4 dpi Upregulated and Downregulated (Cov vs. Mock)",
    x     = "Gene Ratio",
    y     = NULL,
    color = "p.adjust",
    size  = "Gene Count"
  ) +
  theme_light() +
  theme(
    axis.text.y      = element_text(color = "black", size = 9),
    panel.border     = element_rect(color = "black"),
    strip.background = element_rect(color = "black"),
    panel.spacing    = unit(0.8, "lines")
  )

ggsave("figures/kegg_4dpi_up_down.png", plot = p_kegg_4comb, width = 10.5, height = 6.6, dpi = 300)


# ==============================================================================
# 4. 80 dpi: Mock vs. Cov
# ==============================================================================
dat_80dpi <- read_csv("data/80dpi_results_filtered.csv", show_col_types = FALSE) %>%
  rename(Gene = 1) %>%
  left_join(kegg_map, by = "Gene") %>%
  mutate(
    entrez_id = if_else(
      is.na(entrez_id) & str_detect(Gene, "^LOC[0-9]+$"),
      str_remove(Gene, "^LOC"),
      entrez_id
    )
  ) %>%
  filter(!is.na(entrez_id))

up_80dpi <- dat_80dpi %>%
  filter(log2FoldChange >= 1, padj <= 0.05) %>%
  pull(entrez_id)

dn_80dpi <- dat_80dpi %>%
  filter(log2FoldChange <= -1, padj <= 0.05) %>%
  pull(entrez_id)

kegg_80dpi_up <- enrichKEGG(
  gene         = up_80dpi,
  organism     = "maua",
  keyType      = "ncbi-geneid",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

kegg_80dpi_dn <- enrichKEGG(
  gene         = dn_80dpi,
  organism     = "maua",
  keyType      = "ncbi-geneid",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

write_csv(as.data.frame(kegg_80dpi_up), "data/kegg_80dpi_up.csv")
write_csv(as.data.frame(kegg_80dpi_dn), "data/kegg_80dpi_dn.csv")

# 80 dpi UP plot
kegg_80up_df <- as.data.frame(kegg_80dpi_up) %>%
  filter(category != "Human Diseases") %>%
  head(15) %>%
  mutate(
    Direction     = "UP",
    GeneRatio_num = Count / DESize,
    Label         = str_wrap(Description, 45)
  )

p_kegg_80up <- kegg_80up_df %>%
  mutate(Label = fct_reorder(Label, GeneRatio_num)) %>%
  ggplot(aes(x = GeneRatio_num, y = Label, color = p.adjust)) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(3, 8), breaks = c(3, 4, 5)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    title = "KEGG Enrichment: 80 dpi Upregulated (Cov vs. Mock)",
    x     = "Gene Ratio",
    y     = NULL,
    color = "p.adjust",
    size  = "Gene Count"
  ) +
  theme_light() +
  theme(
    axis.text.y  = element_text(color = "black", size = 9.5),
    panel.border = element_rect(color = "black")
  )

ggsave("figures/kegg_80dpi_up.png", plot = p_kegg_80up, width = 9, height = 3.8, dpi = 300)

# 80 dpi DOWN plot
kegg_80dn_df <- as.data.frame(kegg_80dpi_dn) %>%
  filter(category != "Human Diseases") %>%
  head(15) %>%
  mutate(
    Direction     = "DOWN",
    GeneRatio_num = Count / DESize,
    Label         = str_wrap(Description, 45)
  )

p_kegg_80dn <- kegg_80dn_df %>%
  mutate(Label = fct_reorder(Label, GeneRatio_num)) %>%
  ggplot(aes(x = GeneRatio_num, y = Label, color = p.adjust)) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(3, 8), breaks = c(4, 5, 6, 7, 8)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    title = "KEGG Enrichment: 80 dpi Downregulated (Cov vs. Mock)",
    x     = "Gene Ratio",
    y     = NULL,
    color = "p.adjust",
    size  = "Gene Count"
  ) +
  theme_light() +
  theme(
    axis.text.y  = element_text(color = "black", size = 9.5),
    panel.border = element_rect(color = "black")
  )

ggsave("figures/kegg_80dpi_dn.png", plot = p_kegg_80dn, width = 9, height = 4.4, dpi = 300)

# 80 dpi combined UP & DOWN plot
kegg_80comb_df <- bind_rows(kegg_80dn_df, kegg_80up_df) %>%
  mutate(
    Direction = factor(Direction, levels = c("DOWN", "UP")),
    Label     = fct_reorder(Label, GeneRatio_num, .fun = max)
  )

p_kegg_80comb <- kegg_80comb_df %>%
  ggplot(aes(x = GeneRatio_num, y = Label, color = p.adjust)) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(3, 8), breaks = c(3, 4, 5, 6, 7, 8)) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_grid(. ~ Direction) +
  labs(
    title = "KEGG Enrichment: 80 dpi Upregulated and Downregulated (Cov vs. Mock)",
    x     = "Gene Ratio",
    y     = NULL,
    color = "p.adjust",
    size  = "Gene Count"
  ) +
  theme_light() +
  theme(
    axis.text.y      = element_text(color = "black", size = 9),
    panel.border     = element_rect(color = "black"),
    strip.background = element_rect(color = "black"),
    panel.spacing    = unit(0.8, "lines")
  )

ggsave("figures/kegg_80dpi_up_down.png", plot = p_kegg_80comb, width = 10.5, height = 6, dpi = 300)


# ==============================================================================
# 5. 80dpi vs 4dpi cov
# ==============================================================================
dat_80vs4 <- read_csv("data/infectcov_dpi80.csv", show_col_types = FALSE) %>%
  rename(Gene = 1) %>%
  left_join(kegg_map, by = "Gene") %>%
  mutate(
    entrez_id = if_else(
      is.na(entrez_id) & str_detect(Gene, "^LOC[0-9]+$"),
      str_remove(Gene, "^LOC"),
      entrez_id
    )
  ) %>%
  filter(!is.na(entrez_id))

up_80vs4 <- dat_80vs4 %>%
  filter(log2FoldChange >= 1, padj <= 0.05) %>%
  pull(entrez_id)

dn_80vs4 <- dat_80vs4 %>%
  filter(log2FoldChange <= -1, padj <= 0.05) %>%
  pull(entrez_id)

kegg_80vs4_up <- enrichKEGG(
  gene         = up_80vs4,
  organism     = "maua",
  keyType      = "ncbi-geneid",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

kegg_80vs4_dn <- enrichKEGG(
  gene         = dn_80vs4,
  organism     = "maua",
  keyType      = "ncbi-geneid",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

write_csv(as.data.frame(kegg_80vs4_up), "data/kegg_80dpi_vs_4dpi_cov_up.csv")
write_csv(as.data.frame(kegg_80vs4_dn), "data/kegg_80dpi_vs_4dpi_cov_dn.csv")

# 80dpi vs 4dpi cov UP plot
kegg_80vs4_up_df <- as.data.frame(kegg_80vs4_up) %>%
  filter(category != "Human Diseases") %>%
  head(15) %>%
  mutate(
    Direction     = "UP",
    GeneRatio_num = Count / DESize,
    Label         = str_wrap(Description, 45)
  )

p_kegg_80vs4_up <- kegg_80vs4_up_df %>%
  mutate(Label = fct_reorder(Label, GeneRatio_num)) %>%
  ggplot(aes(x = GeneRatio_num, y = Label, color = p.adjust)) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(3, 8), breaks = c(12, 14, 16, 18, 20)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    title = "KEGG Enrichment: 80dpi vs 4dpi cov Upregulated",
    x     = "Gene Ratio",
    y     = NULL,
    color = "p.adjust",
    size  = "Gene Count"
  ) +
  theme_light() +
  theme(
    axis.text.y  = element_text(color = "black", size = 9.5),
    panel.border = element_rect(color = "black")
  )

ggsave("figures/kegg_80dpi_vs_4dpi_cov_up.png", plot = p_kegg_80vs4_up, width = 9, height = 3.6, dpi = 300)

# 80dpi vs 4dpi cov DOWN plot
kegg_80vs4_dn_df <- as.data.frame(kegg_80vs4_dn) %>%
  filter(category != "Human Diseases") %>%
  head(15) %>%
  mutate(
    Direction     = "DOWN",
    GeneRatio_num = Count / DESize,
    Label         = str_wrap(Description, 45)
  )

p_kegg_80vs4_dn <- kegg_80vs4_dn_df %>%
  mutate(Label = fct_reorder(Label, GeneRatio_num)) %>%
  ggplot(aes(x = GeneRatio_num, y = Label, color = p.adjust)) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(3, 8)) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    title = "KEGG Enrichment: 80dpi vs 4dpi cov Downregulated",
    x     = "Gene Ratio",
    y     = NULL,
    color = "p.adjust",
    size  = "Gene Count"
  ) +
  theme_light() +
  theme(
    axis.text.y  = element_text(color = "black", size = 9.5),
    panel.border = element_rect(color = "black")
  )

ggsave("figures/kegg_80dpi_vs_4dpi_cov_dn.png", plot = p_kegg_80vs4_dn, width = 9, height = 6.3, dpi = 300)

# 80dpi vs 4dpi cov combined UP & DOWN plot
kegg_80vs4_comb_df <- bind_rows(kegg_80vs4_dn_df, kegg_80vs4_up_df) %>%
  mutate(
    Direction = factor(Direction, levels = c("DOWN", "UP")),
    Label     = fct_reorder(Label, GeneRatio_num, .fun = max)
  )

p_kegg_80vs4_comb <- kegg_80vs4_comb_df %>%
  ggplot(aes(x = GeneRatio_num, y = Label, color = p.adjust)) +
  geom_point(aes(size = Count)) +
  scale_size(range = c(3, 8)) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_grid(. ~ Direction) +
  labs(
    title = "KEGG Enrichment: 80dpi vs 4dpi cov Upregulated and Downregulated",
    x     = "Gene Ratio",
    y     = NULL,
    color = "p.adjust",
    size  = "Gene Count"
  ) +
  theme_light() +
  theme(
    axis.text.y      = element_text(color = "black", size = 9),
    panel.border     = element_rect(color = "black"),
    strip.background = element_rect(color = "black"),
    panel.spacing    = unit(0.8, "lines")
  )

ggsave("figures/kegg_80dpi_vs_4dpi_cov_up_down.png", plot = p_kegg_80vs4_comb, width = 10.5, height = 7.5, dpi = 300)
