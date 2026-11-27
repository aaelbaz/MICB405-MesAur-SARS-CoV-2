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

# 1)
#==========================
#Filter
EXCL_OLFACTORY_MUCOSAL <- c(
  "Bpifa1", "Bpifb4", "LOC101823870", "LOC101827556", "LOC121138580", "LOC101843093",
  "LOC101843357", "LOC121140164", "LOC101826747", "LOC101842056", "LOC101835304",
  "Umodl1", "LOC101823275", "LOC101842827", "LOC101826429", "Scgb1c1", "Reg3gamma",
  "Gpx6", "LOC106021050", "Gfy", "Muc5b", "Amelx", "Ambn", "Fetub", "CUNH6orf58",
  "LOC101843376", "LOC110342735", "LOC106020982", "LOC101842428", "LOC101827565"
)
# 2) 
EXCL_ERYTHROCYTE <- c(
  "LOC121144049", "LOC121141272", "LOC101833678", "LOC121141276", "LOC121144051", "Alas2"
)
# 3) 
EXCL_MISC_UNCHAR <- c(
  "LOC101831512", "CUNH20orf173", "Smim38", "LOC106021361",
  "LOC121141649", "LOC101838074", "LOC110339469"
)

# 4) 
EXCL_UNCHAR_LOC <- c(
  "LOC121139129", "LOC121140554",
  "LOC106022148", "LOC121140390", "LOC106021785"
)
# 5)
EXCL_GENES <- c(
  EXCL_OLFACTORY_MUCOSAL,
  EXCL_ERYTHROCYTE,
  EXCL_MISC_UNCHAR,
  EXCL_UNCHAR_LOC
)
#==========================
# Sort
  rank_by_xiao_pi <- function(csv_path, basemean_min = 15, top_n = 5) {
    raw_df <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(raw_df)[1] <- "Gene"
    clean_df <- as_tibble(raw_df) |>
      mutate(
        baseMean       = as.numeric(baseMean),
        log2FoldChange = as.numeric(log2FoldChange),
        lfcSE          = as.numeric(lfcSE),
        padj           = as.numeric(padj),
        xiao_pi_score  = abs(log2FoldChange) * (-log10(pmax(padj, 1e-300)))
      ) |>
      filter(
        !Gene %in% EXCL_GENES,
        !startsWith(Gene, "unassigned"),
        baseMean >= basemean_min
      )
    top_up <- clean_df |>
      filter(log2FoldChange > 0) |>
      arrange(desc(xiao_pi_score)) |>
      slice_head(n = top_n) |>
      mutate(Direction = "UP", Rank = row_number())
    top_down <- clean_df |>
      filter(log2FoldChange < 0) |>
      arrange(desc(xiao_pi_score)) |>
      slice_head(n = top_n) |>
      mutate(Direction = "DOWN", Rank = row_number())
    bind_rows(top_up, top_down) |>
      select(Gene, baseMean, log2FoldChange, lfcSE, padj, xiao_pi_score)
  }
# ============================
# Data
top5_4dpi <- rank_by_xiao_pi(
  "data/4dpi_results_filtered.csv",
  basemean_min = 15, top_n = 5
)
top5_80dpi <- rank_by_xiao_pi(
  "data/80dpi_results_filtered.csv",
  basemean_min = 10, top_n = 5
)
top5_int_updated <- rank_by_xiao_pi(
  "data/infectcov_dpi80.csv",
  basemean_min = 20, top_n = 5
)

dir.create("data/filtered and sorted", recursive = TRUE, showWarnings = FALSE)
write_csv(top5_4dpi,        "data/filtered and sorted/4dpi_top5_up_down.csv")
write_csv(top5_80dpi,       "data/filtered and sorted/80dpi_top5_up_down.csv")
write_csv(top5_int_updated, "data/filtered and sorted/80dpi_vs_4dpi_cov_top5_up_down.csv")