# Running DESeq2
# MICB 405 101 (2026W1)


# Loading required packages ---------------------------------------------------

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)


# Global Variable Definition ---------------------------------------------------

sam_metadata <- data.frame()
sam_count_matrix <- data.frame()

rep_metadata <- data.frame()
rep_count_matrix <- data.frame()


# Helper Function for Loading Data (both by sample & by replicate) -------------

loadData <- function(dpi_dir_path, dpi) {
  for (dir in dpi_dir_path) {
    sam_info <- str_split(basename(dir), pattern = "-")[[1]]
    smeta <- str_split(sam_info[1], pattern = "_")[[1]]
    sname <- paste(smeta[1], sam_info[2], sep = "-")
    smeta_row <- c(sname, dpi, smeta[2:length(smeta)])
    sam_metadata <<- rbind(sam_metadata, smeta_row)
    
    rep_files <- list.files(path = dir, pattern = "*.readcounts.out.tab", full.names = TRUE)
    sam_counts <- data.frame()
    for (rep in rep_files) {
      rname <- sub(".readcounts.out.tab", "", basename(rep))
      rdata <- read.delim(rep, header = FALSE, skip = 4, col.names = c("gene_id", "total", "antisense", "sense")) %>%
        dplyr::select(gene_id, sense) %>%
        dplyr::rename(!!rname := sense)
      if (length(sam_counts) == 0) {
        sam_counts <- rdata
      } else {
        sam_counts <- inner_join(sam_counts, rdata, by = "gene_id")
      }
      
      rmeta_row <- c(rname, smeta_row)
      rep_metadata <<- rbind(rep_metadata, rmeta_row)
      if (length(rep_count_matrix) == 0) {
        rep_count_matrix <<- rdata
      } else {
        rep_count_matrix <<- inner_join(rep_count_matrix, rdata, by = "gene_id")
      }
    }
    sam_counts <- sam_counts %>%
      rowwise(gene_id) %>%
      summarize(!!sname := sum(c_across(where(is.numeric))))
    if (length(sam_count_matrix) == 0) {
      sam_count_matrix <<- sam_counts
    } else {
      sam_count_matrix <<- inner_join(sam_count_matrix, sam_counts, by = "gene_id")
    }
  }
  colnames(sam_metadata) <<- c("sname", "dpi", "sex", "infect")
  colnames(rep_metadata) <<- c("rname", "sname", "dpi", "sex", "infect") 
}


# Loading data FROM SCRATCH ----------------------------------------------------

dpi4_samples <- list.files(path = "../gene_counts/4dpi_samples", full.names = TRUE)
dpi80_samples <- list.files(path = "../gene_counts/80dpi_samples", full.names = TRUE)

loadData(dpi4_samples, "4")
loadData(dpi80_samples, "80")


# Saving Data for Faster Loading -----------------------------------------------
write_tsv(sam_metadata, file = "../DESeq2_files/sam_input/sam_metadata.tsv")
write_tsv(sam_count_matrix, file = "../DESeq2_files/sam_input/sam_count_matrix.tsv")

write_tsv(rep_metadata, file = "../DESeq2_files/rep_input/rep_metadata.tsv")
write_tsv(rep_count_matrix, file = "../DESeq2_files/rep_input/rep_count_matrix.tsv")







# Starting from files (above has previously been run) --------------------------

# Loading required packages ----------------------------------------------------

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)





# Workflow by SAMPLE -----------------------------------------------------------

# Loading Data
sam_metadata <- read_tsv("../DESeq2_files/sam_input/sam_metadata.tsv")
sam_count_matrix <- read_tsv("../DESeq2_files/sam_input/sam_count_matrix.tsv")


# Coerce count_matrix & metadata into appropriate form & prep for DESeq2
sam_count_matrix <- sam_count_matrix %>% column_to_rownames(var = "gene_id")
sam_count_matrix <- as.matrix(sam_count_matrix)

sam_metadata <- sam_metadata %>% column_to_rownames(var = "sname")
sam_metadata$dpi <- as.factor(sam_metadata$dpi)
sam_metadata$sex <- as.factor(sam_metadata$sex)
sam_metadata$infect <- as.factor(sam_metadata$infect)


# Create our DESeq2 object.
sam_dds_matrix <- DESeqDataSetFromMatrix(countData = sam_count_matrix, 
                                         colData = sam_metadata,
                                         design = ~dpi+sex+infect)


# Set control condition using the relevel function.
sam_dds_matrix$dpi <- relevel(sam_dds_matrix$dpi, ref = "4")
sam_dds_matrix$infect <- relevel(sam_dds_matrix$infect, ref = "mock")
sam_dds_matrix$sex <- relevel(sam_dds_matrix$sex, ref = "M")


# Run DESeq2
sam_dds <- DESeq(sam_dds_matrix)



# Sanity Check 1: Make a PCA plot with the data

# Define helper function for PCA plotting & saving
sam_plot_and_save_PCA <- function(rld, group, name) {
  PCA_plot <- plotPCA(rld, intgroup = group)
  ggsave(paste("../DESeq2_files/sam_output/figures", name, sep = "/"),
         plot = PCA_plot, width = 10, height = 10)
  PCA_plot
}

# Perform log transformation on our count data
sam_rld <- vst(sam_dds)


# Generate a PCA plot with DESeq2's plotPCA function (via helper)
sam_plot_and_save_PCA(sam_rld, "dpi", "pca_dpi.png")
sam_plot_and_save_PCA(sam_rld, "infect", "pca_infect.png")
sam_plot_and_save_PCA(sam_rld, "sex", "pca_sex.png")
sam_plot_and_save_PCA(sam_rld, c("dpi", "infect"), "pca_dpi_infect.png")
sam_plot_and_save_PCA(sam_rld, c("dpi", "sex"), "pca_dpi_sex.png")
sam_plot_and_save_PCA(sam_rld, c("infect", "sex"), "pca_infect_sex.png")
sam_plot_and_save_PCA(sam_rld, c("dpi", "infect", "sex"), "pca_dpi_infect_sex.png")



# Sanity Check 2: Make a distance matrix with the data

# Calculate distances between samples in our log-transformed data
sam_dists <- dist(t(assay(sam_rld)))

# Convert the output to a matrix
sam_dist_matrix <- as.matrix(sam_dists)

# Remove the column names of our matrix
colnames(sam_dist_matrix) <- NULL

# Set the colour palette for our heatmap
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Generate a heatmap using the pheatmap package
sam_heatmap <- pheatmap(sam_dist_matrix,
                        clustering_distance_rows = sam_dists,
                        clustering_distance_cols = sam_dists, 
                        col = colours)

ggsave("../DESeq2_files/sam_output/figures/heatmap.png", plot = sam_heatmap, width = 12, height = 12)


# Further analyses of differentially expressed genes

# names of the results that DESeq2 calculated
sam_res_names <- resultsNames(sam_dds)

# Get results for each comparison
sam_res_dpi <- results(sam_dds, name = sam_res_names[2]) %>% as.data.frame()
sam_res_sex <- results(sam_dds, name = sam_res_names[3]) %>% as.data.frame()
sam_res_infect <- results(sam_dds, name = sam_res_names[4]) %>% as.data.frame()

# Getting rid of NA values
sam_res_dpi_no_NA <- sam_res_dpi %>% drop_na()
sam_res_sex_no_NA <- sam_res_sex %>% drop_na()
sam_res_infect_no_NA <- sam_res_infect %>% drop_na()

# Filtering for adjusted p-values <=0.05
sam_res_dpi_filtered <- sam_res_dpi_no_NA %>% filter(padj <= 0.05)
sam_res_sex_filtered <- sam_res_sex_no_NA %>% filter(padj <= 0.05)
sam_res_infect_filtered <- sam_res_infect_no_NA %>% filter(padj <= 0.05)


# Pull genes with more than 2x higher/lower expression
sam_res_dpi_final <- sam_res_dpi_filtered %>% filter(abs(log2FoldChange) >= 1)
sam_res_sex_final <- sam_res_sex_filtered %>% filter(abs(log2FoldChange) >= 1)
sam_res_infect_final <- sam_res_infect_filtered %>% filter(abs(log2FoldChange) >= 1)

# Top up-regulated differentially expressed genes.
# dpi genes are top 20 given the number of differentially expressed.
# sex and infect are all up-regulated genes given lower number of differentially expressed.
sam_top_dpi_genes <- sam_res_dpi_final %>%
  filter(log2FoldChange >= 1) %>%
  arrange(desc(log2FoldChange)) %>%
  head(n = 20)

sam_top_sex_genes <- sam_res_sex_final %>%
  filter(log2FoldChange >= 1) %>%
  arrange(desc(log2FoldChange))

sam_top_infect_genes <- sam_res_infect_final %>%
  filter(log2FoldChange >= 1) %>%
  arrange(desc(log2FoldChange))


# Top down-regulated differentially expressed genes. 
# dpi genes are top 20 given the number of differentially expressed.
# sex and infect are all down-regulated genes given lower number of differentially expressed.
sam_bot_dpi_genes <- sam_res_dpi_final %>%
  filter(log2FoldChange <= -1) %>%
  arrange(log2FoldChange) %>%
  head(n = 20)

sam_bot_sex_genes <- sam_res_sex_final %>%
  filter(log2FoldChange <= -1) %>%
  arrange(log2FoldChange)

sam_bot_infect_genes <- sam_res_infect_final %>%
  filter(log2FoldChange <= -1) %>%
  arrange(log2FoldChange)


# Write Results to file
write_csv(sam_res_dpi_final %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/sam_output/data/all_dpi_diff_expressed.csv")
write_csv(sam_res_sex_final %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/sam_output/data/all_sex_diff_expressed.csv")
write_csv(sam_res_infect_final %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/sam_output/data/all_infect_diff_expressed.csv")

write_csv(sam_top_dpi_genes %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/sam_output/data/top_dpi_diff_expressed.csv")
write_csv(sam_top_sex_genes %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/sam_output/data/top_sex_diff_expressed.csv")
write_csv(sam_top_infect_genes %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/sam_output/data/top_infect_diff_expressed.csv")

write_csv(sam_bot_dpi_genes %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/sam_output/data/bot_dpi_diff_expressed.csv")
write_csv(sam_bot_sex_genes %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/sam_output/data/bot_sex_diff_expressed.csv")
write_csv(sam_bot_infect_genes %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/sam_output/data/bot_infect_diff_expressed.csv")







# Workflow by TECHNICAL REPLICATE ----------------------------------------------
# Loading Data
rep_metadata <- read_tsv("../DESeq2_files/rep_input/rep_metadata.tsv")
rep_count_matrix <- read_tsv("../DESeq2_files/rep_input/rep_count_matrix.tsv")


# Coerce count_matrix & metadata into appropriate form & prep for DESeq2
rep_count_matrix <- rep_count_matrix %>% column_to_rownames(var = "gene_id")
rep_count_matrix <- as.matrix(rep_count_matrix)

rep_metadata <- rep_metadata %>% column_to_rownames(var = "rname")
rep_metadata$sname <- as.factor(rep_metadata$sname)
rep_metadata$dpi <- as.factor(rep_metadata$dpi)
rep_metadata$sex <- as.factor(rep_metadata$sex)
rep_metadata$infect <- as.factor(rep_metadata$infect)


# Create our FIRST DESeq2 object (to acquire batch effects)
rep_dds_matrix <- DESeqDataSetFromMatrix(countData = rep_count_matrix,
                                         colData = rep_metadata,
                                         design = ~dpi+sex+infect)


# Set control condition using the relevel function.
rep_dds_matrix$dpi <- relevel(rep_dds_matrix$dpi, ref = "4")
rep_dds_matrix$infect <- relevel(rep_dds_matrix$infect, ref = "mock")
rep_dds_matrix$sex <- relevel(rep_dds_matrix$sex, ref = "M")

rep_dds_matrix$sname <- relevel(rep_dds_matrix$sname, ref = "G1-11_S2")


# Run DESeq2
rep_dds <- DESeq(rep_dds_matrix)



# Sanity Check 1: Make a PCA plot with the data

# Define helper function for PCA plotting & saving
rep_plot_and_save_PCA <- function(rld, group, name) {
  PCA_plot <- plotPCA(rld, intgroup = group)
  ggsave(paste("../DESeq2_files/rep_output/figures", name, sep = "/"), plot = PCA_plot, width = 10, height = 10)
  PCA_plot
}

# Perform log transformation on our count data
rep_rld <- vst(rep_dds)


# Generate a PCA plot with DESeq2's plotPCA function (via helper)
rep_plot_and_save_PCA(rep_rld, "dpi", "without_pc1/pca_dpi.png")
rep_plot_and_save_PCA(rep_rld, "infect", "without_pc1/pca_infect.png")
rep_plot_and_save_PCA(rep_rld, "sex", "without_pc1/pca_sex.png")
rep_plot_and_save_PCA(rep_rld, c("dpi", "infect"), "without_pc1/pca_dpi_infect.png")
rep_plot_and_save_PCA(rep_rld, c("dpi", "sex"), "without_pc1/pca_dpi_sex.png")
rep_plot_and_save_PCA(rep_rld, c("infect", "sex"), "without_pc1/pca_infect_sex.png")
rep_plot_and_save_PCA(rep_rld, c("dpi", "infect", "sex"), "without_pc1/pca_dpi_infect_sex.png")

rep_plot_and_save_PCA(rep_rld, "sname", "without_pc1/pca_sample.png")


# Get PCA data and bind PC1 to metadata
rep_pca_data <- plotPCA(rep_rld, intgroup = "dpi", returnData = TRUE) %>%
  mutate(batch_pc1 = PC1) %>%
  select(batch_pc1, name)

rep_metadata_alt <- rep_metadata %>%
  rownames_to_column(var = "name") %>%
  inner_join(rep_pca_data, by = "name") %>%
  column_to_rownames(var = "name")


# Center & scale PC1 values for use in DESeq2 (rejects numeric columns with large means)
rep_metadata_alt <- rep_metadata_alt %>%
  mutate(batch_pc1 = scale(batch_pc1, center = TRUE, scale = TRUE))


# Create our SECOND DESeq2 object (accounting for batch effects)
rep_dds_matrix_alt <- DESeqDataSetFromMatrix(countData = rep_count_matrix,
                                             colData = rep_metadata_alt,
                                             design = ~dpi+sex+infect+batch_pc1)


# Set control condition using the relevel function.
rep_dds_matrix_alt$dpi <- relevel(rep_dds_matrix_alt$dpi, ref = "4")
rep_dds_matrix_alt$infect <- relevel(rep_dds_matrix_alt$infect, ref = "mock")
rep_dds_matrix_alt$sex <- relevel(rep_dds_matrix_alt$sex, ref = "M")
rep_dds_matrix_alt$sname <- relevel(rep_dds_matrix_alt$sname, ref = "G1-11_S2")
# Not sure what we would relevel PC1 to considering it is not a factor...


# Run DESeq2
rep_dds_alt <- DESeq(rep_dds_matrix_alt)


# Perform log transformation on our count data
rep_rld_alt <- vst(rep_dds_alt)


# Generate a PCA plot with DESeq2's plotPCA function (via helper)
rep_plot_and_save_PCA(rep_rld_alt, "dpi", "with_pc1/pca_dpi.png")
rep_plot_and_save_PCA(rep_rld_alt, "infect", "with_pc1/pca_infect.png")
rep_plot_and_save_PCA(rep_rld_alt, "sex", "with_pc1/pca_sex.png")
rep_plot_and_save_PCA(rep_rld_alt, c("dpi", "infect"), "with_pc1/pca_dpi_infect.png")
rep_plot_and_save_PCA(rep_rld_alt, c("dpi", "sex"), "with_pc1/pca_dpi_sex.png")
rep_plot_and_save_PCA(rep_rld_alt, c("infect", "sex"), "with_pc1/pca_infect_sex.png")
rep_plot_and_save_PCA(rep_rld_alt, c("dpi", "infect", "sex"), "with_pc1/pca_dpi_infect_sex.png")

rep_plot_and_save_PCA(rep_rld_alt, "sname", "with_pc1/pca_sample.png")









### BELOW NOT COMPLETED YET!!!

# Sanity Check 2: Make a distance matrix with the data --------------------

# Calculate distances between samples in our log-transformed data
sample_dists <- dist(t(assay(rld)))

# Convert the output to a matrix
sample_dist_matrix <- as.matrix(sample_dists)

# Remove the column names of our matrix
colnames(sample_dist_matrix) <- NULL

# Set the colour palette for our heatmap
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Generate a heatmap using the pheatmap package
heatmap <- pheatmap(sample_dist_matrix,
                    clustering_distance_rows = sample_dists,
                    clustering_distance_cols = sample_dists, 
                    col = colours)

ggsave("../DESeq2_files/rep_output/heatmap.png", plot = heatmap, width = 12, height = 12)


# Further analyses of differentially expressed genes ----------------------

# names of the results that DESeq2 calculated
result_names <- resultsNames(dds)

# Get results for each comparison
res_dpi <- results(dds, name = result_names[2]) %>% as.data.frame()
res_sex <- results(dds, name = result_names[3]) %>% as.data.frame()
res_infect <- results(dds, name = result_names[4]) %>% as.data.frame()

# Getting rid of NA values
res_dpi_no_NA <- res_dpi %>% drop_na()
res_sex_no_NA <- res_sex %>% drop_na()
res_infect_no_NA <- res_infect %>% drop_na()

# Filtering for adjusted p-values <=0.05
res_dpi_filtered <- res_dpi_no_NA %>% filter(padj <= 0.05)
res_sex_filtered <- res_sex_no_NA %>% filter(padj <= 0.05)
res_infect_filtered <- res_infect_no_NA %>% filter(padj <= 0.05)


# Pull genes with more than 2x higher/lower expression
res_dpi_final <- res_dpi_filtered %>% filter(abs(log2FoldChange) >= 1)
res_sex_final <- res_sex_filtered %>% filter(abs(log2FoldChange) >= 1)
res_infect_final <- res_infect_filtered %>% filter(abs(log2FoldChange) >= 1)

# Top up-regulated differentially expressed genes.
# dpi genes are top 20 given the number of differentially expressed.
# sex and infect are all up-regulated genes given lower number of differentially expressed.
top_dpi_genes <- res_dpi_final %>%
  filter(log2FoldChange >= 1) %>%
  arrange(desc(log2FoldChange)) %>%
  head(n = 20)

top_sex_genes <- res_sex_final %>%
  filter(log2FoldChange >= 1) %>%
  arrange(desc(log2FoldChange))

top_infect_genes <- res_infect_final %>%
  filter(log2FoldChange >= 1) %>%
  arrange(desc(log2FoldChange))


# Top down-regulated differentially expressed genes. 
# dpi genes are top 20 given the number of differentially expressed.
# sex and infect are all down-regulated genes given lower number of differentially expressed.
bot_dpi_genes <- res_dpi_final %>%
  filter(log2FoldChange <= -1) %>%
  arrange(log2FoldChange) %>%
  head(n = 20)

bot_sex_genes <- res_sex_final %>%
  filter(log2FoldChange <= -1) %>%
  arrange(log2FoldChange)

bot_infect_genes <- res_infect_final %>%
  filter(log2FoldChange <= -1) %>%
  arrange(log2FoldChange)


# Write Results to file
write_csv(res_dpi_final %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/rep_output/all_dpi_diff_expressed.csv")
write_csv(res_sex_final %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/rep_output/all_sex_diff_expressed.csv")
write_csv(res_infect_final %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/rep_output/all_infect_diff_expressed.csv")

write_csv(top_dpi_genes %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/rep_output/top_dpi_diff_expressed.csv")
write_csv(top_sex_genes %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/rep_output/top_sex_diff_expressed.csv")
write_csv(top_infect_genes %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/rep_output/top_infect_diff_expressed.csv")

write_csv(bot_dpi_genes %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/rep_output/bot_dpi_diff_expressed.csv")
write_csv(bot_sex_genes %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/rep_output/bot_sex_diff_expressed.csv")
write_csv(bot_infect_genes %>% rownames_to_column(var = "gene_id"), file = "../DESeq2_files/rep_output/bot_infect_diff_expressed.csv")





















# Aside to check sample membership of 4dpi technical replicates in PCA

dpi4_only_reps_meta <- replicate_metadata %>% filter(dpi == "4")
dpi4_only_rnames <- rownames(dpi4_only_reps_meta)

dpi4_only_reps_count_matrix <- as.data.frame(replicate_count_matrix) %>%
  select(all_of(dpi4_only_rnames)) %>%
  as.matrix()

dpi4_only_dds_matrix <- DESeqDataSetFromMatrix(countData = dpi4_only_reps_count_matrix, 
                                               colData = dpi4_only_reps_meta,
                                               design = ~sex+infect)

dpi4_only_dds_matrix$infect <- relevel(dpi4_only_dds_matrix$infect, ref = "mock")
dpi4_only_dds_matrix$sex <- relevel(dpi4_only_dds_matrix$sex, ref = "M")

dpi4_only_dds <- DESeq(dpi4_only_dds_matrix)

dpi4_only_replicate_rld <- vst(dpi4_only_dds)

pca_data <- plotPCA(dpi4_only_replicate_rld, intgroup = c("sname"), returnData = TRUE)
ggplot(pca_data, aes(x = PC1, y = PC2, shape = sname, colour = sname)) +
  geom_point(size = 3) +
  #scale_y_continuous(limit = c(-10, 5)) +
  scale_shape_manual(values=1:16) +
  scale_color_manual(values = c("#5050FF", "#CE3D32", "#749B58", "#F0E685", "#466983", "#BA6338", "#5DB1DD", "#802268", "#6BD76B", "#D595A7", "#924822", "#837B8D", "#C75127", "#D58F5C", "#612A79", "#990033"))


pca_data_filt <- filter(pca_data, PC1 >= 0)
pca_data_filt$name


ggplot(pca_data_filt, aes(x = PC1, y = PC2, shape = sname, colour = sname)) +
  geom_point(size = 3) +
  #scale_y_continuous(limit = c(-10, 5)) +
  scale_shape_manual(values=1:16) +
  scale_color_manual(values = c("#5050FF", "#CE3D32", "#749B58", "#F0E685", "#466983", "#BA6338", "#5DB1DD", "#802268", "#6BD76B", "#D595A7", "#924822", "#837B8D", "#C75127", "#D58F5C", "#612A79", "#990033"))



ggsave("deseq2_files/pca_replicates_sname_split.png", plot = get_last_plot(), width = 20, height = 20)





#### BELOW THIS POINT THE ANALYSIS ONLY INCLUDES PER BIOLOGICAL SAMPLE ####
  # that is, no further analysis of differential expression among individual technical replicates





# Sanity Check 2: Make a distance matrix with the data --------------------

# Calculate distances between samples in our log-transformed data
sample_dists <- dist(t(assay(rld)))

# Convert the output to a matrix
sample_dist_matrix <- as.matrix(sample_dists)

# Remove the column names of our matrix
colnames(sample_dist_matrix) <- NULL

# Set the colour palette for our heatmap
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Generate a heatmap using the pheatmap package
heatmap <- pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists, 
         col = colours)

ggsave("deseq2_files/heatmap.png", plot = heatmap, width = 12, height = 12)


# Further analyses of differentially expressed genes ----------------------

# names of the results that DESeq2 calculated
result_names <- resultsNames(dds)

# Get results for each comparison
res_dpi <- results(dds, name = result_names[2]) %>% as.data.frame()
res_sex <- results(dds, name = result_names[3]) %>% as.data.frame()
res_infect <- results(dds, name = result_names[4]) %>% as.data.frame()

# Getting rid of NA values
res_dpi_no_NA <- res_dpi %>% drop_na()
res_sex_no_NA <- res_sex %>% drop_na()
res_infect_no_NA <- res_infect %>% drop_na()

# Filtering for adjusted p-values <=0.05
res_dpi_filtered <- res_dpi_no_NA %>% filter(padj <= 0.05)
res_sex_filtered <- res_sex_no_NA %>% filter(padj <= 0.05)
res_infect_filtered <- res_infect_no_NA %>% filter(padj <= 0.05)


# Pull genes with more than 2x higher/lower expression
res_dpi_final <- res_dpi_filtered %>% filter(abs(log2FoldChange) >= 1)
res_sex_final <- res_sex_filtered %>% filter(abs(log2FoldChange) >= 1)
res_infect_final <- res_infect_filtered %>% filter(abs(log2FoldChange) >= 1)

# Top up-regulated differentially expressed genes.
  # dpi genes are top 20 given the number of differentially expressed.
  # sex and infect are all up-regulated genes given lower number of differentially expressed.
top_dpi_genes <- res_dpi_final %>%
  filter(log2FoldChange >= 1) %>%
  arrange(desc(log2FoldChange)) %>%
  head(n = 20)

top_sex_genes <- res_sex_final %>%
  filter(log2FoldChange >= 1) %>%
  arrange(desc(log2FoldChange))

top_infect_genes <- res_infect_final %>%
  filter(log2FoldChange >= 1) %>%
  arrange(desc(log2FoldChange))


# Top down-regulated differentially expressed genes. 
  # dpi genes are top 20 given the number of differentially expressed.
  # sex and infect are all down-regulated genes given lower number of differentially expressed.
bot_dpi_genes <- res_dpi_final %>%
  filter(log2FoldChange <= -1) %>%
  arrange(log2FoldChange) %>%
  head(n = 20)

bot_sex_genes <- res_sex_final %>%
  filter(log2FoldChange <= -1) %>%
  arrange(log2FoldChange)

bot_infect_genes <- res_infect_final %>%
  filter(log2FoldChange <= -1) %>%
  arrange(log2FoldChange)


# Write Results to file
write_csv(res_dpi_final %>% rownames_to_column(var = "gene_id"), file = "deseq2_files/all_dpi_diff_expressed.csv")
write_csv(res_sex_final %>% rownames_to_column(var = "gene_id"), file = "deseq2_files/all_sex_diff_expressed.csv")
write_csv(res_infect_final %>% rownames_to_column(var = "gene_id"), file = "deseq2_files/all_infect_diff_expressed.csv")

write_csv(top_dpi_genes %>% rownames_to_column(var = "gene_id"), file = "deseq2_files/top_dpi_diff_expressed.csv")
write_csv(top_sex_genes %>% rownames_to_column(var = "gene_id"), file = "deseq2_files/top_sex_diff_expressed.csv")
write_csv(top_infect_genes %>% rownames_to_column(var = "gene_id"), file = "deseq2_files/top_infect_diff_expressed.csv")

write_csv(bot_dpi_genes %>% rownames_to_column(var = "gene_id"), file = "deseq2_files/bot_dpi_diff_expressed.csv")
write_csv(bot_sex_genes %>% rownames_to_column(var = "gene_id"), file = "deseq2_files/bot_sex_diff_expressed.csv")
write_csv(bot_infect_genes %>% rownames_to_column(var = "gene_id"), file = "deseq2_files/bot_infect_diff_expressed.csv")
