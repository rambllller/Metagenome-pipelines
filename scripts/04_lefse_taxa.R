#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) args[idx + 1] else default
}

input_file <- get_arg("--input")
metadata_file <- get_arg("--metadata")
outdir <- get_arg("--outdir")
group_col <- get_arg("--group-col", "Group")
prevalence_cutoff <- as.numeric(get_arg("--prevalence", "0.10"))
lda_cutoff <- as.numeric(get_arg("--lda", "2"))
kruskal_cutoff <- as.numeric(get_arg("--kruskal", "0.05"))
p_adjust_method <- get_arg("--padj", "none")

if (is.null(input_file) || is.null(metadata_file) || is.null(outdir)) {
  stop("Usage: Rscript 04_lefse_taxa.R --input taxa_abundance.tsv --metadata groups.tsv --outdir outdir [--group-col Group]")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

need_pkgs <- c("data.table", "dplyr", "stringr", "SummarizedExperiment", "S4Vectors", "lefser", "ggplot2")
for (pkg in need_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    if (pkg %in% c("SummarizedExperiment", "S4Vectors", "lefser")) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(lefser)
  library(ggplot2)
})

abund_raw <- fread(input_file, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
meta_raw  <- fread(metadata_file, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)

stopifnot("Taxon_Name" %in% colnames(abund_raw))
stopifnot("SampleID" %in% colnames(meta_raw))
stopifnot(group_col %in% colnames(meta_raw))

sample_cols <- setdiff(colnames(abund_raw), "Taxon_Name")
common_samples <- intersect(sample_cols, meta_raw$SampleID)
if (length(common_samples) < 2) stop("Too few matched samples")

meta_use <- meta_raw %>% filter(SampleID %in% common_samples)
meta_use[[group_col]] <- factor(meta_use[[group_col]])

abund_use <- abund_raw[, c("Taxon_Name", meta_use$SampleID), drop = FALSE]
feature_mat <- as.matrix(abund_use[, -1, drop = FALSE])
storage.mode(feature_mat) <- "numeric"
rownames(feature_mat) <- abund_use$Taxon_Name
feature_mat[!is.finite(feature_mat)] <- 0

feature_mat <- feature_mat[rowSums(feature_mat, na.rm = TRUE) > 0, , drop = FALSE]
prev <- rowMeans(feature_mat > 0, na.rm = TRUE)
feature_mat <- feature_mat[prev >= prevalence_cutoff, , drop = FALSE]
keep_var <- apply(feature_mat, 1, function(x) stats::var(x, na.rm = TRUE) > 0)
feature_mat <- feature_mat[keep_var, , drop = FALSE]

write.table(
  data.frame(Taxon_Name = rownames(feature_mat), feature_mat, check.names = FALSE),
  file = file.path(outdir, "01_filtered_taxa_matrix.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

coldata <- DataFrame(row.names = meta_use$SampleID)
coldata[[group_col]] <- meta_use[[group_col]]
rowdata <- DataFrame(Taxon_Name = rownames(feature_mat), row.names = rownames(feature_mat))
se <- SummarizedExperiment(assays = list(abundance = feature_mat), colData = coldata, rowData = rowdata)
se_relab <- relativeAb(se)

lefse_res <- lefser(
  relab = se_relab,
  classCol = group_col,
  kruskal.threshold = kruskal_cutoff,
  lda.threshold = lda_cutoff,
  method = p_adjust_method,
  trim.names = FALSE,
  checkAbundances = TRUE
)

write.table(lefse_res, file = file.path(outdir, "02_LEFSE_results_raw.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

if (nrow(lefse_res) == 0) {
  message("No significant taxa detected under current thresholds.")
  quit(save = "no")
}

class_levels <- levels(meta_use[[group_col]])
stopifnot(length(class_levels) == 2)
class0 <- class_levels[1]
class1 <- class_levels[2]

lefse_res2 <- lefse_res %>%
  mutate(
    Taxon_Name = features,
    Enriched_Group = ifelse(scores > 0, class1, class0),
    Abs_LDA = abs(scores)
  ) %>%
  arrange(desc(Abs_LDA))

write.table(lefse_res2, file = file.path(outdir, "03_LEFSE_results_with_taxon.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

plot_res <- lefse_res %>% mutate(features = stringr::str_wrap(features, width = 50))
p1 <- suppressWarnings(lefserPlot(plot_res, trim.names = FALSE, title = "LEfSe of taxonomic features")) +
  scale_fill_discrete(name = group_col, labels = c(class0, class1))

ggsave(file.path(outdir, "04_LEFSE_LDA_barplot.pdf"), plot = p1, width = 12, height = max(5, 0.25 * nrow(lefse_res2) + 2))
ggsave(file.path(outdir, "04_LEFSE_LDA_barplot.png"), plot = p1, width = 12, height = max(5, 0.25 * nrow(lefse_res2) + 2), dpi = 300)

plot_df <- lefse_res2 %>%
  slice_head(n = min(30, nrow(lefse_res2))) %>%
  arrange(Abs_LDA) %>%
  mutate(Taxon_Name = factor(stringr::str_wrap(Taxon_Name, 50), levels = stringr::str_wrap(Taxon_Name, 50)))

p2 <- ggplot(plot_df, aes(x = Taxon_Name, y = Abs_LDA, fill = Enriched_Group)) +
  geom_col(width = 0.8) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(title = "Top differential taxa by LEfSe", x = NULL, y = "Absolute LDA score (log10)")

ggsave(file.path(outdir, "05_LEFSE_top30_custom_barplot.pdf"), plot = p2, width = 12, height = 8)
ggsave(file.path(outdir, "05_LEFSE_top30_custom_barplot.png"), plot = p2, width = 12, height = 8, dpi = 300)
