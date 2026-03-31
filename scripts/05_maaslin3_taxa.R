#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) args[idx + 1] else default
}

input_file <- get_arg("--input")
metadata_file <- get_arg("--metadata")
outdir <- get_arg("--outdir")
fixed_effects_str <- get_arg("--fixed-effects", "Group")
reference <- get_arg("--reference", "Group,NO")
max_significance <- as.numeric(get_arg("--max-significance", "0.2"))
min_prevalence <- as.numeric(get_arg("--min-prevalence", "0.1"))
min_abundance <- as.numeric(get_arg("--min-abundance", "0"))
evaluate_only <- get_arg("--evaluate-only", "abundance")

if (is.null(input_file) || is.null(metadata_file) || is.null(outdir)) {
  stop("Usage: Rscript 05_maaslin3_taxa.R --input taxa_abundance.tsv --metadata groups.tsv --outdir outdir [--fixed-effects Group] [--reference Group,NO]")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("maaslin3", quietly = TRUE)) BiocManager::install("biobakery/maaslin3", ask = FALSE, update = FALSE)
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

suppressPackageStartupMessages({
  library(maaslin3)
  library(data.table)
  library(dplyr)
})

abund_raw <- fread(input_file, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
meta_raw  <- fread(metadata_file, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)

stopifnot("Taxon_Name" %in% colnames(abund_raw))
stopifnot("SampleID" %in% colnames(meta_raw))

sample_cols <- setdiff(colnames(abund_raw), "Taxon_Name")
common_samples <- intersect(sample_cols, meta_raw$SampleID)
if (length(common_samples) < 2) stop("Too few matched samples")

meta_use <- meta_raw[match(common_samples, meta_raw$SampleID), , drop = FALSE]
rownames(meta_use) <- meta_use$SampleID
meta_use$SampleID <- NULL

feat <- t(as.matrix(abund_raw[, common_samples, drop = FALSE]))
colnames(feat) <- abund_raw$Taxon_Name
feat <- as.data.frame(feat, check.names = FALSE)
storage.mode(feat) <- "numeric"

fixed_effects <- trimws(strsplit(fixed_effects_str, ",")[[1]])
missing_effects <- setdiff(fixed_effects, colnames(meta_use))
if (length(missing_effects) > 0) stop(paste("Missing metadata columns:", paste(missing_effects, collapse = ", ")))

# Set factor reference if the reference variable exists in metadata
if (!is.null(reference)) {
  ref_parts <- strsplit(reference, ";")[[1]]
  for (rp in ref_parts) {
    x <- strsplit(rp, ",")[[1]]
    if (length(x) == 2 && x[1] %in% colnames(meta_use)) {
      meta_use[[x[1]]] <- factor(meta_use[[x[1]]])
      if (x[2] %in% levels(meta_use[[x[1]]])) {
        meta_use[[x[1]]] <- stats::relevel(meta_use[[x[1]]], ref = x[2])
      }
    }
  }
}

fit_out <- maaslin3(
  input_data = feat,
  input_metadata = meta_use,
  output = outdir,
  fixed_effects = fixed_effects,
  reference = reference,
  normalization = "TSS",
  transform = "LOG",
  standardize = TRUE,
  max_significance = max_significance,
  min_prevalence = min_prevalence,
  min_abundance = min_abundance,
  evaluate_only = evaluate_only,
  augment = TRUE
)

res_file <- file.path(outdir, "all_results.tsv")
res <- fread(res_file, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
write.table(res, file = file.path(outdir, "01_all_results_with_taxon.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

res_abund <- res %>% filter(model == "abundance")
write.table(res_abund, file = file.path(outdir, "02_abundance_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

if ("qval_individual" %in% colnames(res_abund)) {
  res_sig <- res_abund %>% filter(!is.na(qval_individual) & qval_individual <= max_significance)
  if (nrow(res_sig) > 0 && "coef" %in% colnames(res_sig)) {
    res_sig$Direction <- ifelse(res_sig$coef > 0, "Up_in_test_level", ifelse(res_sig$coef < 0, "Down_in_test_level", "No_change"))
  }
  write.table(res_sig, file = file.path(outdir, "03_significant_abundance_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}
