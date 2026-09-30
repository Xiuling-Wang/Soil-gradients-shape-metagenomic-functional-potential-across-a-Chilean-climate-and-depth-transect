#!/usr/bin/env Rscript
# Exploratory KO-level differential-abundance analysis for the Chile iDNA metagenomes.
#
# Scientific object: raw mapped DNA reads aggregated to fractionally assigned KOs.
# This is differential abundance / encoded functional potential, not gene expression.
#
# Two prespecified exploratory questions:
#   1. Surface sites: AZ (n = 2), SG (n = 3), LC (n = 3), NA surface (n = 3).
#   2. NA depths: surface, middle, and deep samples from three mapped pits.
#      The primary exploratory comparison is among the three depth groups;
#      pit is used only in a separate sensitivity-adjusted model.
#
# The script is additive: it never changes historical outputs, figures, raw data,
# or manuscript files. It writes a standalone exploratory result bundle.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(edgeR)
})

options(scipen = 999)
set.seed(42)
invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "UTF-8"), silent = TRUE)))

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run with Rscript so the project root can be resolved from this script.")
}
SCRIPT <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
ROOT <- normalizePath(file.path(dirname(SCRIPT), ".."), mustWork = TRUE)
OUT <- file.path(ROOT, "analysis", "exploratory_KO_differential_abundance_2026-08-29")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
INTERMEDIATE <- file.path(OUT, "intermediate")

SAMPLES <- c("C01", sprintf("C%02d", 3:18))
SURFACE_SITES <- c("AZ", "SG", "LC", "NB")
DEPTH_LEVELS <- c("5_10cm", "20_40cm", "60_80cm")

PATHS <- c(
  counts = file.path(ROOT, "data/processed/metagenome_output/Genecatalog/Nmapped_reads.tsv"),
  annotation = file.path(ROOT, "data/processed/metagenome_output/Genecatalog/eggNog.tsv"),
  metadata = file.path(ROOT, "data/metadata/env_2024_meta.tsv")
)
missing_inputs <- PATHS[!file.exists(PATHS)]
if (length(missing_inputs)) {
  stop("Missing required input(s):\n", paste(missing_inputs, collapse = "\n"))
}

write_tsv_checked <- function(x, path) {
  readr::write_tsv(x, path, na = "NA")
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Output was not written: ", path)
  }
  invisible(path)
}

first_informative <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "-"]
  if (length(x)) x[[1]] else NA_character_
}

safe_label <- function(ko, name) {
  name <- ifelse(is.na(name) | name == "", "uncharacterized", name)
  str_trunc(paste0(ko, " | ", name), width = 52, side = "right")
}

message("Reading and validating metadata...")
metadata <- readr::read_tsv(PATHS[["metadata"]], show_col_types = FALSE, na = "") %>%
  filter(sample_id %in% SAMPLES) %>%
  arrange(match(sample_id, SAMPLES)) %>%
  mutate(
    site = factor(site, levels = SURFACE_SITES),
    depth_label = factor(depth_label, levels = DEPTH_LEVELS),
    site_display = recode(as.character(site), NB = "NA"),
    depth_display = recode(
      as.character(depth_label),
      `5_10cm` = "surface",
      `20_40cm` = "middle",
      `60_80cm` = "deep"
    ),
    pit = factor(site_pit)
  )

if (!identical(metadata$sample_id, SAMPLES)) {
  stop("Metadata sample order does not match the 17-sample metagenome universe.")
}
if (!all(metadata$dna_type == "iDNA")) {
  stop("All metagenome samples must be iDNA; check metadata before modelling.")
}

surface_meta <- metadata %>%
  filter(depth_label == "5_10cm") %>%
  mutate(site_display = factor(site_display, levels = c("AZ", "SG", "LC", "NA")))
if (!identical(surface_meta$sample_id, c("C01", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12"))) {
  stop("The surface subset does not match the frozen sample map.")
}
if (!identical(as.integer(table(surface_meta$site_display)), c(2L, 3L, 3L, 3L))) {
  stop("Unexpected surface group sizes.")
}

depth_meta <- metadata %>%
  filter(site == "NB") %>%
  mutate(
    depth_display = factor(depth_display, levels = c("surface", "middle", "deep")),
    pit = droplevels(pit)
  )
if (!identical(depth_meta$sample_id, c("C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18"))) {
  stop("The NA depth subset does not match the frozen sample map.")
}
pit_depth_counts <- depth_meta %>% count(pit, depth_display, name = "n")
if (nrow(pit_depth_counts) != 9L || any(pit_depth_counts$n != 1L)) {
  stop("The NA three-pit by three-depth map is incomplete; depth model not estimable as specified.")
}

design_map <- bind_rows(
  surface_meta %>%
    transmute(
      analysis = "surface_site",
      sample_id,
      group = as.character(site_display),
      pit = as.character(pit),
      depth = as.character(depth_display),
      model = "~ site",
      biological_unit = "independent metagenome"
    ),
  depth_meta %>%
    transmute(
      analysis = "NA_depth_groups",
      sample_id,
      group = as.character(depth_display),
      pit = as.character(pit),
      depth = as.character(depth_display),
      model = "~ depth",
      biological_unit = "independent metagenome"
    ),
  depth_meta %>%
    transmute(
      analysis = "NA_depth_pit_adjusted_sensitivity",
      sample_id,
      group = as.character(depth_display),
      pit = as.character(pit),
      depth = as.character(depth_display),
      model = "~ pit + depth",
      biological_unit = "one metagenome per pit-depth cell"
    )
)
write_tsv_checked(design_map, file.path(OUT, "KO_DA_design_map.tsv"))

prepared_count_path <- file.path(INTERMEDIATE, "KO_fractional_raw_counts.rds")
prepared_annotation_path <- file.path(INTERMEDIATE, "KO_annotation.rds")
prepared_offsets_path <- file.path(OUT, "KO_DA_full_catalog_TMM_offsets.tsv")
prepared_audit_path <- file.path(OUT, "KO_DA_aggregation_audit.tsv")
use_prepared_inputs <- all(file.exists(c(
  prepared_count_path, prepared_annotation_path, prepared_offsets_path, prepared_audit_path
)))

if (use_prepared_inputs) {
  message("Using staged, validated KO aggregation and full-catalog TMM inputs...")
  ko_counts <- readRDS(prepared_count_path)
  ko_annotation <- readRDS(prepared_annotation_path)
  normalisation_audit <- readr::read_tsv(prepared_offsets_path, show_col_types = FALSE)
  surface_norm <- normalisation_audit %>%
    filter(analysis == "surface_site") %>%
    arrange(match(sample_id, surface_meta$sample_id))
  depth_norm <- normalisation_audit %>%
    filter(analysis == "NA_depth_blocked_by_pit") %>%
    arrange(match(sample_id, depth_meta$sample_id))
  aggregation_audit <- readr::read_tsv(prepared_audit_path, show_col_types = FALSE)
  if (nrow(surface_norm) != 11L || nrow(depth_norm) != 9L || nrow(ko_counts) < 1000L) {
    stop("Prepared KO inputs failed dimension checks; rerun staged preparation.")
  }
} else {
message("Reading full gene counts and estimating full-catalog TMM offsets...")
counts <- data.table::fread(PATHS[["counts"]], select = c("GeneName", SAMPLES))
if (anyDuplicated(counts$GeneName) || anyNA(counts$GeneName)) {
  stop("GeneName values in Nmapped_reads.tsv must be unique and non-missing.")
}
if (any(!is.finite(as.matrix(counts[, ..SAMPLES]))) || any(as.matrix(counts[, ..SAMPLES]) < 0)) {
  stop("Raw mapped-read table contains invalid values.")
}

full_catalog_tmm <- function(sample_ids, analysis_name) {
  mat <- as.matrix(counts[, ..sample_ids])
  storage.mode(mat) <- "double"
  y <- edgeR::DGEList(counts = mat)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  tibble(
    analysis = analysis_name,
    sample_id = sample_ids,
    full_catalog_library_size = y$samples$lib.size,
    full_catalog_TMM_norm_factor = y$samples$norm.factors,
    full_catalog_effective_library_size = y$samples$lib.size * y$samples$norm.factors
  )
}

surface_norm <- full_catalog_tmm(surface_meta$sample_id, "surface_site")
depth_norm <- full_catalog_tmm(depth_meta$sample_id, "NA_depth_blocked_by_pit")
normalisation_audit <- bind_rows(surface_norm, depth_norm)
write_tsv_checked(normalisation_audit, file.path(OUT, "KO_DA_full_catalog_TMM_offsets.tsv"))

message("Reading eggNOG KO mappings and aggregating fractional KO raw counts...")
annotation <- data.table::fread(
  PATHS[["annotation"]],
  select = c("GeneName", "KO", "Name", "Description", "KEGG_Pathway", "KEGG_Module", "BRITE")
)

ko_map <- annotation[
  !is.na(KO) & KO != "" & KO != "-",
  .(KO = trimws(unlist(strsplit(gsub("ko:", "", KO, fixed = TRUE), ",", fixed = TRUE)))),
  by = GeneName
]
ko_map <- unique(ko_map[grepl("^K[0-9]{5}$", KO), .(GeneName, KO)])
ko_map[, n_valid_KO := .N, by = GeneName]
if (!nrow(ko_map) || any(ko_map$n_valid_KO < 1L)) {
  stop("No valid KO mappings were reconstructed.")
}

annotation_lookup <- annotation[, .(GeneName, Name, Description, KEGG_Pathway, KEGG_Module, BRITE)]
ko_annotation <- merge(ko_map[, .(GeneName, KO)], annotation_lookup, by = "GeneName", all.x = TRUE) %>%
  as_tibble() %>%
  group_by(KO) %>%
  summarise(
    KO_name = first_informative(Name),
    description = first_informative(Description),
    KEGG_Pathway = first_informative(KEGG_Pathway),
    KEGG_Module = first_informative(KEGG_Module),
    BRITE = first_informative(BRITE),
    n_catalog_genes = n_distinct(GeneName),
    .groups = "drop"
  )

ko_joined <- merge(
  ko_map[, .(GeneName, KO, n_valid_KO)],
  counts,
  by = "GeneName",
  allow.cartesian = TRUE
)
for (sample_id in SAMPLES) {
  data.table::set(
    ko_joined,
    j = sample_id,
    value = ko_joined[[sample_id]] / ko_joined$n_valid_KO
  )
}
ko_counts <- ko_joined[, lapply(.SD, sum, na.rm = TRUE), by = KO, .SDcols = SAMPLES]
if (nrow(ko_counts) < 1000L || any(!is.finite(as.matrix(ko_counts[, ..SAMPLES])))) {
  stop("Fractional KO count aggregation failed integrity checks.")
}

aggregation_audit <- tibble(
  metric = c(
    "metagenome_samples",
    "catalog_gene_clusters",
    "genes_with_at_least_one_valid_KO",
    "gene_to_KO_rows_after_fractional_mapping",
    "detected_fractional_KO_features"
  ),
  value = c(
    length(SAMPLES),
    nrow(counts),
    data.table::uniqueN(ko_map$GeneName),
    nrow(ko_map),
    nrow(ko_counts)
  )
)
write_tsv_checked(aggregation_audit, file.path(OUT, "KO_DA_aggregation_audit.tsv"))
}

make_group_means <- function(y, meta, group_column) {
  group <- droplevels(meta[[group_column]])
  cpm_mat <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
  summary <- tibble(KO = rownames(cpm_mat))
  for (level in levels(group)) {
    ids <- meta$sample_id[group == level]
    summary[[paste0("mean_CPM_", level)]] <- rowMeans(cpm_mat[, ids, drop = FALSE])
  }
  mean_cols <- grep("^mean_CPM_", names(summary), value = TRUE)
  mean_mat <- as.matrix(summary[, mean_cols, drop = FALSE])
  summary %>%
    mutate(
      highest_mean_group = sub("^mean_CPM_", "", mean_cols[max.col(mean_mat, ties.method = "first")]),
      lowest_mean_group = sub("^mean_CPM_", "", mean_cols[max.col(-mean_mat, ties.method = "first")]),
      max_min_log2_mean_CPM = log2((apply(mean_mat, 1, max) + 0.5) / (apply(mean_mat, 1, min) + 0.5))
    )
}

test_table <- function(qlf, test_type, contrast_label, direction_label = NA_character_) {
  tab <- as.data.frame(qlf$table)
  is_pairwise <- identical(test_type, "pairwise")
  log2fc_values <- if ("logFC" %in% names(tab) && is_pairwise) {
    tab$logFC
  } else {
    rep(NA_real_, nrow(tab))
  }
  result <- tibble(KO = rownames(tab)) %>%
    mutate(
      test_type = test_type,
      contrast = contrast_label,
      log2FC = log2fc_values,
      logCPM = if ("logCPM" %in% names(tab)) tab$logCPM else NA_real_,
      F = tab$F,
      PValue = tab$PValue,
      FDR = p.adjust(tab$PValue, method = "BH"),
      direction = case_when(
        !is_pairwise ~ NA_character_,
        log2FC > 0 ~ paste0("higher in ", direction_label),
        log2FC < 0 ~ paste0("lower in ", direction_label),
        TRUE ~ "no estimated difference"
      )
    )
  result
}

fit_da <- function(analysis_name, meta, formula, group_column, norm_table, contrasts, omnibus_coef) {
  sample_ids <- meta$sample_id
  mat <- as.matrix(ko_counts[, ..sample_ids])
  rownames(mat) <- ko_counts$KO
  storage.mode(mat) <- "double"

  design <- model.matrix(formula, data = meta)
  if (qr(design)$rank != ncol(design)) {
    stop("Design matrix is not full rank for ", analysis_name)
  }
  if (nrow(design) - qr(design)$rank < 2L) {
    stop("Too few residual degrees of freedom for ", analysis_name)
  }

  norm_table <- norm_table %>% arrange(match(sample_id, sample_ids))
  if (!identical(norm_table$sample_id, sample_ids)) {
    stop("TMM offset order mismatch for ", analysis_name)
  }
  y <- edgeR::DGEList(
    counts = mat,
    lib.size = norm_table$full_catalog_library_size
  )
  y$samples$norm.factors <- norm_table$full_catalog_TMM_norm_factor
  keep <- edgeR::filterByExpr(y, design = design, min.count = 10)
  if (sum(keep) < 100L) {
    stop("Unexpectedly few KOs retained after filtering for ", analysis_name)
  }
  y <- y[keep, , keep.lib.sizes = TRUE]
  expected_offset <- log(norm_table$full_catalog_effective_library_size)
  observed_offset <- as.numeric(edgeR::getOffset(y))
  if (!isTRUE(all.equal(observed_offset, expected_offset, tolerance = 1e-10))) {
    stop("Full-catalog TMM offsets were not preserved for ", analysis_name)
  }

  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  omnibus <- edgeR::glmQLFTest(fit, coef = omnibus_coef)
  group_means <- make_group_means(y, meta, group_column)

  test_results <- list(
    test_table(
      omnibus,
      "omnibus",
      case_when(
        analysis_name == "surface_site" ~ "site (4 levels)",
        analysis_name == "NA_depth_pit_adjusted_sensitivity" ~ "depth (3 levels; pit-adjusted)",
        TRUE ~ "depth (3 levels)"
      )
    )
  )
  for (contrast_name in names(contrasts)) {
    contrast_spec <- contrasts[[contrast_name]]
    qlf <- edgeR::glmQLFTest(fit, contrast = contrast_spec$vector)
    test_results[[contrast_name]] <- test_table(
      qlf,
      "pairwise",
      contrast_name,
      contrast_spec$positive_group
    )
  }
  results <- bind_rows(test_results) %>%
    left_join(group_means, by = "KO") %>%
    left_join(ko_annotation, by = "KO") %>%
    mutate(
      analysis = analysis_name,
      n_samples = nrow(meta),
      n_tested_KOs = nrow(y),
      n_groups = nlevels(droplevels(meta[[group_column]])),
      model = deparse(formula),
      normalization = "TMM estimated from the full gene catalog within analysis subset",
      feature_definition = "raw mapped reads fractionally assigned across valid KO labels per gene"
    ) %>%
    select(
      analysis, test_type, contrast, KO, KO_name, description, KEGG_Pathway,
      KEGG_Module, BRITE, n_catalog_genes, n_samples, n_groups, n_tested_KOs,
      model, normalization, feature_definition, log2FC, logCPM, F, PValue, FDR,
      direction, highest_mean_group, lowest_mean_group, max_min_log2_mean_CPM,
      starts_with("mean_CPM_")
    ) %>%
    arrange(test_type, contrast, FDR, PValue, desc(abs(log2FC)))

  model_audit <- tibble(
    analysis = analysis_name,
    formula = deparse(formula),
    coefficient = colnames(design),
    n_samples = nrow(meta),
    residual_df = nrow(design) - qr(design)$rank,
    retained_KOs = nrow(y),
    filtered_out_KOs = nrow(mat) - nrow(y),
    TMM_effective_library_size_min = min(norm_table$full_catalog_effective_library_size),
    TMM_effective_library_size_max = max(norm_table$full_catalog_effective_library_size)
  )

  list(results = results, model_audit = model_audit)
}

message("Fitting surface-site KO differential-abundance model...")
surface_meta <- surface_meta %>%
  mutate(site_display = relevel(site_display, ref = "AZ"))
surface_design <- model.matrix(~ site_display, data = surface_meta)
surface_contrasts <- list(
  SG_vs_AZ = list(vector = c(0, 1, 0, 0), positive_group = "SG"),
  LC_vs_AZ = list(vector = c(0, 0, 1, 0), positive_group = "LC"),
  NA_surface_vs_AZ = list(vector = c(0, 0, 0, 1), positive_group = "NA surface"),
  LC_vs_SG = list(vector = c(0, -1, 1, 0), positive_group = "LC"),
  NA_surface_vs_SG = list(vector = c(0, -1, 0, 1), positive_group = "NA surface"),
  NA_surface_vs_LC = list(vector = c(0, 0, -1, 1), positive_group = "NA surface")
)
if (!identical(colnames(surface_design), c("(Intercept)", "site_displaySG", "site_displayLC", "site_displayNA"))) {
  stop("Unexpected surface-site design encoding.")
}
surface_fit <- fit_da(
  "surface_site",
  surface_meta,
  ~ site_display,
  "site_display",
  surface_norm,
  surface_contrasts,
  omnibus_coef = c(2, 3, 4)
)

message("Fitting NA-depth KO differential-abundance model (surface, middle, deep)...")
depth_meta <- depth_meta %>%
  mutate(
    pit = relevel(pit, ref = "NB1"),
    depth_display = relevel(depth_display, ref = "surface")
  )
depth_group_design <- model.matrix(~ depth_display, data = depth_meta)
if (!identical(colnames(depth_group_design), c("(Intercept)", "depth_displaymiddle", "depth_displaydeep"))) {
  stop("Unexpected NA depth-group model encoding.")
}
depth_group_contrasts <- list(
  middle_vs_surface = list(vector = c(0, 1, 0), positive_group = "middle"),
  deep_vs_surface = list(vector = c(0, 0, 1), positive_group = "deep"),
  deep_vs_middle = list(vector = c(0, -1, 1), positive_group = "deep")
)
depth_group_fit <- fit_da(
  "NA_depth_groups",
  depth_meta,
  ~ depth_display,
  "depth_display",
  depth_norm,
  depth_group_contrasts,
  omnibus_coef = c(2, 3)
)

message("Fitting pit-adjusted NA-depth sensitivity model...")
depth_pit_design <- model.matrix(~ pit + depth_display, data = depth_meta)
if (!identical(colnames(depth_pit_design), c("(Intercept)", "pitNB2", "pitNB3", "depth_displaymiddle", "depth_displaydeep"))) {
  stop("Unexpected pit-adjusted NA depth model encoding.")
}
depth_pit_contrasts <- list(
  middle_vs_surface = list(vector = c(0, 0, 0, 1, 0), positive_group = "middle"),
  deep_vs_surface = list(vector = c(0, 0, 0, 0, 1), positive_group = "deep"),
  deep_vs_middle = list(vector = c(0, 0, 0, -1, 1), positive_group = "deep")
)
depth_pit_fit <- fit_da(
  "NA_depth_pit_adjusted_sensitivity",
  depth_meta,
  ~ pit + depth_display,
  "depth_display",
  depth_norm,
  depth_pit_contrasts,
  omnibus_coef = c(4, 5)
)

all_results <- bind_rows(surface_fit$results, depth_group_fit$results, depth_pit_fit$results)
model_audit <- bind_rows(surface_fit$model_audit, depth_group_fit$model_audit, depth_pit_fit$model_audit)
write_tsv_checked(all_results, file.path(OUT, "KO_DA_all_tests.tsv"))
write_tsv_checked(model_audit, file.path(OUT, "KO_DA_model_audit.tsv"))

summary_table <- all_results %>%
  group_by(analysis, test_type, contrast, n_samples, n_groups, n_tested_KOs, model) %>%
  summarise(
    KOs_FDR_lt_0_05 = sum(FDR < 0.05, na.rm = TRUE),
    KOs_FDR_lt_0_10 = sum(FDR < 0.10, na.rm = TRUE),
    KOs_raw_p_lt_0_05 = sum(PValue < 0.05, na.rm = TRUE),
    min_raw_p = min(PValue, na.rm = TRUE),
    min_FDR = min(FDR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(analysis, test_type, contrast)
write_tsv_checked(summary_table, file.path(OUT, "KO_DA_test_summary.tsv"))

top_features <- all_results %>%
  group_by(analysis, test_type, contrast) %>%
  arrange(FDR, PValue, desc(abs(log2FC)), .by_group = TRUE) %>%
  slice_head(n = 50L) %>%
  ungroup() %>%
  mutate(rank_within_test = ave(PValue, interaction(analysis, test_type, contrast), FUN = rank, ties.method = "first"))
write_tsv_checked(top_features, file.path(OUT, "KO_DA_top50_per_test.tsv"))

# Inspect every FDR-supported unadjusted depth-group hit at the individual
# metagenome level. This diagnostic makes it possible to see whether a group
# result is repeated across pits or is dominated by a particular pit-depth cell.
depth_group_hits <- all_results %>%
  filter(analysis == "NA_depth_groups", test_type == "omnibus", FDR < 0.05) %>%
  arrange(FDR, PValue)

depth_pit_diagnostic <- tibble()
if (nrow(depth_group_hits) > 0L) {
  hit_ids <- depth_group_hits$KO
  depth_count_matrix <- as.matrix(
    as.data.frame(ko_counts)[match(hit_ids, ko_counts$KO), depth_meta$sample_id, drop = FALSE]
  )
  rownames(depth_count_matrix) <- hit_ids
  depth_log_cpm <- edgeR::cpm(
    depth_count_matrix,
    lib.size = depth_norm$full_catalog_effective_library_size,
    log = TRUE,
    prior.count = 2
  )
  depth_pit_diagnostic <- as_tibble(depth_log_cpm, rownames = "KO") %>%
    pivot_longer(-KO, names_to = "sample_id", values_to = "log2_CPM") %>%
    left_join(depth_meta %>% select(sample_id, pit, depth_display), by = "sample_id") %>%
    left_join(depth_group_hits %>% select(KO, KO_name, description, FDR), by = "KO") %>%
    mutate(
      pit = factor(pit, levels = c("NB1", "NB2", "NB3")),
      depth_display = factor(depth_display, levels = c("surface", "middle", "deep")),
      feature_label = safe_label(KO, KO_name)
    )
  write_tsv_checked(depth_pit_diagnostic, file.path(OUT, "KO_DA_NA_depth_hit_pit_diagnostic.tsv"))
}

plot_one_analysis <- function(analysis_name, title) {
  global <- all_results %>%
    filter(analysis == analysis_name, test_type == "omnibus") %>%
    arrange(FDR, PValue)
  pairwise_summary <- summary_table %>%
    filter(analysis == analysis_name, test_type == "pairwise") %>%
    mutate(contrast = factor(contrast, levels = contrast))
  plot_features <- global %>% slice_head(n = 24L) %>%
    mutate(feature_label = safe_label(KO, KO_name))
  mean_cols <- grep("^mean_CPM_", names(plot_features), value = TRUE)
  mean_cols <- mean_cols[vapply(
    plot_features[, mean_cols, drop = FALSE],
    function(x) any(!is.na(x)),
    logical(1)
  )]
  heatmap_df <- plot_features %>%
    select(KO, feature_label, all_of(mean_cols)) %>%
    pivot_longer(all_of(mean_cols), names_to = "group", values_to = "mean_CPM") %>%
    mutate(group = sub("^mean_CPM_", "", group)) %>%
    group_by(KO) %>%
    mutate(
      log_mean_CPM = log2(mean_CPM + 0.5),
      z = as.numeric(scale(log_mean_CPM)),
      z = if_else(is.finite(z), z, 0),
      z = pmax(pmin(z, 2.5), -2.5)
    ) %>%
    ungroup() %>%
    mutate(feature_label = factor(feature_label, levels = rev(unique(plot_features$feature_label))))

  p_bar <- ggplot(pairwise_summary, aes(x = contrast, y = KOs_FDR_lt_0_05)) +
    geom_col(width = 0.68, fill = "#3D6F97") +
    geom_text(aes(label = KOs_FDR_lt_0_05), vjust = -0.35, size = 3.6) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = "Pairwise FDR-supported KOs",
      x = NULL,
      y = "Number of KOs (FDR < 0.05)"
    ) +
    theme_classic(base_size = 11, base_family = "Helvetica") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(face = "bold"))

  global_summary <- summary_table %>%
    filter(analysis == analysis_name, test_type == "omnibus")
  support_note <- sprintf(
    "Omnibus: %d/%d KOs at FDR < 0.05; heatmap shows 24 lowest-FDR KOs (row z-score of group mean CPM).",
    global_summary$KOs_FDR_lt_0_05,
    global_summary$n_tested_KOs
  )
  p_heat <- ggplot(heatmap_df, aes(x = group, y = feature_label, fill = z)) +
    geom_tile(color = "grey84", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-2.5, 2.5), name = "Row z-score\nof mean CPM"
    ) +
    labs(title = "Most differential KOs", subtitle = support_note, x = NULL, y = NULL) +
    theme_minimal(base_size = 10.5, base_family = "Helvetica") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(size = 7.8),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 8.4),
      legend.position = "right"
    )
  (p_heat / p_bar) + plot_layout(heights = c(3.1, 1)) + plot_annotation(
    title = title,
    subtitle = "Exploratory KO-level differential abundance from iDNA metagenomic mapped reads"
  )
}

pdf_path <- file.path(OUT, "KO_DA_exploratory_overview.pdf")
grDevices::cairo_pdf(pdf_path, width = 11.7, height = 8.3, onefile = TRUE)
print(plot_one_analysis("surface_site", "Surface sites: AZ, SG, LC, and NA surface"))
print(plot_one_analysis("NA_depth_groups", "Nahuelbuta depth groups: surface, middle, and deep"))
grDevices::dev.off()
if (!file.exists(pdf_path) || file.info(pdf_path)$size <= 0) {
  stop("Exploratory KO overview PDF was not created.")
}

if (nrow(depth_pit_diagnostic) > 0L) {
  sensitivity_min_fdr <- summary_table %>%
    filter(analysis == "NA_depth_pit_adjusted_sensitivity", test_type == "omnibus") %>%
    pull(min_FDR)
  diagnostic_pdf <- file.path(OUT, "KO_DA_NA_depth_pit_diagnostic.pdf")
  p_depth_diagnostic <- ggplot(
    depth_pit_diagnostic,
    aes(x = depth_display, y = log2_CPM, colour = pit)
  ) +
    geom_point(size = 2.2, alpha = 0.9, position = position_jitter(width = 0.06, height = 0)) +
    facet_wrap(~ feature_label, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = c(NB1 = "#1B9E77", NB2 = "#7570B3", NB3 = "#D95F02")) +
    labs(
      title = "Individual-pit diagnostic for unadjusted NA depth-group KO hits",
      subtitle = sprintf(
        "%d omnibus FDR-supported KOs in ~ depth; the pit-adjusted omnibus sensitivity model has minimum FDR = %.3f.",
        n_distinct(depth_pit_diagnostic$KO), sensitivity_min_fdr
      ),
      x = NULL,
      y = "log2 CPM (full-catalog TMM offset)",
      colour = "Pit"
    ) +
    theme_classic(base_size = 10.5, base_family = "Helvetica") +
    theme(
      strip.text = element_text(size = 8.5),
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
  grDevices::cairo_pdf(diagnostic_pdf, width = 11.7, height = 8.3)
  print(p_depth_diagnostic)
  grDevices::dev.off()
  if (!file.exists(diagnostic_pdf) || file.info(diagnostic_pdf)$size <= 0) {
    stop("NA depth pit diagnostic PDF was not created.")
  }
}

summary_lines <- c(
  "# Exploratory KO-level differential-abundance analysis",
  "",
  "## Scope",
  "",
  "This analysis uses iDNA metagenomic raw mapped reads. It tests differential KO abundance (encoded functional potential), not differential gene expression or transcription.",
  "",
  "- Surface-site model: `~ site`, with AZ (n = 2), SG (n = 3), LC (n = 3), and NA surface (n = 3).",
  "- Nahuelbuta primary depth model: `~ depth`, comparing surface, middle, and deep groups (n = 3 per depth). A separate `~ pit + depth` model is retained only as a sensitivity analysis; it does not redefine pit as the grouping variable.",
  "- KO definition: genes with multiple valid KO labels were divided fractionally across those labels before aggregation.",
  "- Normalisation: TMM factors were estimated from the full gene catalog within each tested subset, then retained as offsets for KO-level edgeR QL models.",
  "- FDR family: Benjamini-Hochberg correction was applied across all retained KOs separately for each omnibus or pairwise test.",
  "",
  "## Test summary",
  "",
  "| Analysis | Test | Contrast | Retained KOs | Raw p < 0.05 | FDR < 0.05 | Minimum raw p | Minimum FDR |",
  "|---|---|---|---:|---:|---:|---:|---:|"
)
summary_lines <- c(
  summary_lines,
  apply(summary_table, 1, function(row) {
    sprintf(
      "| %s | %s | %s | %s | %s | %s | %.3g | %.3g |",
      row[["analysis"]], row[["test_type"]], row[["contrast"]],
      row[["n_tested_KOs"]], row[["KOs_raw_p_lt_0_05"]], row[["KOs_FDR_lt_0_05"]],
      as.numeric(row[["min_raw_p"]]), as.numeric(row[["min_FDR"]])
    )
  }),
  "",
  "## Interpretation boundary",
  "",
  "These outputs are exploratory. They can identify KO candidates and repeated directional patterns for biological interpretation, but they do not establish enzyme activity, process rates, or causal environmental drivers. Any candidate selected for manuscript use should be checked against its raw sample distribution, functional annotation, marker specificity, and the current multivariate evidence.",
  "",
  "## Outputs",
  "",
  "- `KO_DA_all_tests.tsv`: complete KO-by-test table with effect sizes, raw p values, FDR, annotation, and group mean CPM.",
  "- `KO_DA_top50_per_test.tsv`: top 50 KOs per test for quick inspection.",
  "- `KO_DA_NA_depth_hit_pit_diagnostic.tsv` and `KO_DA_NA_depth_pit_diagnostic.pdf`: individual-pit display of unadjusted NA depth-group FDR hits; this is a diagnostic, not a second hypothesis test.",
  "- `KO_DA_test_summary.tsv`: counts of raw-p and FDR-supported KOs by test.",
  "- `KO_DA_exploratory_overview.pdf`: two-page overview; top rows are selected by omnibus FDR/raw p and are a display, not an additional test.",
  "- `KO_DA_design_map.tsv`, `KO_DA_model_audit.tsv`, and `KO_DA_full_catalog_TMM_offsets.tsv`: reproducibility audit records."
)
writeLines(summary_lines, file.path(OUT, "KO_DA_exploratory_summary.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))
message("Complete. Outputs written to: ", OUT)
