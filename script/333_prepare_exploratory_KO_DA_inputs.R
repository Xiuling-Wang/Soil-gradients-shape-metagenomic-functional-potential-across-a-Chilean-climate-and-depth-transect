#!/usr/bin/env Rscript
# Staged input preparation for script/332_exploratory_KO_differential_abundance.R.
#
# Use one phase at a time in constrained execution environments:
#   Rscript script/333_prepare_exploratory_KO_DA_inputs.R annotation
#   Rscript script/333_prepare_exploratory_KO_DA_inputs.R surface_tmm
#   Rscript script/333_prepare_exploratory_KO_DA_inputs.R depth_tmm
#   Rscript script/333_prepare_exploratory_KO_DA_inputs.R aggregate
#
# `all` is provided for ordinary local use, but each phase remains independently
# reproducible and writes only inside the exploratory output bundle.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(edgeR)
})

options(scipen = 999)
set.seed(42)
invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "UTF-8"), silent = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[[1]] else "all"
valid_phases <- c("annotation", "surface_tmm", "depth_tmm", "aggregate", "all")
if (!phase %in% valid_phases) {
  stop("Phase must be one of: ", paste(valid_phases, collapse = ", "))
}

script_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", script_args, value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript.")
SCRIPT <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
ROOT <- normalizePath(file.path(dirname(SCRIPT), ".."), mustWork = TRUE)
OUT <- file.path(ROOT, "analysis", "exploratory_KO_differential_abundance_2026-08-29")
INTERMEDIATE <- file.path(OUT, "intermediate")
dir.create(INTERMEDIATE, recursive = TRUE, showWarnings = FALSE)

SAMPLES <- c("C01", sprintf("C%02d", 3:18))
PATHS <- c(
  counts = file.path(ROOT, "data/processed/metagenome_output/Genecatalog/Nmapped_reads.tsv"),
  annotation = file.path(ROOT, "data/processed/metagenome_output/Genecatalog/eggNog.tsv"),
  metadata = file.path(ROOT, "data/metadata/env_2024_meta.tsv")
)
missing_inputs <- PATHS[!file.exists(PATHS)]
if (length(missing_inputs)) stop("Missing input(s):\n", paste(missing_inputs, collapse = "\n"))

first_informative <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "-"]
  if (length(x)) x[[1]] else NA_character_
}

write_tsv_checked <- function(x, path) {
  readr::write_tsv(x, path, na = "NA")
  if (!file.exists(path) || file.info(path)$size <= 0) stop("Failed to write: ", path)
}

metadata <- readr::read_tsv(PATHS[["metadata"]], show_col_types = FALSE, na = "") %>%
  filter(sample_id %in% SAMPLES) %>%
  arrange(match(sample_id, SAMPLES)) %>%
  mutate(
    site = factor(site, levels = c("AZ", "SG", "LC", "NB")),
    depth_label = factor(depth_label, levels = c("5_10cm", "20_40cm", "60_80cm")),
    site_display = recode(as.character(site), NB = "NA"),
    depth_display = recode(as.character(depth_label), `5_10cm` = "surface", `20_40cm` = "middle", `60_80cm` = "deep"),
    pit = factor(site_pit)
  )
if (!identical(metadata$sample_id, SAMPLES) || !all(metadata$dna_type == "iDNA")) {
  stop("The current metadata does not match the 17 iDNA metagenome universe.")
}

surface_ids <- metadata %>% filter(depth_label == "5_10cm") %>% pull(sample_id)
depth_ids <- metadata %>% filter(site == "NB") %>% pull(sample_id)
depth_check <- metadata %>%
  filter(site == "NB") %>%
  count(site_pit, depth_label, name = "n")
if (length(surface_ids) != 11L || length(depth_ids) != 9L || nrow(depth_check) != 9L || any(depth_check$n != 1L)) {
  stop("Surface or NA-depth sample map failed validation.")
}

map_path <- file.path(INTERMEDIATE, "KO_fractional_map.rds")
streamed_map_path <- file.path(INTERMEDIATE, "KO_fractional_map.tsv")
annotation_path <- file.path(INTERMEDIATE, "KO_annotation.rds")
ko_count_path <- file.path(INTERMEDIATE, "KO_fractional_raw_counts.rds")
tmm_path <- file.path(OUT, "KO_DA_full_catalog_TMM_offsets.tsv")

prepare_annotation <- function() {
  if (file.exists(streamed_map_path)) {
    message("Phase annotation: using streamed GeneName-KO map...")
    ko_map <- data.table::fread(
      streamed_map_path,
      header = FALSE,
      col.names = c("GeneName", "KO", "n_valid_KO")
    )
    ko_map <- unique(ko_map[grepl("^K[0-9]{5}$", KO)])
    if (!nrow(ko_map) || any(ko_map$n_valid_KO < 1L)) {
      stop("Streamed KO map failed integrity checks.")
    }
    ko_annotation <- ko_map[, .(n_catalog_genes = uniqueN(GeneName)), by = KO] %>%
      as_tibble() %>%
      mutate(
        KO_name = NA_character_, description = NA_character_,
        KEGG_Pathway = NA_character_, KEGG_Module = NA_character_, BRITE = NA_character_
      ) %>%
      select(KO, KO_name, description, KEGG_Pathway, KEGG_Module, BRITE, n_catalog_genes)
  } else {
    message("Phase annotation: reading eggNOG annotation table...")
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
    if (!nrow(ko_map)) stop("No valid KO mappings reconstructed.")
    lookup <- annotation[, .(GeneName, Name, Description, KEGG_Pathway, KEGG_Module, BRITE)]
    ko_annotation <- merge(ko_map[, .(GeneName, KO)], lookup, by = "GeneName", all.x = TRUE) %>%
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
  }
  saveRDS(ko_map, map_path, compress = FALSE)
  saveRDS(ko_annotation, annotation_path, compress = FALSE)
  if (!file.exists(map_path) || !file.exists(annotation_path)) stop("KO mapping cache was not written.")
  message("Annotation phase complete: ", nrow(ko_annotation), " KOs.")
}

prepare_tmm <- function(sample_ids, analysis_name) {
  message("Phase ", analysis_name, ": reading full gene count catalog for TMM...")
  raw_counts <- data.table::fread(PATHS[["counts"]], select = c("GeneName", sample_ids))
  if (anyDuplicated(raw_counts$GeneName)) stop("Gene names are not unique in raw count catalog.")
  mat <- as.matrix(raw_counts[, ..sample_ids])
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

write_tmm <- function() {
  existing <- if (file.exists(tmm_path)) readr::read_tsv(tmm_path, show_col_types = FALSE) else tibble()
  if (phase %in% c("surface_tmm", "all")) {
    surface_tmm <- prepare_tmm(surface_ids, "surface_site")
    existing <- existing %>% filter(analysis != "surface_site") %>% bind_rows(surface_tmm)
  }
  if (phase %in% c("depth_tmm", "all")) {
    depth_tmm <- prepare_tmm(depth_ids, "NA_depth_blocked_by_pit")
    existing <- existing %>% filter(analysis != "NA_depth_blocked_by_pit") %>% bind_rows(depth_tmm)
  }
  write_tsv_checked(existing, tmm_path)
}

prepare_aggregation <- function() {
  if (!file.exists(map_path) || !file.exists(annotation_path)) {
    stop("Run the annotation phase before aggregate.")
  }
  message("Phase aggregate: reading raw mapped counts and fractional KO map...")
  ko_map <- readRDS(map_path)
  raw_counts <- data.table::fread(PATHS[["counts"]], select = c("GeneName", SAMPLES))
  if (anyDuplicated(raw_counts$GeneName)) stop("Gene names are not unique in raw count catalog.")
  joined <- merge(
    ko_map[, .(GeneName, KO, n_valid_KO)], raw_counts,
    by = "GeneName", allow.cartesian = TRUE
  )
  for (sample_id in SAMPLES) {
    data.table::set(joined, j = sample_id, value = joined[[sample_id]] / joined$n_valid_KO)
  }
  ko_counts <- joined[, lapply(.SD, sum, na.rm = TRUE), by = KO, .SDcols = SAMPLES]
  if (nrow(ko_counts) < 1000L || any(!is.finite(as.matrix(ko_counts[, ..SAMPLES])))) {
    stop("KO aggregation failed integrity checks.")
  }
  saveRDS(ko_counts, ko_count_path, compress = FALSE)
  aggregation_audit <- tibble(
    metric = c("metagenome_samples", "catalog_gene_clusters", "genes_with_at_least_one_valid_KO", "gene_to_KO_rows_after_fractional_mapping", "detected_fractional_KO_features"),
    value = c(length(SAMPLES), nrow(raw_counts), uniqueN(ko_map$GeneName), nrow(ko_map), nrow(ko_counts))
  )
  write_tsv_checked(aggregation_audit, file.path(OUT, "KO_DA_aggregation_audit.tsv"))
  message("Aggregation phase complete: ", nrow(ko_counts), " detected fractional KOs.")
}

if (phase %in% c("annotation", "all")) prepare_annotation()
if (phase %in% c("surface_tmm", "depth_tmm", "all")) write_tmm()
if (phase %in% c("aggregate", "all")) prepare_aggregation()
