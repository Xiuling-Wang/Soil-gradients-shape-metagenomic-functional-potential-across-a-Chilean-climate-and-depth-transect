#!/usr/bin/env Rscript
# V43 targeted robustness analyses.
#
# This script leaves V42 untouched and writes only to the V43 analysis folder.
# It evaluates (i) robust-CLR/Aitchison sensitivity for KO and MAG composition
# and (ii) MAG marker-carrier prevalence across completeness strata.

suppressPackageStartupMessages({
  library(data.table)
  library(permute)
  library(tidyverse)
  library(vegan)
})

options(scipen = 999)
set.seed(42)
invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "UTF-8"), silent = TRUE)))

ROOT <- normalizePath(Sys.getenv("CHILE_PROJECT_ROOT", getwd()), mustWork = TRUE)
setwd(ROOT)
OUT_ROOT <- file.path(ROOT, "analysis/V43_targeted_robustness_revision_2026-08-27")
ROBUST_DIR <- file.path(OUT_ROOT, "robustness")
SOURCE_DIR <- file.path(OUT_ROOT, "source_data")
dir.create(ROBUST_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

V40_ROOT <- Sys.glob(file.path(ROOT, "*", "V40_targeted_revision_2026-08-26"))
if (length(V40_ROOT) != 1L) {
  stop("Could not resolve the unique V40 reference package")
}

SAMPLES <- c("C01", sprintf("C%02d", 3:18))
DEPTH_ORDER <- c("5_10cm", "20_40cm", "60_80cm")

PATHS <- c(
  counts = "data/processed/metagenome_output/Genecatalog/Nmapped_reads.tsv",
  lengths = "data/processed/metagenome_output/Genecatalog/seq_len_nucl.tsv",
  annotation = "data/processed/metagenome_output/Genecatalog/eggNog.tsv",
  mag_counts = "data/processed/metagenome_output/MAGs/counts_genomes.tsv",
  mag_quality = "data/processed/metagenome_output/MAGs/genome_quality.tsv",
  mag_taxonomy = "data/processed/metagenome_output/MAGs/gtdb_taxonomy.tsv",
  dram_annotations = "data/processed/metagenome_output/MAGs/dram/annotations.tsv",
  metadata = "data/metadata/env_2024_meta.tsv",
  marker_assignments = paste0(
    "analysis/V18_figure_caption_cleanup_2026-08-19/source_data/",
    "Figure5_MAG_marker_assignments.tsv"
  ),
  marker_registry = paste0(
    "analysis/V18_figure_caption_cleanup_2026-08-19/source_data/",
    "Figure5_frozen_primary_marker_registry.tsv"
  ),
  bray_location = file.path(
    V40_ROOT, "phase1_support/statistics/exact_direct_location_permanova.tsv"
  ),
  bray_dispersion = file.path(
    V40_ROOT, "phase1_support/statistics/exact_direct_location_dispersion.tsv"
  )
)

missing_inputs <- PATHS[!file.exists(PATHS)]
if (length(missing_inputs)) {
  stop("Missing V43 input(s): ", paste(missing_inputs, collapse = ", "))
}

write_tsv_checked <- function(x, path) {
  readr::write_tsv(x, path, na = "NA")
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Failed to write output: ", path)
  }
  invisible(path)
}

classification <- function(p, q = p) {
  case_when(
    is.na(p) ~ "effect size only",
    q < 0.05 ~ "FDR-supported",
    p < 0.05 ~ "nominal only; not FDR-supported",
    TRUE ~ "not nominally supported"
  )
}

# ---- Frozen design ---------------------------------------------------------
metadata <- read_tsv(PATHS[["metadata"]], show_col_types = FALSE, na = "") %>%
  filter(sample_id %in% SAMPLES) %>%
  arrange(match(sample_id, SAMPLES)) %>%
  mutate(
    site = factor(site, levels = c("AZ", "SG", "LC", "NB")),
    depth_label = factor(depth_label, levels = DEPTH_ORDER),
    site_pit = factor(site_pit)
  )

stopifnot(
  identical(metadata$sample_id, SAMPLES),
  all(metadata$dna_type == "iDNA")
)

depth_metadata <- metadata %>%
  filter(site == "NB") %>%
  mutate(site_pit = droplevels(site_pit), depth_label = droplevels(depth_label))
stopifnot(
  identical(depth_metadata$sample_id, c("C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18")),
  all(table(depth_metadata$site_pit, depth_metadata$depth_label) == 1L)
)

permutations_of <- function(values) {
  if (length(values) == 1L) return(matrix(values, nrow = 1L))
  do.call(rbind, lapply(seq_along(values), function(index) {
    cbind(values[index], permutations_of(values[-index]))
  }))
}

enumerate_block_permutations <- function(blocks) {
  blocks <- as.character(blocks)
  block_levels <- unique(blocks)
  positions <- lapply(block_levels, function(block) which(blocks == block))
  local_permutations <- lapply(positions, permutations_of)
  combinations <- expand.grid(
    lapply(local_permutations, function(x) seq_len(nrow(x))),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  full <- matrix(NA_integer_, nrow = nrow(combinations), ncol = length(blocks))
  for (row_index in seq_len(nrow(combinations))) {
    permutation <- seq_along(blocks)
    for (block_index in seq_along(block_levels)) {
      destination <- positions[[block_index]]
      source_row <- combinations[row_index, block_index]
      permutation[destination] <- local_permutations[[block_index]][source_row, ]
    }
    full[row_index, ] <- permutation
  }
  full
}

depth_full <- enumerate_block_permutations(depth_metadata$site_pit)
depth_identity <- apply(depth_full, 1, function(x) all(unname(x) == seq_len(nrow(depth_metadata))))
depth_nonidentity <- depth_full[!depth_identity, , drop = FALSE]
stopifnot(nrow(depth_full) == 216L, nrow(depth_nonidentity) == 215L)

# ---- Reconstruct the fractional-KO and MAG abundance matrices -------------
message("Reconstructing the fractional KO profile from the frozen gene catalog...")
counts <- fread(PATHS[["counts"]], select = c("GeneName", SAMPLES))
lengths <- fread(PATHS[["lengths"]], header = FALSE, col.names = c("GeneName", "Length"))
annotation <- fread(PATHS[["annotation"]], select = c("GeneName", "KO"))

gene_rpk <- merge(counts, lengths, by = "GeneName", all.x = TRUE)
gene_rpk[, total_mapped_reads := rowSums(.SD, na.rm = TRUE), .SDcols = SAMPLES]
gene_rpk <- gene_rpk[
  total_mapped_reads > 0 & is.finite(Length) & Length > 0
]
for (sample_id in SAMPLES) {
  set(gene_rpk, j = sample_id, value = gene_rpk[[sample_id]] / (gene_rpk$Length / 1000))
}

ko_map <- annotation[
  !is.na(KO) & KO != "" & KO != "-",
  .(KO = trimws(unlist(strsplit(gsub("ko:", "", KO, fixed = TRUE), ",", fixed = TRUE)))),
  by = GeneName
]
ko_map <- unique(ko_map[grepl("^K[0-9]{5}$", KO)])
ko_map[, n_valid_KO := .N, by = GeneName]
ko_joined <- merge(
  ko_map[, .(GeneName, KO, n_valid_KO)],
  gene_rpk[, c("GeneName", SAMPLES), with = FALSE],
  by = "GeneName",
  allow.cartesian = TRUE
)
for (sample_id in SAMPLES) {
  set(ko_joined, j = sample_id, value = ko_joined[[sample_id]] / ko_joined$n_valid_KO)
}
ko_rpk <- ko_joined[, lapply(.SD, sum, na.rm = TRUE), by = KO, .SDcols = SAMPLES]
stopifnot(nrow(ko_rpk) == 9031L)

rm(counts, lengths, annotation, gene_rpk, ko_map, ko_joined)
invisible(gc())

message("Building the genome-size-normalised MAG profile...")
mag_counts <- fread(PATHS[["mag_counts"]])
setnames(mag_counts, "Sample", "MAG")
mag_quality_size <- fread(PATHS[["mag_quality"]], select = c("Bin Id", "Genome_Size"))
setnames(mag_quality_size, "Bin Id", "MAG")
mag_rpk <- merge(mag_counts, mag_quality_size, by = "MAG", all.x = TRUE)
stopifnot(nrow(mag_rpk) == 228L, !anyNA(mag_rpk$Genome_Size), all(mag_rpk$Genome_Size > 0))
for (sample_id in SAMPLES) {
  set(mag_rpk, j = sample_id, value = mag_rpk[[sample_id]] / (mag_rpk$Genome_Size / 1000))
}

profiles <- list(
  "Fractional KO catalog" = ko_rpk[, c("KO", SAMPLES), with = FALSE],
  "Genome-size-normalised MAG assemblage" = mag_rpk[, c("MAG", SAMPLES), with = FALSE]
)

# ---- Robust-CLR / Aitchison sensitivity -----------------------------------
# vegan's rCLR replaces zero entries by low-rank matrix completion after
# centring log-positive abundances. Euclidean distance in this transformed
# space is used as the Aitchison sensitivity distance.
rclr_distance <- function(feature_table, feature_id, sample_ids) {
  x <- t(as.matrix(feature_table[, ..sample_ids]))
  storage.mode(x) <- "double"
  colnames(x) <- feature_table[[feature_id]]
  rownames(x) <- sample_ids
  x <- x[, colSums(x, na.rm = TRUE) > 0, drop = FALSE]
  if (ncol(x) < 3L) stop("Too few non-zero features for rCLR")
  rank <- as.integer(min(3L, nrow(x) - 1L, ncol(x) - 1L))
  transformed <- decostand(x, method = "rclr", impute = TRUE, ropt = rank)
  if (!all(is.finite(transformed))) stop("rCLR returned non-finite values")
  list(
    distance = dist(transformed, method = "euclidean"),
    n_features = ncol(x),
    rclr_rank = rank,
    zero_fraction = mean(x == 0)
  )
}

run_rclr_test <- function(profile_name, feature_table, subset_name) {
  feature_id <- names(feature_table)[1]
  if (subset_name == "surface_site") {
    subset_meta <- metadata %>% filter(depth_label == "5_10cm")
    group <- droplevels(subset_meta$site)
    permutations <- 999L
    scheme <- "999 unrestricted permutations among 11 surface samples"
  } else {
    subset_meta <- depth_metadata
    group <- droplevels(subset_meta$depth_label)
    permutations <- depth_nonidentity
    scheme <- "exact 215 nonidentity permutations within three pit profiles"
  }

  transformed <- rclr_distance(feature_table, feature_id, subset_meta$sample_id)
  model_data <- data.frame(group = group)
  set.seed(42)
  location_fit <- adonis2(
    transformed$distance ~ group,
    data = model_data,
    permutations = permutations
  )
  dispersion_fit <- betadisper(
    transformed$distance,
    group,
    type = "median",
    bias.adjust = TRUE
  )
  set.seed(42)
  dispersion_test <- permutest(dispersion_fit, permutations = permutations)

  list(
    location = tibble(
      profile = profile_name,
      subset = subset_name,
      transformation = "robust CLR with low-rank zero imputation",
      distance = "Euclidean (Aitchison sensitivity)",
      n_samples = nrow(subset_meta),
      n_features = transformed$n_features,
      rclr_rank = transformed$rclr_rank,
      zero_fraction = transformed$zero_fraction,
      groups = nlevels(group),
      R2 = location_fit$R2[1],
      F = location_fit$F[1],
      p = location_fit$`Pr(>F)`[1],
      permutation_scheme = scheme
    ),
    dispersion = tibble(
      profile = profile_name,
      subset = subset_name,
      transformation = "robust CLR with low-rank zero imputation",
      distance = "Euclidean (Aitchison sensitivity)",
      n_samples = nrow(subset_meta),
      n_features = transformed$n_features,
      rclr_rank = transformed$rclr_rank,
      zero_fraction = transformed$zero_fraction,
      groups = nlevels(group),
      F = dispersion_test$tab[1, "F"],
      p = dispersion_test$tab[1, "Pr(>F)"],
      permutation_scheme = scheme
    )
  )
}

rclr_nested <- imap(profiles, function(profile, profile_name) {
  map(c("surface_site", "nahuelbuta_depth"), ~ run_rclr_test(profile_name, profile, .x))
}) %>% flatten()

rclr_location <- map_dfr(rclr_nested, "location") %>%
  group_by(subset) %>%
  mutate(q_within_subset = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(conclusion = classification(p, q_within_subset)) %>%
  arrange(subset, profile)

rclr_dispersion <- map_dfr(rclr_nested, "dispersion") %>%
  group_by(subset) %>%
  mutate(q_within_subset = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(conclusion = classification(p, q_within_subset)) %>%
  arrange(subset, profile)

write_tsv_checked(rclr_location, file.path(ROBUST_DIR, "rclr_aitchison_location_permanova.tsv"))
write_tsv_checked(rclr_dispersion, file.path(ROBUST_DIR, "rclr_aitchison_dispersion.tsv"))

bray_location <- read_tsv(PATHS[["bray_location"]], show_col_types = FALSE, na = "NA") %>%
  filter(profile %in% names(profiles)) %>%
  transmute(
    profile, subset,
    bray_R2 = R2,
    bray_p = p,
    bray_q = q_within_subset,
    bray_conclusion = conclusion
  )
bray_dispersion <- read_tsv(PATHS[["bray_dispersion"]], show_col_types = FALSE, na = "NA") %>%
  filter(profile %in% names(profiles)) %>%
  transmute(
    profile, subset,
    bray_dispersion_F = F,
    bray_dispersion_p = p,
    bray_dispersion_q = q_within_subset,
    bray_dispersion_conclusion = conclusion
  )

rclr_comparison <- rclr_location %>%
  transmute(
    profile, subset,
    rclr_R2 = R2,
    rclr_p = p,
    rclr_q = q_within_subset,
    rclr_conclusion = conclusion
  ) %>%
  left_join(bray_location, by = c("profile", "subset")) %>%
  left_join(
    rclr_dispersion %>%
      transmute(
        profile, subset,
        rclr_dispersion_F = F,
        rclr_dispersion_p = p,
        rclr_dispersion_q = q_within_subset,
        rclr_dispersion_conclusion = conclusion
      ),
    by = c("profile", "subset")
  ) %>%
  left_join(bray_dispersion, by = c("profile", "subset")) %>%
  mutate(
    location_conclusion_consistent = rclr_conclusion == bray_conclusion,
    dispersion_conclusion_consistent = rclr_dispersion_conclusion == bray_dispersion_conclusion
  )
write_tsv_checked(rclr_comparison, file.path(ROBUST_DIR, "rclr_vs_bray_direct_test_comparison.tsv"))

rm(ko_rpk, mag_rpk, profiles)
invisible(gc())

# ---- Carrier-prevalence completeness sensitivity --------------------------
message("Recomputing carrier prevalence across completeness strata...")
quality <- fread(PATHS[["mag_quality"]])
setnames(quality, "Bin Id", "MAG")
taxonomy <- fread(PATHS[["mag_taxonomy"]])
setnames(taxonomy, "user_genome", "MAG")
assignments <- fread(PATHS[["marker_assignments"]])
registry <- fread(PATHS[["marker_registry"]])
orf_denominators <- fread(PATHS[["dram_annotations"]], select = "fasta")[
  , .(total_ORFs = .N), by = .(MAG = fasta)
]

stopifnot(
  nrow(quality) == 228L,
  nrow(taxonomy) == 228L,
  setequal(quality$MAG, taxonomy$MAG),
  setequal(quality$MAG, orf_denominators$MAG)
)

categories <- unique(
  registry[, .(module, category, module_order, category_order)]
)[order(module_order, category_order)]
assignment_counts <- unique(
  assignments[, .(gene_id, MAG, module, category)]
)[, .(marker_hits = .N), by = .(MAG, module, category)]

mag_universe <- merge(
  quality[, .(MAG, Completeness, Contamination, Total_Coding_Sequences)],
  taxonomy[, .(MAG, Domain, phylum)],
  by = "MAG"
)
mag_universe <- merge(mag_universe, orf_denominators, by = "MAG", all.x = TRUE)
if (anyNA(mag_universe$total_ORFs)) stop("Missing DRAM ORF denominator")

strata <- list(
  all_retained_MAGs = mag_universe,
  completeness_ge80_contamination_lt5 = mag_universe[
    Completeness >= 80 & Contamination < 5
  ],
  completeness_ge90_contamination_lt5 = mag_universe[
    Completeness >= 90 & Contamination < 5
  ]
)

all_grid <- CJ(
  MAG = mag_universe$MAG,
  category_row = seq_len(nrow(categories)),
  unique = TRUE
)
category_key <- copy(categories)
category_key[, category_row := seq_len(.N)]
all_grid <- merge(all_grid, category_key, by = "category_row")
all_grid <- merge(
  all_grid,
  assignment_counts,
  by = c("MAG", "module", "category"),
  all.x = TRUE
)
all_grid[is.na(marker_hits), marker_hits := 0L]
all_grid <- merge(
  all_grid,
  mag_universe[, .(MAG, Domain, phylum, Completeness, Contamination, total_ORFs)],
  by = "MAG"
)
all_grid[, marker_hits_per_1000_ORFs := 1000 * marker_hits / total_ORFs]

carrier_sensitivity <- imap_dfr(strata, function(stratum_mags, stratum_name) {
  all_grid[MAG %in% stratum_mags$MAG,
    .(
      n_MAG = uniqueN(MAG),
      carrier_n = sum(marker_hits > 0),
      carrier_pct = 100 * mean(marker_hits > 0),
      mean_marker_hits_per_1000_ORFs = mean(marker_hits_per_1000_ORFs),
      mean_marker_hits_per_1000_carrier_ORFs = if (
        any(marker_hits > 0)
      ) mean(marker_hits_per_1000_ORFs[marker_hits > 0]) else 0
    ),
    by = .(Domain, phylum, module, category, module_order, category_order)
  ] %>%
    as_tibble() %>%
    mutate(stratum = stratum_name, .before = 1)
}) %>%
  arrange(stratum, desc(n_MAG), phylum, module_order, category_order)

write_tsv_checked(
  carrier_sensitivity,
  file.path(ROBUST_DIR, "MAG_marker_carrier_prevalence_by_quality_stratum.tsv")
)

stratum_counts <- imap_dfr(strata, function(x, name) {
  as_tibble(x) %>%
    summarise(
      stratum = name,
      n_MAG = n(),
      n_bacterial_MAG = sum(Domain == "Bacteria"),
      n_archaeal_MAG = sum(Domain == "Archaea"),
      median_completeness = median(Completeness),
      median_contamination = median(Contamination)
    )
})
write_tsv_checked(stratum_counts, file.path(ROBUST_DIR, "MAG_quality_stratum_counts.tsv"))

focal_phylum_summary <- carrier_sensitivity %>%
  filter(phylum %in% c("Actinomycetota", "Dormibacterota")) %>%
  group_by(stratum, phylum) %>%
  summarise(
    n_MAG = first(n_MAG),
    categories_with_at_least_one_carrier = sum(carrier_n > 0),
    categories_carried_by_at_least_half = sum(carrier_pct >= 50),
    P_acquisition_carrier_pct = carrier_pct[category == "P acquisition"][1],
    .groups = "drop"
  )
write_tsv_checked(
  focal_phylum_summary,
  file.path(ROBUST_DIR, "focal_phylum_carrier_sensitivity_summary.tsv")
)

n_scaffold_audit <- quality[, .(
  n_MAG = .N,
  max_N_contigs = max(N_contigs),
  n_N_contigs_ge1000 = sum(N_contigs >= 1000),
  max_N_scaffolds = max(N_scaffolds),
  n_N_scaffolds_ge1000 = sum(N_scaffolds >= 1000)
)] %>% as_tibble()
write_tsv_checked(n_scaffold_audit, file.path(ROBUST_DIR, "MAG_scaffold_field_audit.tsv"))

# ---- Compact interpretation and reproducibility record --------------------
summary_lines <- c(
  "V43 targeted robustness analysis",
  "",
  "Primary analyses remain the Hellinger-transformed Bray-Curtis tests reported in V42.",
  "The robust-CLR/Aitchison results are a compositional sensitivity analysis and do not replace the prespecified primary endpoint.",
  "Carrier prevalence was recomputed from per-ORF curated marker assignments across all retained MAGs and two stricter completeness/contamination strata.",
  "Table S5A should display N_scaffolds, because the Atlas retention criterion was fewer than 1,000 scaffolds; N_contigs is not the corresponding audit field.",
  "",
  capture.output(print(rclr_comparison)),
  "",
  capture.output(print(focal_phylum_summary))
)
writeLines(summary_lines, file.path(ROBUST_DIR, "V43_robustness_interpretation.txt"), useBytes = TRUE)

writeLines(capture.output(sessionInfo()), file.path(OUT_ROOT, "sessionInfo.txt"), useBytes = TRUE)

manifest <- tibble(
  file = sort(list.files(OUT_ROOT, recursive = TRUE, full.names = TRUE)),
  size_bytes = file.info(file)$size
) %>%
  mutate(file = sub(paste0("^", OUT_ROOT, "/"), "", file))
write_tsv_checked(manifest, file.path(OUT_ROOT, "V43_analysis_manifest.tsv"))

message("V43 targeted robustness analysis completed: ", OUT_ROOT)
