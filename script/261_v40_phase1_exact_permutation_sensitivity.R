# =============================================================================
# V40 Phase 1: exact restricted permutations, environmental sensitivities,
# and TPM / feature-count reproducibility audit.
#
# V39 remains read-only. All outputs are written below the V40 execution folder.
# No figure or manuscript file is generated or modified.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(permute)
  library(tidyverse)
  library(vegan)
})

options(scipen = 999)
set.seed(42)

ROOT <- normalizePath(Sys.getenv("CHILE_PROJECT_ROOT", getwd()), mustWork = TRUE)
setwd(ROOT)
V40 <- file.path(ROOT, "legacy_reference/V40_targeted_revision_2026-08-26")
SUPPORT <- file.path(V40, "phase1_support")
PERM_DIR <- file.path(SUPPORT, "permutations")
STAT_DIR <- file.path(SUPPORT, "statistics")
SENS_DIR <- file.path(SUPPORT, "missingness_sensitivity")
TPM_DIR <- file.path(SUPPORT, "tpm_feature_count")
for (path in c(PERM_DIR, STAT_DIR, SENS_DIR, TPM_DIR)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

PATHS <- c(
  reads = "data/processed/metagenome_output/Genecatalog/Nmapped_reads.tsv",
  lengths = "data/processed/metagenome_output/Genecatalog/seq_len_nucl.tsv",
  annotation = "data/processed/metagenome_output/Genecatalog/eggNog.tsv",
  mag_counts = "data/processed/metagenome_output/MAGs/counts_genomes.tsv",
  mag_quality = "data/processed/metagenome_output/MAGs/genome_quality.tsv",
  metadata = "data/metadata/env_2024_meta.tsv",
  marker_tpm = paste0(
    "analysis/targeted_scientific_restructuring_2026-07-13/",
    "curated_module_validation/curated_feature_sample_TPM.tsv"
  ),
  old_module_tests = paste0(
    "analysis/targeted_scientific_restructuring_2026-07-13/",
    "curated_module_validation/curated_module_PERMANOVA_dispersion.tsv"
  ),
  current_direct_location = paste0(
    "analysis/V39_supplementary_evidence_pruning_2026-08-25/",
    "source_data/direct_location_permanova.tsv"
  ),
  current_direct_dispersion = paste0(
    "analysis/V39_supplementary_evidence_pruning_2026-08-25/",
    "source_data/direct_location_dispersion.tsv"
  ),
  current_soil_axis = paste0(
    "analysis/V05_comment_review_2026-08-17/figure3_cross_resolution/",
    "cross_resolution_soil_axis_models.tsv"
  ),
  current_gene_ko_procrustes = paste0(
    "analysis/V39_supplementary_evidence_pruning_2026-08-25/",
    "source_data/whole_gene_ko_procrustes.tsv"
  ),
  current_function_taxonomy = paste0(
    "analysis/V05_comment_review_2026-08-17/figure3_cross_resolution/",
    "function_taxonomy_cross_resolution.tsv"
  )
)

missing_inputs <- PATHS[!file.exists(PATHS)]
if (length(missing_inputs)) {
  stop("Missing V40 Phase 1 input(s): ", paste(missing_inputs, collapse = ", "))
}

SAMPLES <- c("C01", sprintf("C%02d", 3:18))
ENV_VARS <- c(
  "pH", "Conductivity", "moisture", "C", "N", "Feo",
  "Alo", "Mno", "Sio", "NH4", "NO3", "Pt"
)
NB2_MISSING_VARS <- c("C", "N", "Feo", "Alo", "Sio", "Pt")
COMPLETE_VARS <- setdiff(ENV_VARS, NB2_MISSING_VARS)
DEPTH_LEVELS <- c("5_10cm", "20_40cm", "60_80cm")
PROFILE_LEVELS <- c(
  "Whole gene catalog",
  "Fractional KO catalog",
  "Genome-size-normalised MAG assemblage"
)

write_tsv_checked <- function(x, path) {
  readr::write_tsv(x, path, na = "NA")
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Failed to write output: ", path)
  }
  invisible(path)
}

subset_dist <- function(distance, sample_ids) {
  matrix <- as.matrix(distance)
  as.dist(matrix[sample_ids, sample_ids, drop = FALSE])
}

classification <- function(p, q = p) {
  case_when(
    is.na(p) ~ "effect size only",
    q < 0.05 ~ "FDR-supported",
    p < 0.05 ~ "nominal only; not FDR-supported",
    TRUE ~ "not nominally supported"
  )
}

# ---- Canonical metadata and design validation -------------------------------
metadata <- read_tsv(PATHS[["metadata"]], show_col_types = FALSE, na = "") %>%
  filter(sample_id %in% SAMPLES) %>%
  arrange(match(sample_id, SAMPLES)) %>%
  mutate(
    site = factor(site, levels = c("AZ", "SG", "LC", "NB")),
    depth_label = factor(depth_label, levels = DEPTH_LEVELS),
    site_pit = factor(site_pit)
  )

if (!identical(metadata$sample_id, SAMPLES)) {
  stop("Canonical metadata do not match the frozen 17-metagenome order.")
}
if (!all(metadata$dna_type == "iDNA")) {
  stop("Canonical metadata do not identify all 17 metagenomes as iDNA.")
}

depth_metadata <- metadata %>%
  filter(site == "NB") %>%
  mutate(site_pit = droplevels(site_pit), depth_label = droplevels(depth_label))
expected_depth_order <- c("C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18")
if (!identical(depth_metadata$sample_id, expected_depth_order)) {
  stop("Unexpected Nahuelbuta depth-sample order.")
}
if (!all(table(depth_metadata$site_pit, depth_metadata$depth_label) == 1L)) {
  stop("Nahuelbuta is not a complete 3-pit by 3-depth block design.")
}

# ---- Explicit finite permutation spaces ------------------------------------
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
  colnames(full) <- paste0("position_", seq_len(ncol(full)))
  full
}

write_permutation_matrix <- function(matrix, sample_ids, path) {
  observed <- apply(matrix, 1, function(x) all(unname(x) == seq_along(sample_ids)))
  table <- as_tibble(matrix) %>%
    mutate(
      permutation_id = seq_len(n()),
      is_observed_arrangement = observed,
      .before = 1
    )
  write_tsv_checked(table, path)
}

depth_full <- enumerate_block_permutations(depth_metadata$site_pit)
depth_identity <- apply(
  depth_full,
  1,
  function(x) all(unname(x) == seq_len(nrow(depth_metadata)))
)
depth_nonidentity <- depth_full[!depth_identity, , drop = FALSE]

depth_control_999 <- how(nperm = 999)
setBlocks(depth_control_999) <- depth_metadata$site_pit
requested_generated <- shuffleSet(
  n = nrow(depth_metadata),
  nset = 999,
  control = depth_control_999
)

matrix_keys <- function(matrix) apply(matrix, 1, paste, collapse = ",")
depth_position_key <- depth_metadata %>%
  transmute(
    position = row_number(),
    sample_id,
    site_pit = as.character(site_pit),
    depth_label = as.character(depth_label)
  )

write_tsv_checked(depth_position_key, file.path(PERM_DIR, "pit3_position_key.tsv"))
write_permutation_matrix(
  depth_full,
  depth_metadata$sample_id,
  file.path(PERM_DIR, "pit3_complete_space_216.tsv")
)
write_permutation_matrix(
  depth_nonidentity,
  depth_metadata$sample_id,
  file.path(PERM_DIR, "pit3_exact_nonidentity_215.tsv")
)
write_permutation_matrix(
  requested_generated,
  depth_metadata$sample_id,
  file.path(PERM_DIR, "pit3_requested999_actual_matrix.tsv")
)

permutation_audit <- tibble(
  design = "three complete Nahuelbuta pit profiles",
  requested_permutations = 999L,
  theoretical_arrangements_including_observed = nrow(depth_full),
  theoretical_nonidentity_permutations = nrow(depth_nonidentity),
  generated_rows_from_how_shuffleSet = nrow(requested_generated),
  unique_generated_rows = n_distinct(matrix_keys(requested_generated)),
  generated_contains_observed = any(matrix_keys(requested_generated) == paste(seq_len(9), collapse = ",")),
  duplicate_generated_rows = nrow(requested_generated) - n_distinct(matrix_keys(requested_generated)),
  generated_set_equals_explicit_nonidentity_set = setequal(
    matrix_keys(requested_generated),
    matrix_keys(depth_nonidentity)
  )
)

if (!identical(nrow(depth_full), 216L) || !identical(nrow(depth_nonidentity), 215L)) {
  stop("The explicit three-pit permutation space is not 216/215.")
}
if (!all(permutation_audit$unique_generated_rows == 215L) ||
    permutation_audit$generated_contains_observed ||
    permutation_audit$duplicate_generated_rows != 0L ||
    !permutation_audit$generated_set_equals_explicit_nonidentity_set) {
  stop("The requested-999 restricted permutation audit failed.")
}
write_tsv_checked(permutation_audit, file.path(PERM_DIR, "pit3_permutation_audit.tsv"))

# A separate exact space is needed after excluding NB2: 3! x 3! = 36.
depth_without_nb2 <- depth_metadata %>% filter(site_pit != "NB2") %>% droplevels()
depth_two_pit_full <- enumerate_block_permutations(depth_without_nb2$site_pit)
depth_two_pit_identity <- apply(
  depth_two_pit_full,
  1,
  function(x) all(unname(x) == seq_len(nrow(depth_without_nb2)))
)
depth_two_pit_nonidentity <- depth_two_pit_full[!depth_two_pit_identity, , drop = FALSE]
if (nrow(depth_two_pit_full) != 36L || nrow(depth_two_pit_nonidentity) != 35L) {
  stop("The explicit two-pit sensitivity space is not 36/35.")
}
write_tsv_checked(
  depth_without_nb2 %>%
    transmute(
      position = row_number(), sample_id,
      site_pit = as.character(site_pit), depth_label = as.character(depth_label)
    ),
  file.path(PERM_DIR, "pit2_exclude_NB2_position_key.tsv")
)
write_permutation_matrix(
  depth_two_pit_full,
  depth_without_nb2$sample_id,
  file.path(PERM_DIR, "pit2_exclude_NB2_complete_space_36.tsv")
)
write_permutation_matrix(
  depth_two_pit_nonidentity,
  depth_without_nb2$sample_id,
  file.path(PERM_DIR, "pit2_exclude_NB2_exact_nonidentity_35.tsv")
)

# ---- Gene catalog, TPM audit, KO profile and MAG profile --------------------
message("Loading the 2.42-million-feature gene catalog once...")
counts <- fread(PATHS[["reads"]], select = c("GeneName", SAMPLES))
lengths <- fread(
  PATHS[["lengths"]],
  header = FALSE,
  col.names = c("GeneName", "Length")
)
annotation <- fread(PATHS[["annotation"]], select = c("GeneName", "KO"))

catalog_gene_clusters <- nrow(counts)
duplicate_count_ids <- counts[, .N, by = GeneName][N > 1, .N]
duplicate_length_ids <- lengths[, .N, by = GeneName][N > 1, .N]

gene_with_length <- merge(counts, lengths, by = "GeneName", all.x = TRUE)
unmatched_length_features <- sum(is.na(gene_with_length$Length))
invalid_length_features <- sum(!is.na(gene_with_length$Length) &
                                 (!is.finite(gene_with_length$Length) | gene_with_length$Length <= 0))
gene_with_length[, total_mapped_reads := rowSums(.SD, na.rm = TRUE), .SDcols = SAMPLES]
all_zero_features <- sum(gene_with_length$total_mapped_reads == 0)
nonzero_valid_features <- sum(
  gene_with_length$total_mapped_reads > 0 &
    is.finite(gene_with_length$Length) & gene_with_length$Length > 0
)

gene_rpk <- gene_with_length[
  total_mapped_reads > 0 & is.finite(Length) & Length > 0
]
for (sample_id in SAMPLES) {
  set(
    gene_rpk,
    j = sample_id,
    value = gene_rpk[[sample_id]] / (gene_rpk$Length / 1000)
  )
}
rpk_library_sums <- colSums(as.matrix(gene_rpk[, ..SAMPLES]), na.rm = TRUE)
tpm_library_sums <- rpk_library_sums / (rpk_library_sums / 1e6)

tpm_feature_audit <- tibble(
  metric = c(
    "catalog_gene_clusters",
    "count_rows_with_valid_length",
    "unmatched_length_features",
    "invalid_nonpositive_length_features",
    "all_zero_or_unmapped_features_across_17_samples",
    "composition_matrix_features_after_all_zero_filter",
    "duplicate_gene_ids_in_count_table",
    "duplicate_gene_ids_in_length_table"
  ),
  value = c(
    catalog_gene_clusters,
    sum(is.finite(gene_with_length$Length) & gene_with_length$Length > 0),
    unmatched_length_features,
    invalid_length_features,
    all_zero_features,
    nonzero_valid_features,
    duplicate_count_ids,
    duplicate_length_ids
  ),
  evidence = c(
    "row count of Nmapped_reads.tsv",
    "GeneName join to seq_len_nucl.tsv with finite positive length",
    "count-table GeneName without a length record",
    "joined length is non-finite or <= 0",
    "row sum of mapped reads equals zero across all 17 metagenomes",
    "valid-length features with at least one mapped read",
    "duplicated GeneName groups in Nmapped_reads.tsv",
    "duplicated GeneName groups in seq_len_nucl.tsv"
  )
)

if (catalog_gene_clusters != 2422720L || all_zero_features != 170L ||
    nonzero_valid_features != 2422550L) {
  stop("The 2,422,720 -> 2,422,550 feature-count audit did not reproduce.")
}
if (unmatched_length_features != 0L || invalid_length_features != 0L ||
    duplicate_count_ids != 0L || duplicate_length_ids != 0L) {
  stop("Unexpected GeneName/length integrity problem.")
}
if (max(abs(tpm_library_sums - 1e6)) > 1e-6) {
  stop("TPM library sums do not equal one million within tolerance.")
}

write_tsv_checked(tpm_feature_audit, file.path(TPM_DIR, "tpm_feature_count_audit.tsv"))
write_tsv_checked(
  tibble(
    sample_id = SAMPLES,
    rpk_library_sum = as.numeric(rpk_library_sums),
    implied_tpm_library_sum = as.numeric(tpm_library_sums),
    absolute_error_from_one_million = abs(as.numeric(tpm_library_sums) - 1e6)
  ),
  file.path(TPM_DIR, "tpm_library_sum_validation.tsv")
)
writeLines(
  c(
    "For each gene cluster and sample, mapped reads were divided by gene length in kilobases to obtain RPK.",
    "Each sample's RPK values were then divided by the sample RPK sum and multiplied by 1,000,000.",
    "The resulting quantity is mathematically TPM, but because it derives from metagenomic DNA it should be described as length- and library-normalised DNA gene abundance, not transcript expression.",
    "The catalog contains 2,422,720 gene clusters; 170 had zero mapped reads across all 17 metagenomes, leaving 2,422,550 features in composition analyses."
  ),
  file.path(TPM_DIR, "methods_table_note_ready_text.txt"),
  useBytes = TRUE
)

profile_from_rpk <- function(table, label) {
  matrix <- as.matrix(table[, ..SAMPLES])
  storage.mode(matrix) <- "double"
  community <- decostand(t(matrix), method = "hellinger")
  rownames(community) <- SAMPLES
  distance <- vegdist(community, method = "bray")
  list(label = label, n_features = nrow(table), distance = distance)
}

whole_gene_profile <- profile_from_rpk(gene_rpk, "Whole gene catalog")

message("Aggregating fractionally assigned KO abundances...")
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
ko_rpk <- ko_joined[
  , lapply(.SD, sum, na.rm = TRUE),
  by = KO,
  .SDcols = SAMPLES
]
ko_profile <- profile_from_rpk(ko_rpk, "Fractional KO catalog")

ko_feature_audit <- tibble(
  metric = c(
    "genes_with_at_least_one_valid_KO",
    "genes_with_multiple_valid_KOs",
    "detected_fractional_KO_features"
  ),
  value = c(
    uniqueN(ko_map$GeneName),
    uniqueN(ko_map[n_valid_KO > 1]$GeneName),
    ko_profile$n_features
  )
)
write_tsv_checked(ko_feature_audit, file.path(TPM_DIR, "ko_feature_count_audit.tsv"))

rm(annotation, counts, lengths, gene_with_length, ko_joined, ko_map, ko_rpk, gene_rpk)
invisible(gc())

message("Building the genome-size-normalised MAG profile...")
mag_counts <- fread(PATHS[["mag_counts"]])
setnames(mag_counts, "Sample", "MAG")
mag_quality <- fread(PATHS[["mag_quality"]], select = c("Bin Id", "Genome_Size"))
setnames(mag_quality, "Bin Id", "MAG")
mag_rpk <- merge(mag_counts, mag_quality, by = "MAG", all.x = TRUE)
if (anyNA(mag_rpk$Genome_Size) || any(mag_rpk$Genome_Size <= 0)) {
  stop("Missing or invalid MAG genome size.")
}
for (sample_id in SAMPLES) {
  set(
    mag_rpk,
    j = sample_id,
    value = mag_rpk[[sample_id]] / (mag_rpk$Genome_Size / 1000)
  )
}
mag_profile <- profile_from_rpk(
  mag_rpk[, c("MAG", SAMPLES), with = FALSE],
  "Genome-size-normalised MAG assemblage"
)
rm(mag_counts, mag_quality, mag_rpk)
invisible(gc())

profiles <- list(whole_gene_profile, ko_profile, mag_profile)

# ---- Exact direct PERMANOVA and PERMDISP -----------------------------------
run_direct_test <- function(profile, subset) {
  profile_label <- profile$label
  profile_features <- profile$n_features
  if (subset == "surface_site") {
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
  distance <- subset_dist(profile$distance, subset_meta$sample_id)
  model_data <- data.frame(group = group)
  set.seed(42)
  location_fit <- adonis2(distance ~ group, data = model_data, permutations = permutations)
  dispersion_fit <- betadisper(distance, group, type = "median", bias.adjust = TRUE)
  set.seed(42)
  dispersion_test <- permutest(dispersion_fit, permutations = permutations)
  list(
    location = tibble(
      profile = profile_label,
      subset,
      n_samples = nrow(subset_meta),
      n_features = profile_features,
      groups = nlevels(group),
      R2 = location_fit$R2[1],
      F = location_fit$F[1],
      p = location_fit$`Pr(>F)`[1],
      permutation_scheme = scheme
    ),
    dispersion = tibble(
      profile = profile_label,
      subset,
      n_samples = nrow(subset_meta),
      groups = nlevels(group),
      F = dispersion_test$tab[1, "F"],
      p = dispersion_test$tab[1, "Pr(>F)"],
      permutation_scheme = scheme
    )
  )
}

direct_nested <- map(profiles, function(profile) {
  map(c("surface_site", "nahuelbuta_depth"), ~ run_direct_test(profile, .x))
}) %>% flatten()

direct_location <- map_dfr(direct_nested, "location") %>%
  mutate(
    profile = factor(profile, levels = PROFILE_LEVELS),
    q_within_subset = ave(p, subset, FUN = function(x) p.adjust(x, method = "BH"))
  ) %>%
  arrange(subset, profile) %>%
  mutate(profile = as.character(profile), conclusion = classification(p, q_within_subset))
direct_dispersion <- map_dfr(direct_nested, "dispersion") %>%
  mutate(
    profile = factor(profile, levels = PROFILE_LEVELS),
    q_within_subset = ave(p, subset, FUN = function(x) p.adjust(x, method = "BH"))
  ) %>%
  arrange(subset, profile) %>%
  mutate(profile = as.character(profile), conclusion = classification(p, q_within_subset))

write_tsv_checked(direct_location, file.path(STAT_DIR, "exact_direct_location_permanova.tsv"))
write_tsv_checked(direct_dispersion, file.path(STAT_DIR, "exact_direct_location_dispersion.tsv"))

current_location <- read_tsv(PATHS[["current_direct_location"]], show_col_types = FALSE, na = "NA")
current_dispersion <- read_tsv(PATHS[["current_direct_dispersion"]], show_col_types = FALSE, na = "NA")
direct_comparison <- bind_rows(
  direct_location %>%
    select(profile, subset, statistic_new = F, p_new = p, q_new = q_within_subset, conclusion_new = conclusion) %>%
    left_join(
      current_location %>%
        transmute(
          profile, subset, statistic_current_v39 = F,
          p_current_v39 = p, q_current_v39 = q_within_subset,
          conclusion_current_v39 = classification(p, q_within_subset)
        ),
      by = c("profile", "subset")
    ) %>%
    mutate(test = "PERMANOVA", .before = 1),
  direct_dispersion %>%
    select(profile, subset, statistic_new = F, p_new = p, q_new = q_within_subset, conclusion_new = conclusion) %>%
    left_join(
      current_dispersion %>%
        transmute(
          profile, subset, statistic_current_v39 = F,
          p_current_v39 = p, q_current_v39 = q_within_subset,
          conclusion_current_v39 = classification(p, q_within_subset)
        ),
      by = c("profile", "subset")
    ) %>%
    mutate(test = "PERMDISP", .before = 1)
) %>%
  mutate(
    statistic_difference = statistic_new - statistic_current_v39,
    p_difference = p_new - p_current_v39,
    q_difference = q_new - q_current_v39,
    scientific_conclusion_changed = conclusion_new != conclusion_current_v39
  )
write_tsv_checked(direct_comparison, file.path(STAT_DIR, "exact_direct_vs_v39_comparison.tsv"))

# ---- Environmental scenarios and whole-gene/KO dbRDA -----------------------
impute_depth_matched <- function(data) {
  out <- data
  log <- list()
  for (variable in ENV_VARS) {
    for (index in which(is.na(out[[variable]]))) {
      if (as.character(out$site[index]) != "NB") {
        stop("Primary depth-matched imputation encountered non-Nahuelbuta missingness.")
      }
      donors <- out %>%
        filter(
          site == "NB",
          depth_label == out$depth_label[index],
          !is.na(.data[[variable]])
        ) %>%
        pull(all_of(variable))
      value <- mean(donors)
      out[[variable]][index] <- value
      log[[length(log) + 1L]] <- tibble(
        sample_id = out$sample_id[index],
        site_pit = as.character(out$site_pit[index]),
        depth_label = as.character(out$depth_label[index]),
        variable,
        imputed_value = value,
        donor_rule = "other observed Nahuelbuta pits at the same depth"
      )
    }
  }
  if (anyNA(out[, ENV_VARS])) stop("Primary imputation left missing values.")
  list(data = out, log = bind_rows(log))
}

observed_missingness <- map_dfr(ENV_VARS, function(variable) {
  tibble(
    variable,
    n_missing = sum(is.na(metadata[[variable]])),
    missing_samples = paste(metadata$sample_id[is.na(metadata[[variable]])], collapse = ";")
  )
})
write_tsv_checked(observed_missingness, file.path(SENS_DIR, "environment_missingness.tsv"))

expected_missing <- observed_missingness %>% filter(n_missing > 0)
if (!setequal(expected_missing$variable, NB2_MISSING_VARS) ||
    !all(expected_missing$n_missing == 3L) ||
    !all(expected_missing$missing_samples == "C11;C14;C17")) {
  stop("Environmental missingness does not match the pre-specified NB2 block.")
}

primary_imputation <- impute_depth_matched(metadata)
write_tsv_checked(primary_imputation$log, file.path(SENS_DIR, "primary_depth_matched_imputation_log.tsv"))

scenario_inputs <- list(
  primary_depth_matched = list(data = primary_imputation$data, variables = ENV_VARS),
  complete_variable = list(data = metadata, variables = COMPLETE_VARS),
  exclude_NB2 = list(
    data = metadata %>% filter(as.character(site_pit) != "NB2") %>% droplevels(),
    variables = ENV_VARS
  )
)

fit_scenario <- function(specification, name) {
  data <- specification$data
  variables <- specification$variables
  if (anyNA(data[, variables])) {
    stop("Scenario ", name, " contains missing predictor values.")
  }
  pca <- prcomp(data[, variables], center = TRUE, scale. = TRUE)
  variance <- pca$sdev^2 / sum(pca$sdev^2)
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  colnames(scores) <- c("SoilPC1", "SoilPC2")
  rownames(scores) <- data$sample_id

  scores_output <- bind_cols(
    tibble(
      scenario = name,
      sample_id = data$sample_id,
      site = as.character(data$site),
      site_pit = as.character(data$site_pit),
      depth_label = as.character(data$depth_label)
    ),
    as_tibble(scores)
  )
  loadings_output <- as_tibble(pca$rotation, rownames = "variable") %>%
    mutate(scenario = name, .before = 1)

  fit_profile <- function(profile, subset) {
    subset_data <- if (subset == "all") {
      data
    } else if (subset == "surface") {
      data %>% filter(depth_label == "5_10cm")
    } else {
      data %>% filter(site == "NB")
    }
    sample_ids <- subset_data$sample_id
    response_distance <- subset_dist(profile$distance, sample_ids)
    predictors <- scores[sample_ids, , drop = FALSE]
    model <- dbrda(
      response_distance ~ SoilPC1 + SoilPC2,
      data = predictors
    )
    r2 <- RsquareAdj(model)
    p_value <- NA_real_
    scheme <- "effect size only; no mixed-design permutation"
    exact_rows <- NA_integer_
    if (subset == "surface") {
      set.seed(42)
      p_value <- anova.cca(model, permutations = 999L)$`Pr(>F)`[1]
      scheme <- paste0("999 unrestricted permutations among ", nrow(subset_data), " surface samples")
      exact_rows <- 999L
    } else if (subset == "depth") {
      if (name == "exclude_NB2") {
        permutations <- depth_two_pit_nonidentity
      } else {
        permutations <- depth_nonidentity
      }
      set.seed(42)
      p_value <- anova.cca(model, permutations = permutations)$`Pr(>F)`[1]
      scheme <- paste0(
        "exact ", nrow(permutations),
        " nonidentity permutations within ", nlevels(droplevels(subset_data$site_pit)),
        " pit profiles"
      )
      exact_rows <- nrow(permutations)
    }
    tibble(
      scenario = name,
      profile = profile$label,
      subset,
      n_samples = nrow(subset_data),
      n_environment_variables = length(variables),
      variables = paste(variables, collapse = ";"),
      constrained_R2 = r2$r.squared,
      adjusted_R2 = r2$adj.r.squared,
      p = p_value,
      PC1_variance = variance[1],
      PC2_variance = variance[2],
      permutation_scheme = scheme,
      permutation_rows = exact_rows
    )
  }

  models <- map_dfr(list(whole_gene_profile, ko_profile), function(profile) {
    map_dfr(c("all", "surface", "depth"), ~ fit_profile(profile, .x))
  })
  list(scores = scores_output, loadings = loadings_output, models = models)
}

scenario_results <- imap(scenario_inputs, fit_scenario)
scenario_scores <- map_dfr(scenario_results, "scores")
scenario_loadings <- map_dfr(scenario_results, "loadings")
scenario_models <- map_dfr(scenario_results, "models") %>%
  group_by(scenario, subset) %>%
  mutate(q_across_gene_and_KO = if (all(is.na(p))) {
    rep(NA_real_, n())
  } else {
    p.adjust(p, method = "BH")
  }) %>%
  ungroup() %>%
  mutate(conclusion = classification(p, q_across_gene_and_KO))

write_tsv_checked(scenario_scores, file.path(SENS_DIR, "soil_pca_scores_by_scenario.tsv"))
write_tsv_checked(scenario_loadings, file.path(SENS_DIR, "soil_pca_loadings_by_scenario.tsv"))
write_tsv_checked(scenario_models, file.path(SENS_DIR, "whole_gene_KO_dbrda_sensitivity.tsv"))

primary_scores <- scenario_scores %>%
  filter(scenario == "primary_depth_matched") %>%
  select(sample_id, SoilPC1, SoilPC2)
ordination_comparison <- map_dfr(setdiff(names(scenario_inputs), "primary_depth_matched"), function(name) {
  current_scores <- scenario_scores %>%
    filter(scenario == name) %>%
    select(sample_id, SoilPC1, SoilPC2)
  common <- intersect(primary_scores$sample_id, current_scores$sample_id)
  primary_matrix <- primary_scores %>%
    filter(sample_id %in% common) %>%
    arrange(match(sample_id, common)) %>%
    select(SoilPC1, SoilPC2) %>% as.matrix()
  current_matrix <- current_scores %>%
    filter(sample_id %in% common) %>%
    arrange(match(sample_id, common)) %>%
    select(SoilPC1, SoilPC2) %>% as.matrix()
  fit <- procrustes(primary_matrix, current_matrix, symmetric = TRUE)
  tibble(
    scenario = name,
    n_common_samples = length(common),
    two_axis_procrustes_r = sqrt(1 - fit$ss),
    spearman_pairwise_distance_r = cor(
      as.vector(dist(primary_matrix)),
      as.vector(dist(current_matrix)),
      method = "spearman"
    ),
    interpretation = "PCA sign and rotation are arbitrary; Procrustes and pairwise-distance correlation assess sample ordering."
  )
})
write_tsv_checked(
  ordination_comparison,
  file.path(SENS_DIR, "soil_pca_ordination_ordering_comparison.tsv")
)

current_axis <- read_tsv(PATHS[["current_soil_axis"]], show_col_types = FALSE, na = "NA")
primary_vs_current <- scenario_models %>%
  filter(scenario == "primary_depth_matched") %>%
  mutate(
    subset_current = recode(
      subset,
      all = "all_17",
      surface = "surface_11",
      depth = "nahuelbuta_depth_9"
    )
  ) %>%
  left_join(
    current_axis %>%
      transmute(
        profile, subset_current = subset,
        current_R2 = R2,
        current_adjusted_R2 = adjusted_R2,
        current_p = p
      ),
    by = c("profile", "subset_current")
  ) %>%
  mutate(
    R2_difference = constrained_R2 - current_R2,
    adjusted_R2_difference = adjusted_R2 - current_adjusted_R2,
    p_difference = p - current_p
  )
write_tsv_checked(
  primary_vs_current,
  file.path(SENS_DIR, "primary_soil_axis_vs_current_values.tsv")
)

primary_conclusions <- scenario_models %>%
  filter(scenario == "primary_depth_matched") %>%
  select(profile, subset, primary_conclusion = conclusion)
sensitivity_conclusion_comparison <- scenario_models %>%
  filter(scenario != "primary_depth_matched") %>%
  left_join(primary_conclusions, by = c("profile", "subset")) %>%
  mutate(conclusion_changed_from_primary = conclusion != primary_conclusion)
write_tsv_checked(
  sensitivity_conclusion_comparison,
  file.path(SENS_DIR, "sensitivity_conclusion_comparison.tsv")
)

# ---- Exact module PERMANOVA and PERMDISP -----------------------------------
message("Recomputing all marker-module depth tests with the exact 215-row set...")
feature_tpm <- read_tsv(PATHS[["marker_tpm"]], show_col_types = FALSE, na = "")
old_module_tests <- read_tsv(PATHS[["old_module_tests"]], show_col_types = FALSE, na = "")
module_order <- c(
  "nitrogen_transformations",
  "pedogenic_resource_acquisition",
  "alternative_energy_chemolithotrophy",
  "stress_response",
  "cazy_substrate_processing"
)

run_exact_module_depth <- function(variant, module) {
  dat <- feature_tpm %>%
    filter(
      analysis_variant == variant,
      .data$module == !!module,
      sample_id %in% depth_metadata$sample_id
    )
  comm <- dat %>%
    select(sample_id, marker, TPM) %>%
    group_by(sample_id, marker) %>%
    summarise(TPM = sum(TPM), .groups = "drop") %>%
    pivot_wider(names_from = marker, values_from = TPM, values_fill = 0) %>%
    arrange(match(sample_id, depth_metadata$sample_id))
  if (!identical(comm$sample_id, depth_metadata$sample_id)) {
    stop("Module matrix lost a Nahuelbuta sample: ", variant, " / ", module)
  }
  comm_matrix <- as.matrix(select(comm, -sample_id))
  comm_matrix <- comm_matrix[, colSums(comm_matrix) > 0, drop = FALSE]
  distance <- vegdist(decostand(comm_matrix, method = "hellinger"), method = "bray")
  group <- depth_metadata$depth_label
  set.seed(42)
  ad <- adonis2(
    distance ~ group,
    data = data.frame(group = group),
    permutations = depth_nonidentity
  )
  dispersion <- betadisper(distance, group, type = "median", add = "lingoes")
  set.seed(42)
  dispersion_test <- permutest(dispersion, permutations = depth_nonidentity)
  old_row <- old_module_tests %>%
    filter(analysis_variant == variant, .data$module == !!module, subset == "NA_depth") %>%
    slice(1)
  tibble(
    analysis_variant = variant,
    predefined_primary = variant == "primary_marker_fractional",
    module,
    tier = old_row$tier,
    subset = "NA_depth",
    n_samples = 9L,
    group_sizes = "surface:3;middle:3;deep:3",
    n_features_registered = old_row$n_features_registered,
    n_features_detected = ncol(comm_matrix),
    distance = "Bray-Curtis on Hellinger-transformed TPM",
    PERMANOVA_Df = ad$Df[1],
    PERMANOVA_R2 = ad$R2[1],
    PERMANOVA_F = ad$F[1],
    PERMANOVA_p = ad$`Pr(>F)`[1],
    dispersion_F = dispersion_test$tab[1, "F"],
    dispersion_p = dispersion_test$tab[1, "Pr(>F)"],
    dispersion_distance_correction = "Lingoes",
    permutation_scheme = "exact 215 nonidentity permutations within site_pit",
    status = "complete"
  )
}

exact_module_depth <- crossing(
  analysis_variant = unique(feature_tpm$analysis_variant),
  module = module_order
) %>%
  pmap_dfr(~ run_exact_module_depth(..1, ..2))

surface_module_rows <- old_module_tests %>%
  filter(subset == "surface_climate") %>%
  mutate(permutation_scheme = "999 unrestricted permutations") %>%
  select(any_of(names(exact_module_depth)))

module_exact_combined <- bind_rows(surface_module_rows, exact_module_depth) %>%
  group_by(analysis_variant, tier) %>%
  mutate(
    PERMANOVA_FDR_within_tier = p.adjust(PERMANOVA_p, method = "BH"),
    dispersion_FDR_within_tier = p.adjust(dispersion_p, method = "BH")
  ) %>%
  ungroup() %>%
  group_by(analysis_variant) %>%
  mutate(
    PERMANOVA_FDR_global = p.adjust(PERMANOVA_p, method = "BH"),
    dispersion_FDR_global = p.adjust(dispersion_p, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    module = factor(module, levels = module_order),
    location_conclusion = classification(PERMANOVA_p, PERMANOVA_FDR_within_tier),
    dispersion_conclusion = classification(dispersion_p, dispersion_FDR_within_tier)
  ) %>%
  arrange(analysis_variant, subset, module) %>%
  mutate(module = as.character(module))

write_tsv_checked(
  module_exact_combined,
  file.path(STAT_DIR, "exact_curated_module_PERMANOVA_dispersion.tsv")
)

module_comparison <- exact_module_depth %>%
  left_join(
    old_module_tests %>%
      filter(subset == "NA_depth") %>%
      select(
        analysis_variant, module,
        old_R2 = PERMANOVA_R2,
        old_F = PERMANOVA_F,
        old_p = PERMANOVA_p,
        old_dispersion_F = dispersion_F,
        old_dispersion_p = dispersion_p
      ),
    by = c("analysis_variant", "module")
  ) %>%
  mutate(
    R2_difference = PERMANOVA_R2 - old_R2,
    F_difference = PERMANOVA_F - old_F,
    p_difference = PERMANOVA_p - old_p,
    dispersion_F_difference = dispersion_F - old_dispersion_F,
    dispersion_p_difference = dispersion_p - old_dispersion_p
  )
write_tsv_checked(
  module_comparison,
  file.path(STAT_DIR, "exact_module_depth_vs_original_unrestricted_source.tsv")
)

# ---- Exact cross-resolution Procrustes tests -------------------------------
run_exact_procrustes <- function(profile_x, profile_y) {
  x_distance <- subset_dist(profile_x$distance, depth_metadata$sample_id)
  y_distance <- subset_dist(profile_y$distance, depth_metadata$sample_id)
  x_coordinates <- wcmdscale(x_distance, k = 2, eig = TRUE)$points[, 1:2, drop = FALSE]
  y_coordinates <- wcmdscale(y_distance, k = 2, eig = TRUE)$points[, 1:2, drop = FALSE]
  set.seed(42)
  test <- protest(
    x_coordinates,
    y_coordinates,
    permutations = depth_nonidentity,
    symmetric = TRUE
  )
  tibble(
    profile_x = profile_x$label,
    profile_y = profile_y$label,
    subset = "nahuelbuta_depth_9",
    n_samples = 9L,
    procrustes_r = unname(test$t0),
    p = test$signif,
    permutation_scheme = "exact 215 nonidentity permutations within site_pit"
  )
}

exact_procrustes <- bind_rows(
  run_exact_procrustes(whole_gene_profile, ko_profile),
  run_exact_procrustes(whole_gene_profile, mag_profile),
  run_exact_procrustes(ko_profile, mag_profile)
) %>%
  mutate(
    q_across_three_pairs = p.adjust(p, method = "BH"),
    conclusion = classification(p, q_across_three_pairs)
  )
write_tsv_checked(exact_procrustes, file.path(STAT_DIR, "exact_depth_procrustes.tsv"))

current_procrustes <- bind_rows(
  read_tsv(PATHS[["current_gene_ko_procrustes"]], show_col_types = FALSE, na = "NA"),
  read_tsv(PATHS[["current_function_taxonomy"]], show_col_types = FALSE, na = "NA")
) %>%
  filter(subset == "nahuelbuta_depth_9") %>%
  mutate(across(
    c(profile_x, profile_y),
    ~ recode(.x, "Recovered MAG taxonomy" = "Genome-size-normalised MAG assemblage")
  ))
procrustes_comparison <- exact_procrustes %>%
  left_join(
    current_procrustes %>%
      select(
        profile_x, profile_y,
        current_v39_r = procrustes_r,
        current_v39_p = p
      ),
    by = c("profile_x", "profile_y")
  ) %>%
  mutate(
    r_difference = procrustes_r - current_v39_r,
    p_difference = p - current_v39_p,
    scientific_conclusion_changed = classification(p) != classification(current_v39_p)
  )
write_tsv_checked(
  procrustes_comparison,
  file.path(STAT_DIR, "exact_depth_procrustes_vs_v39.tsv")
)

# ---- Compact Phase 1 numerical summary -------------------------------------
phase1_numeric_summary <- bind_rows(
  permutation_audit %>%
    transmute(
      domain = "permutation",
      item = "three-pit exact space",
      result = paste0(
        theoretical_arrangements_including_observed, " total; ",
        theoretical_nonidentity_permutations, " nonidentity; requested-999 generated ",
        generated_rows_from_how_shuffleSet, " unique nonidentity rows"
      )
    ),
  tpm_feature_audit %>%
    filter(metric %in% c(
      "catalog_gene_clusters",
      "all_zero_or_unmapped_features_across_17_samples",
      "composition_matrix_features_after_all_zero_filter"
    )) %>%
    transmute(domain = "feature_count", item = metric, result = as.character(value)),
  ordination_comparison %>%
    transmute(
      domain = "missingness_sensitivity",
      item = scenario,
      result = sprintf(
        "two-axis Procrustes r=%.6f; pairwise-distance Spearman r=%.6f",
        two_axis_procrustes_r,
        spearman_pairwise_distance_r
      )
    )
)
write_tsv_checked(phase1_numeric_summary, file.path(SUPPORT, "V40_PHASE1_NUMERIC_SUMMARY.tsv"))

capture.output(sessionInfo(), file = file.path(SUPPORT, "R_sessionInfo.txt"))
message("PASS: V40 Phase 1 exact-permutation, sensitivity, and TPM audit completed.")
