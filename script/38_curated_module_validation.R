# =============================================================================
# Curated functional-module validation for targeted scientific restructuring.
#
# Marker identities were frozen before inspecting group-level results. The
# three primary modules are nitrogen transformations, pedogenic resource
# acquisition, and alternative-energy chemolithotrophy. Stress response and
# CAZy substrate-processing potential are supplementary modules.
#
# Samples are independent biological replicates. There are no technical
# replicates and no cross-depth pairing is imposed. All outputs are tables or
# text reports; this script does not modify manuscripts or render figures.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(vegan)
})

try(Sys.setlocale("LC_ALL", "en_US.UTF-8"), silent = TRUE)
set.seed(42)

PROJECT <- normalizePath(Sys.getenv("CHILE_PROJECT_ROOT", getwd()), mustWork = TRUE)
GENE_DIR <- file.path(PROJECT, "data/processed/metagenome_output/Genecatalog")
META_FILE <- file.path(PROJECT, "data/metadata/env_2024_meta.tsv")
OUT <- file.path(
  PROJECT,
  "analysis/targeted_scientific_restructuring_2026-07-13/curated_module_validation"
)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(PROJECT, "script/_functional_marker_registry.R"))

samples <- c(
  "C01", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10",
  "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18"
)

module_order <- c(
  "nitrogen_transformations",
  "pedogenic_resource_acquisition",
  "alternative_energy_chemolithotrophy",
  "stress_response",
  "cazy_substrate_processing"
)

meta <- read_tsv(META_FILE, show_col_types = FALSE) %>%
  filter(sample_id %in% samples) %>%
  arrange(match(sample_id, samples)) %>%
  mutate(
    site_display = factor(
      recode(site, "NB" = "NA"),
      levels = c("AZ", "SG", "LC", "NA")
    ),
    depth_display = factor(
      recode(
        depth_label,
        "5_10cm" = "surface",
        "20_40cm" = "middle",
        "60_80cm" = "deep"
      ),
      levels = c("surface", "middle", "deep")
    )
  )

stopifnot(
  identical(meta$sample_id, samples),
  sum(meta$depth_label == "5_10cm") == 11,
  sum(meta$site == "NB") == 9,
  identical(as.integer(table(meta$site[meta$depth_label == "5_10cm"])), c(2L, 3L, 3L, 3L)),
  identical(as.integer(table(meta$depth_label[meta$site == "NB"])), c(3L, 3L, 3L))
)

write_tsv(
  functional_marker_registry %>% arrange(match(module, module_order), category, KO),
  file.path(OUT, "curated_marker_registry.tsv")
)
write_tsv(legacy_marker_corrections, file.path(OUT, "legacy_marker_corrections.tsv"))

# ---- Annotation parsing and biological filters -----------------------------
message("Reading selected eggNOG annotation columns...")
annotation <- fread(
  file.path(GENE_DIR, "eggNog.tsv"),
  select = c("GeneName", "Description", "Name", "KO", "CAZy"),
  na.strings = c("", "NA")
)

target_kos <- unique(functional_marker_registry$KO)
ko_pattern <- paste(target_kos, collapse = "|")

ko_annotation <- annotation[
  !is.na(KO) & KO != "-" & grepl(ko_pattern, KO),
  c("GeneName", "Description", "Name", "KO")
] %>%
  as_tibble() %>%
  separate_rows(KO, sep = ",") %>%
  mutate(KO = str_remove(str_trim(KO), "^ko:")) %>%
  filter(KO %in% target_kos) %>%
  distinct(GeneName, KO, .keep_all = TRUE) %>%
  inner_join(functional_marker_registry, by = "KO") %>%
  mutate(
    filter_pass = functional_marker_filter_pass(
      gene_filter,
      Name,
      Description
    )
  )

marker_filter_audit <- functional_marker_registry %>%
  left_join(
    ko_annotation %>%
      group_by(module, category, KO) %>%
      summarise(
        candidate_genes = n_distinct(GeneName),
        retained_genes = n_distinct(GeneName[filter_pass]),
        excluded_by_name_filter = n_distinct(GeneName[!filter_pass]),
        .groups = "drop"
      ),
    by = c("module", "category", "KO")
  ) %>%
  mutate(
    across(
      c(candidate_genes, retained_genes, excluded_by_name_filter),
      ~ replace_na(as.integer(.x), 0L)
    )
  ) %>%
  arrange(match(module, module_order), category, KO)

write_tsv(marker_filter_audit, file.path(OUT, "marker_filter_audit.tsv"))

ko_mapping <- ko_annotation %>%
  filter(filter_pass) %>%
  transmute(
    GeneName,
    module,
    tier,
    category,
    marker = KO,
    marker_label = symbol,
    strict_core,
    specificity_note,
    source_basis,
    mapping_basis = str_c("KO; filter=", gene_filter)
  ) %>%
  distinct()

# CAZy genes are retained only when all targeted family assignments for a gene
# point to one substrate category. This avoids arbitrary first-family choices.
cazy_lookup <- imap_dfr(
  cazy_substrate_families,
  ~ tibble(marker = .x, category = .y)
)

cazy_candidates <- annotation[
  !is.na(CAZy) & CAZy != "-",
  c("GeneName", "CAZy")
] %>%
  as_tibble() %>%
  separate_rows(CAZy, sep = ",") %>%
  transmute(
    GeneName,
    marker = str_remove(str_trim(CAZy), "_[0-9]+$")
  ) %>%
  inner_join(cazy_lookup, by = "marker") %>%
  distinct(GeneName, marker, category)

cazy_gene_status <- cazy_candidates %>%
  group_by(GeneName) %>%
  summarise(
    targeted_families = n_distinct(marker),
    targeted_substrates = n_distinct(category),
    substrate_labels = str_c(sort(unique(category)), collapse = ";"),
    retain_conservative = targeted_substrates == 1,
    .groups = "drop"
  )

cazy_mapping <- cazy_candidates %>%
  inner_join(cazy_gene_status, by = "GeneName") %>%
  filter(retain_conservative) %>%
  transmute(
    GeneName,
    module = "cazy_substrate_processing",
    tier = "supplementary",
    category,
    marker,
    marker_label = marker,
    strict_core = TRUE,
    specificity_note = "Targeted CAZy family; cross-substrate multi-annotations excluded",
    source_basis = "Frozen project CAZy substrate-family registry",
    mapping_basis = "CAZy family; unique targeted substrate"
  ) %>%
  distinct()

cazy_mapping_audit_observed <- cazy_candidates %>%
  left_join(cazy_gene_status, by = "GeneName") %>%
  group_by(marker, category) %>%
  summarise(
    candidate_genes = n_distinct(GeneName),
    retained_genes = n_distinct(GeneName[retain_conservative]),
    excluded_cross_substrate_genes = n_distinct(GeneName[!retain_conservative]),
    .groups = "drop"
  )

cazy_mapping_audit <- cazy_lookup %>%
  left_join(cazy_mapping_audit_observed, by = c("marker", "category")) %>%
  mutate(
    across(
      c(candidate_genes, retained_genes, excluded_cross_substrate_genes),
      ~ replace_na(as.integer(.x), 0L)
    )
  ) %>%
  arrange(category, marker)

write_tsv(cazy_mapping_audit, file.path(OUT, "cazy_mapping_audit.tsv"))

feature_mapping <- bind_rows(ko_mapping, cazy_mapping) %>%
  mutate(module = factor(module, levels = module_order)) %>%
  arrange(module, category, marker, GeneName)

mapping_summary <- feature_mapping %>%
  group_by(module, tier, category) %>%
  summarise(
    detected_markers = n_distinct(marker),
    mapped_genes = n_distinct(GeneName),
    strict_mapped_genes = n_distinct(GeneName[strict_core]),
    .groups = "drop"
  ) %>%
  arrange(module, category)
write_tsv(mapping_summary, file.path(OUT, "feature_mapping_summary.tsv"))

# ---- Global gene-catalog TPM denominator -----------------------------------
message("Reading counts and lengths for the full-catalog RPK denominator...")
counts <- fread(
  file.path(GENE_DIR, "Nmapped_reads.tsv"),
  select = c("GeneName", samples)
)
lengths <- fread(
  file.path(GENE_DIR, "seq_len_nucl.tsv"),
  header = FALSE,
  col.names = c("GeneName", "Length")
)

gene_rpk <- merge(counts, lengths, by = "GeneName")
gene_rpk <- gene_rpk[is.finite(Length) & Length > 0]
for (sample_id in samples) {
  set(
    gene_rpk,
    j = sample_id,
    value = gene_rpk[[sample_id]] / (gene_rpk$Length / 1000)
  )
}
global_rpk <- colSums(as.matrix(gene_rpk[, ..samples]), na.rm = TRUE)
stopifnot(all(is.finite(global_rpk)), all(global_rpk > 0))
gene_rpk[, Length := NULL]
rm(counts, lengths, annotation, ko_annotation, cazy_candidates)
invisible(gc())

mapped_genes <- unique(as.character(feature_mapping$GeneName))
mapped_rpk <- gene_rpk[GeneName %chin% mapped_genes]
rm(gene_rpk)
invisible(gc())

mapping_dt <- as.data.table(feature_mapping)
mapping_dt[, module := as.character(module)]
mapped <- merge(mapping_dt, mapped_rpk, by = "GeneName", allow.cartesian = TRUE)
stopifnot(nrow(mapped) > 0)
rm(mapped_rpk, mapping_dt)
invisible(gc())

make_feature_tpm <- function(resolution, assignment, strict_only, variant_name) {
  x <- copy(mapped)
  if (strict_only) x <- x[strict_core == TRUE]
  if (nrow(x) == 0) return(tibble())

  if (assignment == "fractional") {
    x[, assignment_weight := 1 / .N, by = .(GeneName, module)]
  } else {
    x[, assignment_weight := 1]
  }

  for (sample_id in samples) {
    set(x, j = sample_id, value = x[[sample_id]] * x$assignment_weight)
  }

  if (resolution == "marker") {
    by_cols <- c("module", "tier", "category", "marker", "marker_label")
  } else {
    x[, `:=`(marker = category, marker_label = category)]
    by_cols <- c("module", "tier", "category", "marker", "marker_label")
  }

  result <- x[
    , lapply(.SD, sum, na.rm = TRUE),
    by = by_cols,
    .SDcols = samples
  ]
  for (sample_id in samples) {
    set(
      result,
      j = sample_id,
      value = result[[sample_id]] / (global_rpk[[sample_id]] / 1e6)
    )
  }

  as_tibble(result) %>%
    pivot_longer(all_of(samples), names_to = "sample_id", values_to = "TPM") %>%
    left_join(
      meta %>% select(sample_id, site, site_display, depth_label, depth_display),
      by = "sample_id"
    ) %>%
    mutate(
      analysis_variant = variant_name,
      resolution = resolution,
      assignment = assignment,
      strict_only = strict_only
    )
}

variant_spec <- tribble(
  ~analysis_variant, ~resolution, ~assignment, ~strict_only, ~predefined_primary,
  "primary_marker_fractional", "marker", "fractional", FALSE, TRUE,
  "sensitivity_category_fractional", "category", "fractional", FALSE, FALSE,
  "sensitivity_marker_full", "marker", "full", FALSE, FALSE,
  "sensitivity_marker_fractional_strict", "marker", "fractional", TRUE, FALSE
)

message("Aggregating curated marker and category TPM tables...")
feature_tpm <- pmap_dfr(
  variant_spec,
  \(analysis_variant, resolution, assignment, strict_only, predefined_primary) {
    make_feature_tpm(resolution, assignment, strict_only, analysis_variant)
  }
) %>%
  mutate(
    module = factor(module, levels = module_order),
    predefined_primary = analysis_variant == "primary_marker_fractional"
  ) %>%
  arrange(analysis_variant, module, category, marker, sample_id)

write_tsv(feature_tpm, file.path(OUT, "curated_feature_sample_TPM.tsv"))

# ---- Uniform multivariate tests --------------------------------------------
expected_feature_count <- function(variant_name, module_name) {
  spec <- variant_spec %>% filter(analysis_variant == variant_name)
  stopifnot(nrow(spec) == 1)

  if (module_name == "cazy_substrate_processing") {
    registry <- cazy_lookup
  } else {
    registry <- functional_marker_registry %>%
      filter(.data$module == module_name)
    if (spec$strict_only) registry <- registry %>% filter(strict_core)
    registry <- registry %>% rename(marker = KO)
  }

  if (spec$resolution == "marker") {
    n_distinct(registry$marker)
  } else {
    n_distinct(registry$category)
  }
}

empty_multivariate_result <- function(
    variant_name, module_name, tier_name, subset_name, subset_meta,
    n_features_registered, status) {
  tibble(
    analysis_variant = variant_name,
    predefined_primary = variant_name == "primary_marker_fractional",
    module = module_name,
    tier = tier_name,
    subset = subset_name,
    n_samples = nrow(subset_meta),
    group_sizes = paste(
      names(table(subset_meta$analysis_group)),
      as.integer(table(subset_meta$analysis_group)),
      sep = ":",
      collapse = ";"
    ),
    n_features_registered = n_features_registered,
    n_features_detected = NA_integer_,
    distance = "Bray-Curtis on Hellinger-transformed TPM",
    PERMANOVA_Df = NA_real_,
    PERMANOVA_R2 = NA_real_,
    PERMANOVA_F = NA_real_,
    PERMANOVA_p = NA_real_,
    dispersion_F = NA_real_,
    dispersion_p = NA_real_,
    dispersion_distance_correction = "Lingoes",
    status = status
  )
}

run_multivariate_test <- function(variant_name, module_name, subset_name) {
  subset_meta <- if (subset_name == "surface_climate") {
    meta %>%
      filter(depth_label == "5_10cm") %>%
      mutate(analysis_group = droplevels(site_display))
  } else {
    meta %>%
      filter(site == "NB") %>%
      mutate(analysis_group = droplevels(depth_display))
  }
  subset_meta <- subset_meta %>% arrange(match(sample_id, samples))

  dat <- feature_tpm %>%
    filter(
      analysis_variant == variant_name,
      as.character(module) == module_name,
      sample_id %in% subset_meta$sample_id
    )
  tier_name <- unique(dat$tier)
  if (length(tier_name) == 0) {
    tier_name <- if_else(module_name %in% module_order[1:3], "primary", "supplementary")
  }
  n_registered <- expected_feature_count(variant_name, module_name)
  if (nrow(dat) == 0) {
    return(empty_multivariate_result(
      variant_name, module_name, tier_name, subset_name, subset_meta,
      n_registered, "no mapped features"
    ))
  }

  comm <- dat %>%
    select(sample_id, marker, TPM) %>%
    group_by(sample_id, marker) %>%
    summarise(TPM = sum(TPM), .groups = "drop") %>%
    pivot_wider(names_from = marker, values_from = TPM, values_fill = 0) %>%
    arrange(match(sample_id, subset_meta$sample_id))
  stopifnot(identical(comm$sample_id, subset_meta$sample_id))

  comm_mat <- as.matrix(select(comm, -sample_id))
  comm_mat <- comm_mat[, colSums(comm_mat) > 0, drop = FALSE]
  if (ncol(comm_mat) < 2) {
    result <- empty_multivariate_result(
      variant_name, module_name, tier_name, subset_name, subset_meta,
      n_registered, "fewer than two detected features"
    )
    result$n_features_detected <- ncol(comm_mat)
    return(result)
  }
  if (any(rowSums(comm_mat) <= 0)) {
    result <- empty_multivariate_result(
      variant_name, module_name, tier_name, subset_name, subset_meta,
      n_registered, "one or more samples have zero module abundance"
    )
    result$n_features_detected <- ncol(comm_mat)
    return(result)
  }

  comm_hellinger <- decostand(comm_mat, method = "hellinger")
  distance <- vegdist(comm_hellinger, method = "bray")
  if (any(!is.finite(as.vector(distance))) || max(distance) == 0) {
    result <- empty_multivariate_result(
      variant_name, module_name, tier_name, subset_name, subset_meta,
      n_registered, "non-finite or invariant distance matrix"
    )
    result$n_features_detected <- ncol(comm_mat)
    return(result)
  }

  tryCatch({
    set.seed(42)
    ad <- adonis2(
      distance ~ analysis_group,
      data = subset_meta,
      permutations = 999
    )
    dispersion <- betadisper(
      distance,
      subset_meta$analysis_group,
      type = "median",
      add = "lingoes"
    )
    set.seed(42)
    dispersion_perm <- permutest(dispersion, permutations = 999)

    tibble(
      analysis_variant = variant_name,
      predefined_primary = variant_name == "primary_marker_fractional",
      module = module_name,
      tier = tier_name,
      subset = subset_name,
      n_samples = nrow(subset_meta),
      group_sizes = paste(
        names(table(subset_meta$analysis_group)),
        as.integer(table(subset_meta$analysis_group)),
        sep = ":",
        collapse = ";"
      ),
      n_features_registered = n_registered,
      n_features_detected = ncol(comm_mat),
      distance = "Bray-Curtis on Hellinger-transformed TPM",
      PERMANOVA_Df = ad$Df[1],
      PERMANOVA_R2 = ad$R2[1],
      PERMANOVA_F = ad$F[1],
      PERMANOVA_p = ad$`Pr(>F)`[1],
      dispersion_F = dispersion_perm$tab[1, "F"],
      dispersion_p = dispersion_perm$tab[1, "Pr(>F)"],
      dispersion_distance_correction = "Lingoes",
      status = "complete"
    )
  }, error = function(e) {
    result <- empty_multivariate_result(
      variant_name, module_name, tier_name, subset_name, subset_meta,
      n_registered, str_c("test failed: ", conditionMessage(e))
    )
    result$n_features_detected <- ncol(comm_mat)
    result
  })
}

message("Running uniform PERMANOVA and dispersion tests...")
multivariate_tests <- crossing(
  analysis_variant = variant_spec$analysis_variant,
  module = module_order,
  subset = c("surface_climate", "NA_depth")
) %>%
  pmap_dfr(\(analysis_variant, module, subset) {
    run_multivariate_test(analysis_variant, module, subset)
  }) %>%
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
    interpretation = case_when(
      status != "complete" ~ status,
      PERMANOVA_FDR_within_tier < 0.05 & dispersion_FDR_within_tier < 0.05 ~
        "FDR-supported location signal with an FDR-supported dispersion caveat",
      PERMANOVA_FDR_within_tier < 0.05 & dispersion_p < 0.05 ~
        "FDR-supported location signal with a nominal dispersion caveat",
      PERMANOVA_FDR_within_tier < 0.05 ~
        "FDR-supported location signal; dispersion test not significant",
      PERMANOVA_p < 0.05 & dispersion_p < 0.05 ~
        "nominal location signal with a dispersion caveat",
      PERMANOVA_p < 0.05 ~
        "nominal location signal only; dispersion test not significant",
      dispersion_FDR_within_tier < 0.05 ~
        "no location signal; FDR-supported multivariate dispersion difference",
      dispersion_p < 0.05 ~
        "no location signal; nominal multivariate dispersion difference only",
      TRUE ~ "no nominal multivariate location or dispersion signal"
    ),
    module = factor(module, levels = module_order)
  ) %>%
  arrange(analysis_variant, subset, module)

write_tsv(
  multivariate_tests,
  file.path(OUT, "curated_module_PERMANOVA_dispersion.tsv")
)

# ---- Module-total and category-level univariate tests ----------------------
kw_effect <- function(values, groups) {
  keep <- is.finite(values) & !is.na(groups)
  values <- values[keep]
  groups <- droplevels(factor(groups[keep]))
  n <- length(values)
  k <- nlevels(groups)
  if (n < 4 || k < 2 || any(table(groups) < 2)) {
    return(tibble(
      H = NA_real_, p_value = NA_real_, epsilon2 = NA_real_, n = n, k = k,
      status = "insufficient replication"
    ))
  }
  test <- kruskal.test(values ~ groups)
  H <- unname(test$statistic)
  tibble(
    H = H,
    p_value = test$p.value,
    epsilon2 = max(min((H - k + 1) / (n - k), 1), 0),
    n = n,
    k = k,
    status = "complete"
  )
}

test_long_values <- function(dat, value_col, grouping_cols) {
  dat %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(
      surface_climate = list(
        kw_effect(
          .data[[value_col]][depth_label == "5_10cm"],
          site_display[depth_label == "5_10cm"]
        )
      ),
      NA_depth = list(
        kw_effect(
          .data[[value_col]][site == "NB"],
          depth_display[site == "NB"]
        )
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(
      c(surface_climate, NA_depth),
      names_to = "subset",
      values_to = "test"
    ) %>%
    unnest(test)
}

primary_marker_tpm <- feature_tpm %>%
  filter(analysis_variant == "primary_marker_fractional")

module_totals <- primary_marker_tpm %>%
  group_by(module, tier, sample_id, site, site_display, depth_label, depth_display) %>%
  summarise(module_total_TPM = sum(TPM), .groups = "drop")

module_total_tests <- test_long_values(
  module_totals,
  "module_total_TPM",
  c("module", "tier")
) %>%
  group_by(tier) %>%
  mutate(FDR_within_tier = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  arrange(subset, factor(module, levels = module_order))

category_scores <- feature_tpm %>%
  filter(analysis_variant == "sensitivity_category_fractional") %>%
  group_by(
    module, tier, category, sample_id, site, site_display,
    depth_label, depth_display
  ) %>%
  summarise(category_TPM = sum(TPM), .groups = "drop")

category_tests <- test_long_values(
  category_scores,
  "category_TPM",
  c("module", "tier", "category")
) %>%
  group_by(module, subset) %>%
  mutate(FDR_within_module_subset = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  group_by(tier, subset) %>%
  mutate(FDR_within_tier_subset = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  group_by(tier) %>%
  mutate(FDR_within_tier_all_subsets = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(FDR_global_all_categories = p.adjust(p_value, method = "BH")) %>%
  arrange(subset, factor(module, levels = module_order), category)

write_tsv(module_totals, file.path(OUT, "curated_module_total_sample_TPM.tsv"))
write_tsv(module_total_tests, file.path(OUT, "curated_module_total_omnibus_stats.tsv"))
write_tsv(category_scores, file.path(OUT, "curated_category_sample_TPM.tsv"))
write_tsv(category_tests, file.path(OUT, "curated_category_omnibus_stats.tsv"))

# ---- Machine-readable and human-readable audit summaries ------------------
primary_results <- multivariate_tests %>%
  filter(analysis_variant == "primary_marker_fractional") %>%
  mutate(module = as.character(module))

result_lines <- primary_results %>%
  transmute(
    line = sprintf(
      "- %s, %s: R2 = %s, PERMANOVA p = %s, tier-FDR q = %s; dispersion p = %s (%s).",
      module,
      subset,
      if_else(is.na(PERMANOVA_R2), "NA", formatC(PERMANOVA_R2, digits = 3, format = "f")),
      if_else(is.na(PERMANOVA_p), "NA", formatC(PERMANOVA_p, digits = 3, format = "f")),
      if_else(
        is.na(PERMANOVA_FDR_within_tier),
        "NA",
        formatC(PERMANOVA_FDR_within_tier, digits = 3, format = "f")
      ),
      if_else(is.na(dispersion_p), "NA", formatC(dispersion_p, digits = 3, format = "f")),
      interpretation
    )
  ) %>%
  pull(line)

writeLines(
  c(
    "# Curated Functional-Module Validation",
    "",
    "## Frozen design",
    "",
    "- Marker identities were curated and frozen before group-level tests were run.",
    "- Primary modules: nitrogen transformations; pedogenic resource acquisition; alternative-energy chemolithotrophy.",
    "- Supplementary modules: stress response; targeted CAZy substrate-processing potential.",
    "- Primary multivariate endpoint: marker-level TPM with fractional assignment of multi-marker genes.",
    "- Distance: Bray-Curtis after Hellinger transformation; 999 unrestricted permutations.",
    "- PERMANOVA is paired with median-based betadisper permutation tests using a Lingoes correction.",
    "- BH-FDR is calculated separately for the six primary-module tests and four supplementary-module tests, and also across all ten tests.",
    "- Surface-climate comparison: 11 independent biological replicates (AZ 2; SG 3; LC 3; NA 3).",
    "- Nahuelbuta depth comparison: 9 independent biological replicates (surface 3; middle 3; deep 3).",
    "- No technical replicates or cross-depth pairing are assumed.",
    "",
    "## Evidence boundary",
    "",
    "- DNA marker abundance represents encoded functional potential, not activity, process direction, rate, or flux.",
    "- Broad hydrogenase and sulfur-redox homologs remain directionally unresolved.",
    "- Putative cox genes require explicit annotation names; a coxL active-site motif check was not possible because catalog protein sequences were unavailable.",
    "- K10944-K10946 were not used as confirmatory ammonia-oxidation markers because the KO assignments overlap methane-monooxygenase homologs; K10535 is retained as a hydroxylamine-oxidation marker.",
    "- Cross-substrate CAZy multi-annotations were excluded rather than assigned arbitrarily.",
    "",
    "## Predefined primary-variant results",
    "",
    result_lines,
    "",
    "## Interpretation rule",
    "",
    "- Figure 3 and manuscript decisions must use the frozen primary variant first; sensitivity variants may qualify robustness but cannot rescue a failed primary test.",
    "- Module-total and category-level tests are separate endpoint families and are reported in their own FDR-controlled tables.",
    "- No manuscript or figure was modified by this analysis.",
    "",
    "## Outputs",
    "",
    "- curated_marker_registry.tsv",
    "- legacy_marker_corrections.tsv",
    "- marker_filter_audit.tsv",
    "- cazy_mapping_audit.tsv",
    "- feature_mapping_summary.tsv",
    "- curated_feature_sample_TPM.tsv",
    "- curated_module_PERMANOVA_dispersion.tsv",
    "- curated_module_total_sample_TPM.tsv",
    "- curated_module_total_omnibus_stats.tsv",
    "- curated_category_sample_TPM.tsv",
    "- curated_category_omnibus_stats.tsv",
    "- R_sessionInfo.txt"
  ),
  file.path(OUT, "README_curated_module_validation.md")
)

capture.output(sessionInfo(), file = file.path(OUT, "R_sessionInfo.txt"))

message("Completed curated module validation: ", OUT)
print(primary_results %>% select(
  module, tier, subset, n_features_detected, PERMANOVA_R2,
  PERMANOVA_p, PERMANOVA_FDR_within_tier, dispersion_p, interpretation
))
