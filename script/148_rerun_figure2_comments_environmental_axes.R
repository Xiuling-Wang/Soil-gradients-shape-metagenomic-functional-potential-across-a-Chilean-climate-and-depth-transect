# =============================================================================
# Figure 2 comment-review reanalysis
#
# This is a candidate analysis. It does not overwrite the V04 manuscript or
# any existing figure. The primary candidate summarises the 12 correlated soil
# measurements with predictor-only PCA axes, then tests those fixed axes in the
# gene-catalog response. The historical five-variable/HP workflow is retained
# only as a clearly labelled diagnostic sensitivity.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(permute)
  library(rdacca.hp)
  library(patchwork)
  library(cowplot)
  library(ggrepel)
  library(ggtext)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplotify)
  library(grid)
})

options(scipen = 999)
set.seed(42)
source("script/_theme_sbb.R")

# Arial is the author-preferred manuscript font and is available locally.
SBB_FONT <- "Arial"

OUTDIR <- Sys.getenv(
  "CHILE_COMMENT_REVIEW_OUTDIR",
  "analysis/V05_comment_review_2026-08-17/figure2_reanalysis"
)
PANEL_DIR <- file.path(OUTDIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

PATHS <- c(
  reads = "data/processed/metagenome_output/Genecatalog/Nmapped_reads.tsv",
  lengths = "data/processed/metagenome_output/Genecatalog/seq_len_nucl.tsv",
  mag_counts = "data/processed/metagenome_output/MAGs/counts_genomes.tsv",
  mag_quality = "data/processed/metagenome_output/MAGs/genome_quality.tsv",
  metadata = "data/metadata/env_2024_meta.tsv"
)

missing_inputs <- PATHS[!file.exists(PATHS)]
if (length(missing_inputs)) {
  stop("Missing required inputs: ", paste(missing_inputs, collapse = ", "))
}

EXPECTED_SAMPLES <- c("C01", sprintf("C%02d", 3:18))
ENV_VARS <- c(
  "pH", "Conductivity", "moisture", "C", "N", "Feo",
  "Alo", "Mno", "Sio", "NH4", "NO3", "Pt"
)
CURRENT_FIVE <- c("Alo", "Conductivity", "NO3", "moisture", "Mno")

PRETTY_ENV <- c(
  pH = "pH",
  Conductivity = "Conductivity",
  moisture = "Moisture",
  C = "Total C",
  N = "Total N",
  Feo = "Fe\u2092\u2093",
  Alo = "Al\u2092\u2093",
  Mno = "Mn\u2092\u2093",
  Sio = "Si\u2092\u2093",
  NH4 = "NH\u2084\u207a",
  NO3 = "NO\u2083\u207b",
  Pt = "Total P"
)

SITE_LEVELS <- c("AZ", "SG", "LC", "NB")
SITE_LABELS <- c(AZ = "AZ", SG = "SG", LC = "LC", NB = "NA")
SITE_COLORS <- pal_site[SITE_LEVELS]
DEPTH_LEVELS <- c("5_10cm", "20_40cm", "60_80cm")
DEPTH_LABELS <- c(
  `5_10cm` = "5-10 cm",
  `20_40cm` = "20-40 cm",
  `60_80cm` = "60-80 cm"
)
DEPTH_SHAPES <- c(`5_10cm` = 21, `20_40cm` = 22, `60_80cm` = 24)

write_tsv_checked <- function(x, filename) {
  path <- file.path(OUTDIR, filename)
  readr::write_tsv(x, path, na = "NA")
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Failed to write output: ", path)
  }
  path
}

format_p <- function(p) {
  if (is.na(p)) return("not tested")
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}

message("Reading and aligning the 17 metagenomes...")
reads <- read_tsv(PATHS[["reads"]], show_col_types = FALSE)
lengths <- read_tsv(
  PATHS[["lengths"]],
  col_names = c("GeneName", "Length"),
  show_col_types = FALSE
)
metadata_raw <- read_tsv(
  PATHS[["metadata"]],
  show_col_types = FALSE,
  na = c("", "NA")
)

if (!all(EXPECTED_SAMPLES %in% names(reads))) {
  stop("The gene-count table does not contain the expected 17 samples.")
}
if ("C02" %in% EXPECTED_SAMPLES) stop("C02 must not enter the metagenomic set.")

sample_cols <- EXPECTED_SAMPLES
metadata <- metadata_raw %>%
  filter(sample_id %in% sample_cols) %>%
  arrange(match(sample_id, sample_cols)) %>%
  mutate(
    site = factor(site, levels = SITE_LEVELS),
    depth_label = factor(depth_label, levels = DEPTH_LEVELS),
    site_pit = as.character(site_pit)
  )

if (!identical(metadata$sample_id, sample_cols)) {
  stop("Metadata and gene-count sample order do not match.")
}

expected_nb_map <- tribble(
  ~site_pit, ~depth_label, ~sample_id,
  "NB1", "5_10cm", "C10",
  "NB1", "20_40cm", "C13",
  "NB1", "60_80cm", "C16",
  "NB2", "5_10cm", "C11",
  "NB2", "20_40cm", "C14",
  "NB2", "60_80cm", "C17",
  "NB3", "5_10cm", "C12",
  "NB3", "20_40cm", "C15",
  "NB3", "60_80cm", "C18"
)
observed_nb_map <- metadata %>%
  filter(site == "NB") %>%
  transmute(site_pit, depth_label = as.character(depth_label), sample_id) %>%
  arrange(site_pit, match(depth_label, DEPTH_LEVELS))
expected_nb_map <- expected_nb_map %>%
  arrange(site_pit, match(depth_label, DEPTH_LEVELS))
if (!identical(observed_nb_map, expected_nb_map)) {
  stop("Nahuelbuta pit-depth map differs from the frozen design decision.")
}

design_map <- tribble(
  ~field, ~value,
  "Experimental unit", "Distinct soil sample",
  "Biological replicate", "Three pits per represented site/depth cell, except AZ metagenomes n = 2",
  "Technical replicate", "None",
  "Response", "Full gene-catalog TPM profile, Hellinger transformed",
  "Primary candidate predictor", "Soil PC1 and PC2 from 12 z-scaled physicochemical measurements",
  "Surface test", "11 surface metagenomes; unrestricted permutations",
  "Nahuelbuta test", "Nine samples in three pit-resolved depth profiles; permutations restricted within site_pit",
  "All-17 role", "Configuration and effect-size visualisation; no global permutation p-value promoted",
  "Missingness", "Design-based, depth-matched imputation within Nahuelbuta; alternatives reported as sensitivity",
  "Historical five-variable role", "Response-guided diagnostic sensitivity; not the candidate primary model"
)
write_tsv_checked(design_map, "design_map.tsv")
write_tsv_checked(observed_nb_map, "nahuelbuta_pit_depth_map.tsv")

missingness <- map_dfr(ENV_VARS, function(variable) {
  tibble(
    variable = variable,
    n_missing = sum(is.na(metadata[[variable]])),
    missing_samples = paste(metadata$sample_id[is.na(metadata[[variable]])], collapse = ";")
  )
})
write_tsv_checked(missingness, "environment_missingness.tsv")

impute_environment <- function(data, strategy) {
  out <- data
  log_rows <- list()
  for (variable in ENV_VARS) {
    missing_index <- which(is.na(out[[variable]]))
    if (!length(missing_index)) next

    for (index in missing_index) {
      if (strategy == "depth_matched_mean") {
        if (as.character(out$site[index]) != "NB") {
          stop("Depth-matched imputation encountered non-Nahuelbuta missingness.")
        }
        donors <- out %>%
          filter(
            site == "NB",
            depth_label == out$depth_label[index],
            !is.na(.data[[variable]])
          ) %>%
          pull(all_of(variable))
        value <- mean(donors)
        donor_rule <- "Other observed Nahuelbuta pits at the same depth"
      } else if (strategy == "global_mean") {
        value <- mean(out[[variable]], na.rm = TRUE)
        donor_rule <- "Overall observed mean"
      } else if (strategy == "global_median") {
        value <- median(out[[variable]], na.rm = TRUE)
        donor_rule <- "Overall observed median"
      } else {
        stop("Unknown imputation strategy: ", strategy)
      }

      out[[variable]][index] <- value
      log_rows[[length(log_rows) + 1L]] <- tibble(
        strategy = strategy,
        sample_id = out$sample_id[index],
        site_pit = out$site_pit[index],
        depth_label = as.character(out$depth_label[index]),
        variable = variable,
        imputed_value = value,
        donor_rule = donor_rule
      )
    }
  }
  if (anyNA(out[, ENV_VARS])) stop("Environmental imputation left missing values.")
  list(data = out, log = bind_rows(log_rows))
}

primary_imputation <- impute_environment(metadata, "depth_matched_mean")
metadata_imp <- primary_imputation$data
write_tsv_checked(primary_imputation$log, "environment_imputation_log.tsv")

message("Computing full gene-catalog TPM and Hellinger profiles...")
reads_with_length <- reads %>%
  select(GeneName, all_of(sample_cols)) %>%
  inner_join(lengths, by = "GeneName") %>%
  filter(Length > 0)

rpk <- sweep(
  as.matrix(reads_with_length[, sample_cols]),
  1,
  reads_with_length$Length / 1000,
  "/"
)
tpm <- sweep(rpk, 2, colSums(rpk, na.rm = TRUE) / 1e6, "/")
gene_profile <- t(tpm)
rownames(gene_profile) <- sample_cols
gene_profile <- gene_profile[, colSums(gene_profile, na.rm = TRUE) > 0, drop = FALSE]
gene_hell <- decostand(gene_profile, method = "hellinger")
gene_dist <- vegdist(gene_hell, method = "bray")

message("Building predictor-only soil axes from all 12 measurements...")
soil_pca <- prcomp(
  metadata_imp[, ENV_VARS],
  center = TRUE,
  scale. = TRUE
)
soil_variance <- soil_pca$sdev^2 / sum(soil_pca$sdev^2)
soil_scores <- as.data.frame(soil_pca$x[, 1:4, drop = FALSE])
colnames(soil_scores) <- paste0("SoilPC", 1:4)
rownames(soil_scores) <- metadata_imp$sample_id

pca_variance_tbl <- tibble(
  axis = paste0("SoilPC", seq_along(soil_variance)),
  variance_fraction = soil_variance,
  cumulative_fraction = cumsum(soil_variance)
)
pca_loadings_tbl <- as_tibble(
  soil_pca$rotation,
  rownames = "variable"
) %>%
  mutate(variable_label = unname(PRETTY_ENV[variable]), .after = variable)
pca_scores_tbl <- bind_cols(
  tibble(sample_id = metadata_imp$sample_id),
  as_tibble(soil_pca$x)
) %>%
  left_join(
    metadata_imp %>% select(sample_id, site, site_pit, depth_label),
    by = "sample_id"
  )
write_tsv_checked(pca_variance_tbl, "soil_pca_variance.tsv")
write_tsv_checked(pca_loadings_tbl, "soil_pca_loadings.tsv")
write_tsv_checked(pca_scores_tbl, "soil_pca_scores.tsv")

fit_soil_axis_model <- function(indices, label, permutation = NULL) {
  response <- gene_hell[indices, , drop = FALSE]
  predictors <- soil_scores[indices, c("SoilPC1", "SoilPC2"), drop = FALSE]
  formula <- response ~ SoilPC1 + SoilPC2
  model <- dbrda(formula, data = predictors, distance = "bray")
  r2 <- RsquareAdj(model)
  constrained_eigen <- eigenvals(model, model = "constrained")
  axis_fraction <- constrained_eigen / sum(constrained_eigen)

  p_value <- NA_real_
  if (!is.null(permutation)) {
    set.seed(42)
    p_value <- anova.cca(model, permutations = permutation)$`Pr(>F)`[1]
  }

  list(
    model = model,
    table = tibble(
      subset = label,
      n_samples = length(indices),
      n_predictors = 2L,
      residual_df = length(indices) - 3L,
      R2 = r2$r.squared,
      adjusted_R2 = r2$adj.r.squared,
      p = p_value,
      dbRDA1_fraction = axis_fraction[1],
      dbRDA2_fraction = axis_fraction[2],
      max_VIF = max(vif.cca(model)),
      permutation_scheme = case_when(
        is.null(permutation) ~ "not tested for the mixed all-17 configuration",
        is.numeric(permutation) ~ "999 unrestricted permutations",
        TRUE ~ "complete/Monte Carlo permutations restricted within site_pit"
      )
    )
  )
}

index_all <- seq_len(nrow(metadata_imp))
index_surface <- which(metadata_imp$depth_label == "5_10cm")
index_depth <- which(metadata_imp$site == "NB")
depth_permutation <- how(
  nperm = 999,
  blocks = factor(metadata_imp$site_pit[index_depth])
)

soil_model_all <- fit_soil_axis_model(index_all, "all_17", NULL)
soil_model_surface <- fit_soil_axis_model(index_surface, "surface_11", 999)
soil_model_depth <- fit_soil_axis_model(
  index_depth,
  "nahuelbuta_depth_9",
  depth_permutation
)
soil_model_tbl <- bind_rows(
  soil_model_all$table,
  soil_model_surface$table,
  soil_model_depth$table
)
write_tsv_checked(soil_model_tbl, "soil_axis_gene_catalog_models.tsv")

message("Running function-taxonomy concordance with design-aware subset tests...")
mag_counts <- read_tsv(PATHS[["mag_counts"]], show_col_types = FALSE) %>%
  rename(MAG = Sample)
mag_quality <- read_tsv(PATHS[["mag_quality"]], show_col_types = FALSE) %>%
  transmute(MAG = `Bin Id`, Genome_Size)
mag_with_size <- mag_counts %>% left_join(mag_quality, by = "MAG")
if (any(is.na(mag_with_size$Genome_Size)) || any(mag_with_size$Genome_Size <= 0)) {
  stop("Missing or invalid MAG genome sizes.")
}
mag_rpk <- sweep(
  as.matrix(mag_with_size[, sample_cols]),
  1,
  mag_with_size$Genome_Size / 1000,
  "/"
)
mag_relative <- sweep(mag_rpk, 2, colSums(mag_rpk, na.rm = TRUE), "/")
mag_hell <- decostand(t(mag_relative), method = "hellinger")
rownames(mag_hell) <- sample_cols

run_procrustes <- function(indices, label, permutation = NULL, keep_coordinates = FALSE) {
  gene_subset_dist <- vegdist(gene_hell[indices, , drop = FALSE], method = "bray")
  mag_subset_dist <- vegdist(mag_hell[indices, , drop = FALSE], method = "bray")
  gene_pcoa <- wcmdscale(gene_subset_dist, k = 2, eig = TRUE)
  mag_pcoa <- wcmdscale(mag_subset_dist, k = 2, eig = TRUE)
  gene_xy <- gene_pcoa$points[, 1:2, drop = FALSE]
  mag_xy <- mag_pcoa$points[, 1:2, drop = FALSE]

  proc <- procrustes(mag_xy, gene_xy, symmetric = TRUE)
  statistic <- sqrt(1 - proc$ss)
  p_value <- NA_real_
  if (!is.null(permutation)) {
    set.seed(42)
    proc_test <- protest(
      mag_xy,
      gene_xy,
      permutations = permutation,
      symmetric = TRUE
    )
    statistic <- proc_test$t0
    p_value <- proc_test$signif
  }

  coordinates <- NULL
  if (keep_coordinates) {
    coordinates <- tibble(
      sample_id = rownames(proc$X),
      mag_axis1 = proc$X[, 1],
      mag_axis2 = proc$X[, 2],
      gene_axis1 = proc$Yrot[, 1],
      gene_axis2 = proc$Yrot[, 2]
    ) %>%
      left_join(
        metadata_imp %>% select(sample_id, site, site_pit, depth_label),
        by = "sample_id"
      )
  }

  list(
    table = tibble(
      subset = label,
      n_samples = length(indices),
      procrustes_r = statistic,
      p = p_value,
      permutation_scheme = case_when(
        is.null(permutation) ~ "effect size only",
        is.numeric(permutation) ~ "999 unrestricted permutations",
        TRUE ~ "complete/Monte Carlo permutations restricted within site_pit"
      )
    ),
    coordinates = coordinates
  )
}

proc_all <- run_procrustes(index_all, "all_17", NULL, keep_coordinates = TRUE)
proc_surface <- run_procrustes(index_surface, "surface_11", 999)
proc_depth <- run_procrustes(
  index_depth,
  "nahuelbuta_depth_9",
  depth_permutation
)
procrustes_tbl <- bind_rows(
  proc_all$table,
  proc_surface$table,
  proc_depth$table
)
write_tsv_checked(procrustes_tbl, "function_taxonomy_procrustes.tsv")
write_tsv_checked(proc_all$coordinates, "function_taxonomy_procrustes_coordinates.tsv")

fit_taxonomy_conditioned <- function(indices, label, permutation = NULL) {
  gene_subset_dist <- vegdist(gene_hell[indices, , drop = FALSE], method = "bray")
  mag_subset_dist <- vegdist(mag_hell[indices, , drop = FALSE], method = "bray")
  tax_pcoa <- wcmdscale(mag_subset_dist, k = 2, eig = TRUE)
  tax_axes <- as.data.frame(tax_pcoa$points[, 1:2, drop = FALSE])
  colnames(tax_axes) <- c("TaxPCoA1", "TaxPCoA2")
  model_data <- cbind(
    soil_scores[indices, c("SoilPC1", "SoilPC2"), drop = FALSE],
    tax_axes
  )
  model <- dbrda(
    gene_subset_dist ~ SoilPC1 + SoilPC2 + Condition(TaxPCoA1 + TaxPCoA2),
    data = model_data
  )
  r2 <- RsquareAdj(model)
  p_value <- NA_real_
  if (!is.null(permutation)) {
    set.seed(42)
    p_value <- anova.cca(model, permutations = permutation)$`Pr(>F)`[1]
  }
  positive_eigen <- tax_pcoa$eig[tax_pcoa$eig > 0]
  tibble(
    subset = label,
    n_samples = length(indices),
    taxonomic_axes_conditioned = 2L,
    taxonomic_variance_fraction_conditioned = sum(positive_eigen[1:2]) / sum(positive_eigen),
    unique_environment_R2 = r2$r.squared,
    unique_environment_adjusted_R2 = r2$adj.r.squared,
    p = p_value
  )
}

taxonomy_conditioned_tbl <- bind_rows(
  fit_taxonomy_conditioned(index_all, "all_17", NULL),
  fit_taxonomy_conditioned(index_surface, "surface_11", 999),
  fit_taxonomy_conditioned(index_depth, "nahuelbuta_depth_9", depth_permutation)
)
write_tsv_checked(
  taxonomy_conditioned_tbl,
  "soil_axes_conditioned_on_mag_taxonomy.tsv"
)

message("Running imputation sensitivity for the predictor-only soil axes...")
imputation_sensitivity <- map_dfr(
  c("depth_matched_mean", "global_mean", "global_median"),
  function(strategy) {
    imputed <- impute_environment(metadata, strategy)$data
    pca <- prcomp(imputed[, ENV_VARS], center = TRUE, scale. = TRUE)
    scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
    colnames(scores) <- c("SoilPC1", "SoilPC2")
    rownames(scores) <- imputed$sample_id

    run_one <- function(indices, subset, permutation = NULL) {
      response <- gene_hell[indices, , drop = FALSE]
      model <- dbrda(
        response ~ SoilPC1 + SoilPC2,
        data = scores[indices, , drop = FALSE],
        distance = "bray"
      )
      r2 <- RsquareAdj(model)
      p_value <- NA_real_
      if (!is.null(permutation)) {
        set.seed(42)
        p_value <- anova.cca(model, permutations = permutation)$`Pr(>F)`[1]
      }
      tibble(
        strategy = strategy,
        subset = subset,
        R2 = r2$r.squared,
        adjusted_R2 = r2$adj.r.squared,
        p = p_value,
        PC1_variance = pca$sdev[1]^2 / sum(pca$sdev^2),
        PC2_variance = pca$sdev[2]^2 / sum(pca$sdev^2)
      )
    }

    bind_rows(
      run_one(index_all, "all_17", NULL),
      run_one(index_surface, "surface_11", 999),
      run_one(index_depth, "nahuelbuta_depth_9", depth_permutation)
    )
  }
)
write_tsv_checked(imputation_sensitivity, "soil_axis_imputation_sensitivity.tsv")

message("Recomputing historical-model diagnostics with design-based imputation...")
calc_vif <- function(data) {
  map_dbl(names(data), function(variable) {
    others <- setdiff(names(data), variable)
    fit <- lm(reformulate(others, response = variable), data = data)
    1 / (1 - summary(fit)$r.squared)
  })
}

vif_full <- tibble(
  model = "Full 12-variable set",
  variable = ENV_VARS,
  VIF = calc_vif(as.data.frame(metadata_imp[, ENV_VARS]))
)
vif_five <- tibble(
  model = "Historical five-variable set",
  variable = CURRENT_FIVE,
  VIF = calc_vif(as.data.frame(metadata_imp[, CURRENT_FIVE]))
)
vif_tbl <- bind_rows(vif_full, vif_five) %>%
  mutate(
    variable_label = unname(PRETTY_ENV[variable]),
    flag = case_when(VIF > 10 ~ "> 10", VIF > 5 ~ "5-10", TRUE ~ "< 5")
  )
write_tsv_checked(vif_tbl, "predictor_vif_diagnostics.tsv")

fit_effect_size_model <- function(vars, label) {
  predictors <- as.data.frame(scale(metadata_imp[, vars, drop = FALSE]))
  response <- gene_hell
  model <- dbrda(response ~ ., data = predictors, distance = "bray")
  r2 <- RsquareAdj(model)
  max_vif <- max(vif.cca(model))
  tibble(
    model = label,
    n_predictors = length(vars),
    R2 = r2$r.squared,
    adjusted_R2 = r2$adj.r.squared,
    max_VIF = max_vif
  )
}

model_comparison_tbl <- bind_rows(
  fit_effect_size_model(ENV_VARS, "Full 12-variable set"),
  fit_effect_size_model(CURRENT_FIVE, "Historical five-variable set"),
  soil_model_all$table %>%
    transmute(
      model = "Predictor-only soil PC1-PC2",
      n_predictors,
      R2,
      adjusted_R2,
      max_VIF
    )
)
write_tsv_checked(model_comparison_tbl, "model_effect_size_comparison.tsv")

historical_env <- as.data.frame(scale(metadata_imp[, CURRENT_FIVE, drop = FALSE]))
rownames(historical_env) <- metadata_imp$sample_id
set.seed(42)
historical_hp <- rdacca.hp(
  gene_dist,
  historical_env,
  method = "dbRDA",
  type = "adjR2"
)

extract_hp_percent <- function(hp_result) {
  table <- as.data.frame(hp_result$Hier.part)
  column <- intersect(c("I.perc(%)", "I.perc"), names(table))
  if (length(column) != 1L) stop("Could not identify rdacca.hp percentage column.")
  values <- table[[column]]
  names(values) <- rownames(table)
  values
}

hp_point <- extract_hp_percent(historical_hp)[CURRENT_FIVE]

# Stratified cluster bootstrap: resample the eight non-Nahuelbuta singleton
# pits and the three complete Nahuelbuta pit profiles separately. This keeps
# the original 8 + (3 x 3) mixed-design composition in every bootstrap sample.
singleton_clusters <- metadata_imp %>%
  filter(site != "NB") %>%
  pull(site_pit) %>%
  unique()
profile_clusters <- metadata_imp %>%
  filter(site == "NB") %>%
  pull(site_pit) %>%
  unique()
cluster_index <- split(seq_len(nrow(metadata_imp)), metadata_imp$site_pit)
gene_dist_matrix <- as.matrix(gene_dist)

bootstrap_target <- 999L
bootstrap_attempt <- 0L
bootstrap_success <- 0L
bootstrap_values <- matrix(
  NA_real_,
  nrow = bootstrap_target,
  ncol = length(CURRENT_FIVE),
  dimnames = list(NULL, CURRENT_FIVE)
)
set.seed(4201)
while (bootstrap_success < bootstrap_target && bootstrap_attempt < 5000L) {
  bootstrap_attempt <- bootstrap_attempt + 1L
  sampled_singletons <- sample(
    singleton_clusters,
    length(singleton_clusters),
    replace = TRUE
  )
  sampled_profiles <- sample(
    profile_clusters,
    length(profile_clusters),
    replace = TRUE
  )
  sampled_clusters <- c(sampled_singletons, sampled_profiles)
  sampled_indices <- unlist(cluster_index[sampled_clusters], use.names = FALSE)

  boot_names <- make.unique(paste0(metadata_imp$sample_id[sampled_indices], "_b"))
  boot_distance_matrix <- gene_dist_matrix[sampled_indices, sampled_indices, drop = FALSE]
  rownames(boot_distance_matrix) <- boot_names
  colnames(boot_distance_matrix) <- boot_names
  boot_environment <- historical_env[sampled_indices, , drop = FALSE]
  rownames(boot_environment) <- boot_names

  boot_hp <- tryCatch(
    suppressWarnings(
      rdacca.hp(
        as.dist(boot_distance_matrix),
        boot_environment,
        method = "dbRDA",
        type = "adjR2"
      )
    ),
    error = function(error) NULL
  )
  if (is.null(boot_hp)) next
  values <- tryCatch(extract_hp_percent(boot_hp), error = function(error) NULL)
  if (is.null(values) || !all(CURRENT_FIVE %in% names(values))) next
  values <- values[CURRENT_FIVE]
  if (!all(is.finite(values))) next

  bootstrap_success <- bootstrap_success + 1L
  bootstrap_values[bootstrap_success, ] <- values
}
if (bootstrap_success < bootstrap_target) {
  stop(
    "Only ", bootstrap_success, " successful cluster bootstraps after ",
    bootstrap_attempt, " attempts."
  )
}

hp_bootstrap_tbl <- as_tibble(bootstrap_values) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "I_perc") %>%
  group_by(variable) %>%
  summarise(
    point_I_perc = unname(hp_point[variable[1]]),
    bootstrap_mean_I_perc = mean(I_perc),
    ci_low_95 = quantile(I_perc, 0.025, names = FALSE),
    ci_high_95 = quantile(I_perc, 0.975, names = FALSE),
    n_success = n(),
    n_attempts = bootstrap_attempt,
    .groups = "drop"
  ) %>%
  mutate(variable_label = unname(PRETTY_ENV[variable])) %>%
  arrange(match(variable, CURRENT_FIVE))
write_tsv_checked(hp_bootstrap_tbl, "historical_five_variable_cluster_bootstrap.tsv")

message("Building Figure 2 candidate...")
site_scores <- scores(
  soil_model_all$model,
  display = "sites",
  scaling = 2
) %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(
    metadata_imp %>% select(sample_id, site, site_pit, depth_label),
    by = "sample_id"
  ) %>%
  mutate(
    site = factor(site, levels = SITE_LEVELS),
    depth_label = factor(depth_label, levels = DEPTH_LEVELS)
  )

bp_scores <- scores(
  soil_model_all$model,
  display = "bp",
  scaling = 2
) %>%
  as_tibble(rownames = "axis") %>%
  rename(xend = dbRDA1, yend = dbRDA2)
bp_scale <- 0.78 * max(abs(c(site_scores$dbRDA1, site_scores$dbRDA2))) /
  max(sqrt(bp_scores$xend^2 + bp_scores$yend^2))
bp_scores <- bp_scores %>%
  mutate(
    xend = xend * bp_scale,
    yend = yend * bp_scale,
    label = recode(axis, SoilPC1 = "Soil PC1", SoilPC2 = "Soil PC2")
  )

# Export the exact ordination geometry so later visual refinements can reuse
# the frozen analysis without recalculating or changing the fitted model.
write_tsv_checked(site_scores, "soil_axis_dbrda_coordinates.tsv")
write_tsv_checked(bp_scores, "soil_axis_dbrda_vectors.tsv")

surface_stats <- soil_model_surface$table
depth_stats <- soil_model_depth$table
all_stats <- soil_model_all$table
panel_a_subtitle <- sprintf(
  paste0(
    "All 17: R<sup>2</sup> = %.3f, adj. R<sup>2</sup> = %.3f; ",
    "surface <i>p</i> = %.3f; NA-depth <i>p</i> = %.3f"
  ),
  all_stats$R2,
  all_stats$adjusted_R2,
  surface_stats$p,
  depth_stats$p
)

panel_a <- ggplot(site_scores, aes(dbRDA1, dbRDA2)) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3, linetype = "dashed") +
  geom_path(
    data = site_scores %>%
      filter(site == "NB") %>%
      arrange(site_pit, depth_label),
    aes(group = site_pit),
    color = "grey68",
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  geom_segment(
    data = bp_scores,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    inherit.aes = FALSE,
    color = "#303030",
    linewidth = 0.65,
    arrow = arrow(length = unit(0.17, "cm"), type = "closed")
  ) +
  geom_text_repel(
    data = bp_scores,
    aes(x = xend * 1.12, y = yend * 1.12, label = label),
    inherit.aes = FALSE,
    family = SBB_FONT,
    color = "black",
    size = 2.7,
    fontface = "bold",
    segment.color = "grey45",
    segment.size = 0.25,
    box.padding = 0.25,
    max.overlaps = Inf
  ) +
  geom_point(
    aes(fill = site, shape = depth_label),
    size = 2.35,
    stroke = 0.45,
    color = "grey15"
  ) +
  scale_fill_manual(
    name = "Site",
    values = SITE_COLORS,
    breaks = SITE_LEVELS,
    labels = SITE_LABELS
  ) +
  scale_shape_manual(
    name = "Depth",
    values = DEPTH_SHAPES,
    labels = DEPTH_LABELS
  ) +
  labs(
    title = "Whole-gene catalog across orthogonal soil axes",
    subtitle = panel_a_subtitle,
    x = sprintf("dbRDA1 (%.1f%%)", 100 * all_stats$dbRDA1_fraction),
    y = sprintf("dbRDA2 (%.1f%%)", 100 * all_stats$dbRDA2_fraction)
  ) +
  coord_equal() +
  theme_sbb(base_size = 8) +
  theme(
    plot.tag = element_text(family = SBB_FONT, face = "bold", size = 11, color = "black"),
    plot.tag.position = "topleft",
    plot.subtitle = ggtext::element_markdown(size = 6.8, lineheight = 0.98),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = -2),
    legend.key.width = unit(0.42, "cm")
  ) +
  guides(
    fill = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 21, size = 2.6)),
    shape = guide_legend(order = 2, nrow = 1, override.aes = list(fill = "white", size = 2.6))
  ) +
  labs(tag = "A")

loading_matrix <- soil_pca$rotation[ENV_VARS, c("PC1", "PC2"), drop = FALSE]
rownames(loading_matrix) <- unname(PRETTY_ENV[ENV_VARS])
colnames(loading_matrix) <- c(
  "PC1",
  "PC2"
)
loading_limit <- max(abs(loading_matrix))
loading_col_fun <- sbb_div_col_fun(c(-loading_limit, loading_limit))
loading_heatmap <- Heatmap(
  loading_matrix,
  name = "Loading",
  col = loading_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = gpar(fontfamily = SBB_FONT, fontsize = 6.5, col = "black"),
  column_names_gp = gpar(fontfamily = SBB_FONT, fontsize = 6.6, col = "black"),
  column_names_side = "bottom",
  column_names_rot = 0,
  rect_gp = gpar(col = "grey67", lwd = 0.55),
  border = "grey45",
  width = unit(22, "mm"),
  height = unit(23, "mm"),
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "leftcenter",
    legend_width = unit(24, "mm"),
    title_gp = gpar(fontfamily = SBB_FONT, fontsize = 7.0, fontface = "bold"),
    labels_gp = gpar(fontfamily = SBB_FONT, fontsize = 6.5),
    at = c(-loading_limit, 0, loading_limit)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- loading_matrix[i, j]
    grid.text(
      sprintf("%.2f", value),
      x,
      y,
      gp = gpar(
        fontfamily = SBB_FONT,
        fontsize = 5.8,
        col = ifelse(abs(value) >= 0.31, "white", "black")
      )
    )
  }
)
loading_grob <- sbb_heatmap_grob(
  loading_heatmap,
  heatmap_legend_side = "bottom",
  padding = unit(c(3, 2, 5, 2), "mm")
)
panel_b <- cowplot::ggdraw() +
  cowplot::draw_plot(ggplotify::as.ggplot(loading_grob), x = 0, y = 0.08, width = 1, height = 0.68) +
  cowplot::draw_label(
    "B",
    x = 0.00,
    y = 0.995,
    hjust = 0,
    vjust = 1,
    fontfamily = SBB_FONT,
    fontface = "bold",
    size = 11,
    color = "black"
  ) +
  cowplot::draw_label(
    "Soil-axis loadings (12 properties)",
    x = 0.10,
    y = 0.995,
    hjust = 0,
    vjust = 1,
    fontfamily = SBB_FONT,
    fontface = "bold",
    size = 9,
    color = "black"
  ) +
  cowplot::draw_label(
    sprintf("PC1 %.1f%%; PC2 %.1f%%", 100 * soil_variance[1], 100 * soil_variance[2]),
    x = 0.53,
    y = 0.92,
    hjust = 0.5,
    vjust = 1,
    fontfamily = SBB_FONT,
    size = 6.8,
    color = "black"
  )

proc_coordinates <- proc_all$coordinates %>%
  mutate(site = factor(site, levels = SITE_LEVELS))
proc_text <- sprintf(
  paste0(
    "All 17: <i>r</i> = %.3f<br>",
    "Surface: <i>r</i> = %.3f, <i>p</i> = %.3f<br>",
    "NA depths: <i>r</i> = %.3f, <i>p</i> = %.3f"
  ),
  proc_all$table$procrustes_r,
  proc_surface$table$procrustes_r,
  proc_surface$table$p,
  proc_depth$table$procrustes_r,
  proc_depth$table$p
)

panel_c <- ggplot(proc_coordinates) +
  geom_segment(
    aes(x = mag_axis1, y = mag_axis2, xend = gene_axis1, yend = gene_axis2),
    color = "grey60",
    linewidth = 0.42
  ) +
  geom_point(
    aes(mag_axis1, mag_axis2, color = site),
    shape = 21,
    fill = "white",
    size = 2.25,
    stroke = 0.75,
    show.legend = FALSE
  ) +
  geom_point(
    aes(gene_axis1, gene_axis2, fill = site),
    shape = 21,
    color = "grey15",
    size = 2.15,
    stroke = 0.4,
    show.legend = FALSE
  ) +
  ggtext::geom_richtext(
    data = tibble(x = -Inf, y = -Inf, label = proc_text),
    aes(x = x, y = y, label = label),
    hjust = -0.04,
    vjust = -0.12,
    inherit.aes = FALSE,
    family = SBB_FONT,
    size = 2.65,
    lineheight = 1.02,
    fill = alpha("white", 0.88),
    label.color = NA,
    label.padding = unit(1.2, "mm")
  ) +
  scale_color_manual(values = SITE_COLORS) +
  scale_fill_manual(values = SITE_COLORS) +
  labs(
    title = "Function-taxonomy concordance",
    subtitle = "Hollow: MAGs; filled: gene catalog",
    x = "Procrustes axis 1",
    y = "Procrustes axis 2"
  ) +
  coord_equal(clip = "off") +
  theme_sbb(base_size = 8) +
  theme(
    plot.tag = element_text(family = SBB_FONT, face = "bold", size = 11, color = "black"),
    plot.tag.position = "topleft",
    plot.subtitle = element_text(size = 7.0)
  ) +
  labs(tag = "C")

figure_2_candidate <- panel_a + panel_b + panel_c +
  plot_layout(widths = c(1.18, 0.98, 1.08), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  file.path(OUTDIR, "Figure_2_candidate_environmental_axes.pdf"),
  figure_2_candidate,
  width = 18.3,
  height = 7.3,
  units = "cm",
  device = cairo_pdf,
  bg = "white"
)
ggsave(file.path(PANEL_DIR, "Figure_2_A.pdf"), panel_a, width = 7.3, height = 6.8, units = "cm", device = cairo_pdf, bg = "white")
ggsave(file.path(PANEL_DIR, "Figure_2_B.pdf"), panel_b, width = 5.4, height = 6.8, units = "cm", device = cairo_pdf, bg = "white")
ggsave(file.path(PANEL_DIR, "Figure_2_C.pdf"), panel_c, width = 6.8, height = 6.8, units = "cm", device = cairo_pdf, bg = "white")

message("Building supplementary model-diagnostic candidate...")
cor_matrix <- cor(
  metadata_imp[, ENV_VARS],
  method = "spearman",
  use = "pairwise.complete.obs"
)
correlation_pairs <- as.data.frame(as.table(cor_matrix), stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rename(variable_1 = Var1, variable_2 = Var2, spearman_r = Freq) %>%
  mutate(
    variable_1_label = unname(PRETTY_ENV[variable_1]),
    variable_2_label = unname(PRETTY_ENV[variable_2]),
    order_1 = match(variable_1, ENV_VARS),
    order_2 = match(variable_2, ENV_VARS)
  ) %>%
  filter(order_1 > order_2) %>%
  transmute(
    variable_1,
    variable_1_label,
    variable_2,
    variable_2_label,
    n_samples = nrow(metadata_imp),
    spearman_r,
    data_rule = "depth-matched imputation for C11/C14/C17 chemistry block"
  )
write_tsv_checked(correlation_pairs, "environment_spearman_correlations.tsv")
rownames(cor_matrix) <- unname(PRETTY_ENV[ENV_VARS])
colnames(cor_matrix) <- unname(PRETTY_ENV[ENV_VARS])
cor_display <- cor_matrix
cor_display[upper.tri(cor_display, diag = TRUE)] <- NA_real_
cor_heatmap <- Heatmap(
  cor_display,
  name = "Spearman r",
  column_title = "Environmental collinearity",
  column_title_gp = gpar(fontfamily = SBB_FONT, fontsize = 10, fontface = "bold"),
  col = sbb_div_col_fun(c(-1, 1)),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = gpar(fontfamily = SBB_FONT, fontsize = 7.0, col = "black"),
  column_names_gp = gpar(fontfamily = SBB_FONT, fontsize = 6.5, col = "black"),
  column_names_rot = 55,
  column_names_side = "bottom",
  rect_gp = gpar(col = "grey72", lwd = 0.45),
  border = "grey45",
  na_col = "white",
  width = unit(60, "mm"),
  height = unit(60, "mm"),
  heatmap_legend_param = list(
    title_gp = gpar(fontfamily = SBB_FONT, fontsize = 7.0, fontface = "bold"),
    labels_gp = gpar(fontfamily = SBB_FONT, fontsize = 6.5),
    at = c(-1, -0.5, 0, 0.5, 1)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (i > j) {
      value <- cor_matrix[i, j]
      grid.text(
        sprintf("%.1f", value),
        x,
        y,
        gp = gpar(
          fontfamily = SBB_FONT,
          fontsize = 5.8,
          col = ifelse(abs(value) >= 0.65, "white", "black")
        )
      )
    }
  }
)
correlation_pdf <- file.path(OUTDIR, "Figure_S2A_candidate_environmental_collinearity.pdf")
grDevices::cairo_pdf(correlation_pdf, width = 4.8, height = 4.7, family = SBB_FONT)
ComplexHeatmap::draw(
  cor_heatmap,
  heatmap_legend_side = "right",
  padding = unit(c(16, 10, 10, 12), "mm")
)
dev.off()

model_long <- model_comparison_tbl %>%
  select(model, R2, adjusted_R2) %>%
  pivot_longer(c(R2, adjusted_R2), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric, R2 = "R<sup>2</sup>", adjusted_R2 = "adj. R<sup>2</sup>"),
    model = factor(
      model,
      levels = c(
        "Full 12-variable set",
        "Historical five-variable set",
        "Predictor-only soil PC1-PC2"
      )
    )
  )
panel_s_b <- ggplot(model_long, aes(model, value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68, color = "grey55", linewidth = 0.25) +
  geom_text(
    aes(label = sprintf("%.2f", value)),
    position = position_dodge(width = 0.76),
    vjust = -0.35,
    size = 2.55,
    color = "black"
  ) +
  scale_fill_manual(values = c("R<sup>2</sup>" = pal_seq[4], "adj. R<sup>2</sup>" = pal_site[["NB"]])) +
  scale_y_continuous(limits = c(0, 0.95), expand = expansion(mult = c(0.01, 0.04))) +
  labs(
    title = "All-17 effect-size comparison",
    subtitle = "R<sup>2</sup> and adjusted R<sup>2</sup>; subset tests are reported separately",
    x = NULL,
    y = "Explained variation",
    fill = NULL
  ) +
  theme_sbb(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 24, hjust = 1, size = 7.0),
    legend.position = "top",
    legend.text = ggtext::element_markdown(),
    plot.subtitle = ggtext::element_markdown(size = 7.0)
  )

vif_plot_data <- vif_tbl %>%
  mutate(
    variable_label = fct_reorder(variable_label, VIF),
    model = recode(
      model,
      `Full 12-variable set` = "Full 12 variables",
      `Historical five-variable set` = "Five variables"
    ),
    model = factor(model, levels = c("Full 12 variables", "Five variables"))
  )
panel_s_c <- ggplot(vif_plot_data, aes(VIF, variable_label, color = flag)) +
  geom_vline(xintercept = 5, color = "grey55", linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = 10, color = "grey25", linewidth = 0.4, linetype = "dashed") +
  geom_segment(aes(x = 1, xend = VIF, yend = variable_label), color = "grey70", linewidth = 0.5) +
  geom_point(size = 2.2) +
  facet_wrap(~model, scales = "free_y", nrow = 1) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 50, 100, 300)) +
  scale_color_manual(values = c("< 5" = "#3D6F97", "5-10" = "#D39A34", "> 10" = "#B84F42")) +
  labs(title = "Variance-inflation diagnostics", x = "VIF (log scale)", y = NULL, color = "VIF") +
  theme_sbb(base_size = 8) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 7.6),
    legend.position = "top"
  )

hp_plot_data <- hp_bootstrap_tbl %>%
  mutate(variable_label = factor(variable_label, levels = rev(unname(PRETTY_ENV[CURRENT_FIVE]))))
panel_s_d <- ggplot(hp_plot_data, aes(y = variable_label)) +
  geom_segment(
    aes(x = ci_low_95, xend = ci_high_95, yend = variable_label),
    color = "grey38",
    linewidth = 0.7
  ) +
  geom_point(
    aes(x = bootstrap_mean_I_perc),
    shape = 21,
    fill = "white",
    color = "grey15",
    size = 2.3,
    stroke = 0.65
  ) +
  geom_point(
    aes(x = point_I_perc),
    shape = 21,
    fill = pal_site[["NB"]],
    color = "white",
    size = 2.6,
    stroke = 0.35
  ) +
  geom_text(
    aes(x = ci_high_95, label = sprintf("%.1f-%.1f", ci_low_95, ci_high_95)),
    hjust = -0.08,
    size = 2.4,
    color = "black"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.28))) +
  labs(
    title = "Five-variable HP sensitivity",
    subtitle = "Whole-gene catalog\nFilled: point; open: pit-profile bootstrap mean",
    x = "Within-model contribution to adj. R<sup>2</sup> (%)",
    y = NULL
  ) +
  theme_sbb(base_size = 8) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = ggtext::element_markdown(),
    plot.subtitle = element_text(size = 6.8)
  )

figure_s2_candidate <- panel_s_b / (panel_s_c | panel_s_d) +
  plot_layout(heights = c(0.72, 1.00)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(family = SBB_FONT, face = "bold", size = 12, color = "black"),
      plot.tag.position = "topleft",
      plot.tag.location = "plot"
    )
  )
ggsave(
  file.path(OUTDIR, "Figure_S2_candidate_model_diagnostics.pdf"),
  figure_s2_candidate,
  width = 18.3,
  height = 14.0,
  units = "cm",
  device = cairo_pdf,
  bg = "white"
)

review_lines <- c(
  "# Figure 2 comment-review reanalysis",
  "",
  "## Status",
  "",
  "This directory contains candidate outputs only. V04 and its commented bilingual review copy were not modified.",
  "",
  "## Design used",
  "",
  "- Final metagenomic sample set: C01 and C03-C18 (17 biological metagenomes; C02 excluded).",
  "- Nahuelbuta: three pit-resolved profiles (NB1-NB3), each sampled at 5-10, 20-40 and 60-80 cm; no technical replicates.",
  "- Surface test: 11 surface samples with unrestricted permutations.",
  "- Depth test: nine Nahuelbuta samples with permutations restricted within site_pit.",
  "",
  "## Candidate primary analysis",
  "",
  sprintf(
    "- Soil PC1 and PC2 explain %.1f%% and %.1f%% of the variance among 12 z-scaled physicochemical measurements (%.1f%% cumulative).",
    100 * soil_variance[1], 100 * soil_variance[2], 100 * sum(soil_variance[1:2])
  ),
  sprintf(
    "- All-17 visual/effect-size model: R2 = %.3f, adjusted R2 = %.3f; no global permutation p-value is promoted.",
    all_stats$R2, all_stats$adjusted_R2
  ),
  sprintf(
    "- Surface model: R2 = %.3f, adjusted R2 = %.3f, p = %.3f.",
    surface_stats$R2, surface_stats$adjusted_R2, surface_stats$p
  ),
  sprintf(
    "- Nahuelbuta pit-restricted model: R2 = %.3f, adjusted R2 = %.3f, p = %.3f.",
    depth_stats$R2, depth_stats$adjusted_R2, depth_stats$p
  ),
  "",
  "## Function-taxonomy concordance",
  "",
  sprintf("- All 17 effect size: Procrustes r = %.3f.", proc_all$table$procrustes_r),
  sprintf(
    "- Surface: r = %.3f, p = %.3f; Nahuelbuta pit-restricted: r = %.3f, p = %.3f.",
    proc_surface$table$procrustes_r,
    proc_surface$table$p,
    proc_depth$table$procrustes_r,
    proc_depth$table$p
  ),
  "",
  "## Interpretation",
  "",
  "The predictor-only soil-axis model is less response-optimised than the historical five-variable model and therefore explains less all-17 variation, but it removes the circular variable-ranking step, reduces the depth model from five predictors to two orthogonal axes, and gives every measured soil property an explicit ecological loading. The surface association remains supported, whereas the depth association is not supported under this lower-dimensional test. The function-taxonomy correspondence is supported in both design-aware subsets.",
  "",
  "The historical five-variable model and hierarchical partitioning remain available as supplementary diagnostics. Their new intervals use a stratified pit/profile-aware cluster bootstrap rather than independent resampling of 17 rows.",
  "",
  "## Candidate manuscript consequence",
  "",
  "Adopting this analysis would require a targeted V05 rewrite of Methods 2.5-2.6, Results 3.2-3.3, the Figure 2/3 legends, and the corresponding Discussion sentences. That manuscript mutation remains pending author acceptance of the candidate analysis."
)
writeLines(review_lines, file.path(OUTDIR, "ANALYSIS_README.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(OUTDIR, "R_sessionInfo.txt"))

message("Candidate analysis complete: ", normalizePath(OUTDIR))
