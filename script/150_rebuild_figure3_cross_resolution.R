# =============================================================================
# V05 Figure 3: cross-resolution consistency under the response-independent
# soil-axis model.
#
# This script applies the soil PC1-PC2 model accepted for V05 to the whole gene
# catalog and the fractionally assigned KO profile. It also recomputes design-aware
# concordance with recovered-MAG taxonomy and taxonomy-conditioned sensitivity.
# No legacy V04 output is overwritten.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(vegan)
  library(permute)
  library(patchwork)
  library(ggtext)
})

options(scipen = 999)
set.seed(42)
source("script/_theme_sbb.R")
SBB_FONT <- "Arial"

IN_DIR <- "analysis/V05_comment_review_2026-08-17/figure2_reanalysis"
OUT_DIR <- "analysis/V05_comment_review_2026-08-17/figure3_cross_resolution"
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

PATHS <- c(
  reads = "data/processed/metagenome_output/Genecatalog/Nmapped_reads.tsv",
  lengths = "data/processed/metagenome_output/Genecatalog/seq_len_nucl.tsv",
  annotation = "data/processed/metagenome_output/Genecatalog/eggNog.tsv",
  mag_counts = "data/processed/metagenome_output/MAGs/counts_genomes.tsv",
  mag_quality = "data/processed/metagenome_output/MAGs/genome_quality.tsv",
  metadata = "data/metadata/env_2024_meta.tsv",
  soil_scores = file.path(IN_DIR, "soil_pca_scores.tsv")
)
missing_inputs <- PATHS[!file.exists(PATHS)]
if (length(missing_inputs)) {
  stop("Missing required inputs: ", paste(missing_inputs, collapse = ", "))
}

SAMPLES <- c("C01", sprintf("C%02d", 3:18))
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
PROFILE_LEVELS <- c("Whole gene catalog", "Fractional KO catalog")
PROFILE_LABELS <- c(
  `Whole gene catalog` = "Whole gene catalog",
  `Fractional KO catalog` = "KO catalog"
)
PROFILE_COLORS <- c(
  `Whole gene catalog` = SBB_GENERIC_BLUE,
  `Fractional KO catalog` = pal_category[["terracotta"]]
)
SUBSET_LEVELS <- c("all_17", "surface_11", "nahuelbuta_depth_9")
SUBSET_LABELS <- c(
  all_17 = "All 17",
  surface_11 = "Surface",
  nahuelbuta_depth_9 = "NA depth"
)

write_tsv_checked <- function(x, filename) {
  path <- file.path(OUT_DIR, filename)
  readr::write_tsv(x, path, na = "NA")
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Failed to write output: ", path)
  }
  path
}

format_p <- function(value) {
  case_when(
    is.na(value) ~ "effect size",
    value < 0.001 ~ "italic(p) < 0.001",
    TRUE ~ sprintf("italic(p) == %.3f", value)
  )
}

message("Reading the frozen soil scores and design map...")
metadata <- read_tsv(PATHS[["metadata"]], show_col_types = FALSE, na = c("", "NA")) %>%
  filter(sample_id %in% SAMPLES) %>%
  arrange(match(sample_id, SAMPLES)) %>%
  mutate(
    site = factor(site, levels = SITE_LEVELS),
    depth_label = factor(depth_label, levels = DEPTH_LEVELS),
    site_pit = as.character(site_pit)
  )
soil_scores <- read_tsv(PATHS[["soil_scores"]], show_col_types = FALSE) %>%
  filter(sample_id %in% SAMPLES) %>%
  arrange(match(sample_id, SAMPLES)) %>%
  rename(SoilPC1 = PC1, SoilPC2 = PC2)
if (!identical(metadata$sample_id, SAMPLES) || !identical(soil_scores$sample_id, SAMPLES)) {
  stop("Metadata, soil scores, and the frozen 17-sample order do not agree.")
}
soil_model_data <- soil_scores %>%
  select(sample_id, SoilPC1, SoilPC2) %>%
  column_to_rownames("sample_id") %>%
  as.data.frame()

index_all <- seq_along(SAMPLES)
index_surface <- which(metadata$depth_label == "5_10cm")
index_depth <- which(metadata$site == "NB")
depth_control <- how(
  nperm = 999,
  blocks = factor(metadata$site_pit[index_depth])
)

permutation_for <- function(subset) {
  if (subset == "all_17") return(NULL)
  if (subset == "surface_11") return(999)
  if (subset == "nahuelbuta_depth_9") return(depth_control)
  stop("Unknown subset: ", subset)
}

indices_for <- function(subset) {
  if (subset == "all_17") return(index_all)
  if (subset == "surface_11") return(index_surface)
  if (subset == "nahuelbuta_depth_9") return(index_depth)
  stop("Unknown subset: ", subset)
}

permutation_label <- function(subset) {
  if (subset == "all_17") return("effect size only; no mixed-design permutation")
  if (subset == "surface_11") return("999 unrestricted permutations")
  "complete/Monte Carlo permutations restricted within site_pit"
}

message("Constructing whole-gene and fractionally assigned KO feature profiles...")
counts <- fread(PATHS[["reads"]], select = c("GeneName", SAMPLES))
lengths <- fread(
  PATHS[["lengths"]],
  header = FALSE,
  col.names = c("GeneName", "Length")
)
annotation <- fread(PATHS[["annotation"]], select = c("GeneName", "KO"))

gene_rpk <- merge(counts, lengths, by = "GeneName")
gene_rpk <- gene_rpk[Length > 0]
gene_rpk[, total_reads := rowSums(.SD, na.rm = TRUE), .SDcols = SAMPLES]
gene_rpk <- gene_rpk[total_reads > 0]
for (sample_id in SAMPLES) {
  set(
    gene_rpk,
    j = sample_id,
    value = gene_rpk[[sample_id]] / (gene_rpk$Length / 1000)
  )
}
gene_rpk[, c("Length", "total_reads") := NULL]
rm(counts, lengths)
invisible(gc())

profile_from_table <- function(table, label) {
  feature_matrix <- as.matrix(table[, ..SAMPLES])
  storage.mode(feature_matrix) <- "double"
  feature_matrix <- feature_matrix[rowSums(feature_matrix, na.rm = TRUE) > 0, , drop = FALSE]
  community <- decostand(t(feature_matrix), method = "hellinger")
  rownames(community) <- SAMPLES
  list(
    label = label,
    n_features = ncol(community),
    community = community,
    distance = vegdist(community, method = "bray")
  )
}

whole_gene_profile <- profile_from_table(gene_rpk, "Whole gene catalog")

ko_map <- annotation[
  !is.na(KO) & KO != "" & KO != "-",
  .(KO = trimws(unlist(strsplit(gsub("ko:", "", KO, fixed = TRUE), ",", fixed = TRUE)))),
  by = GeneName
]
ko_map <- unique(ko_map[grepl("^K[0-9]{5}$", KO)])
ko_map[, n_valid_KO := .N, by = GeneName]
fractional_joined <- merge(ko_map[, .(GeneName, KO, n_valid_KO)], gene_rpk, by = "GeneName")
for (sample_id in SAMPLES) {
  set(
    fractional_joined,
    j = sample_id,
    value = fractional_joined[[sample_id]] / fractional_joined$n_valid_KO
  )
}
fractional_ko_rpk <- fractional_joined[
  , lapply(.SD, sum, na.rm = TRUE),
  by = KO,
  .SDcols = SAMPLES
]
fractional_ko_profile <- profile_from_table(fractional_ko_rpk, "Fractional KO catalog")

ko_summary <- tibble(
  metric = c(
    "nonzero genes in whole-gene profile",
    "genes with at least one valid KO",
    "genes with multiple valid KOs",
    "detected fractionally assigned KO features"
  ),
  value = c(
    whole_gene_profile$n_features,
    uniqueN(ko_map$GeneName),
    uniqueN(ko_map[n_valid_KO > 1]$GeneName),
    fractional_ko_profile$n_features
  )
)
write_tsv_checked(ko_summary, "ko_annotation_summary.tsv")

profiles <- list(whole_gene_profile, fractional_ko_profile)
rm(annotation, fractional_joined, fractional_ko_rpk, ko_map, gene_rpk)
invisible(gc())

fit_soil_axis <- function(profile, subset) {
  profile_label <- profile$label
  profile_features <- profile$n_features
  indices <- indices_for(subset)
  predictors <- soil_model_data[indices, , drop = FALSE]
  response <- profile$community[indices, , drop = FALSE]
  model <- dbrda(
    response ~ SoilPC1 + SoilPC2,
    data = predictors,
    distance = "bray"
  )
  r2 <- RsquareAdj(model)
  p_value <- NA_real_
  permutation <- permutation_for(subset)
  if (!is.null(permutation)) {
    set.seed(42)
    p_value <- anova.cca(model, permutations = permutation)$`Pr(>F)`[1]
  }
  tibble(
    profile = profile_label,
    subset = subset,
    n_samples = length(indices),
    n_features = profile_features,
    R2 = r2$r.squared,
    adjusted_R2 = r2$adj.r.squared,
    p = p_value,
    permutation_scheme = permutation_label(subset)
  )
}

soil_axis_models <- map_dfr(profiles, function(profile) {
  map_dfr(SUBSET_LEVELS, ~ fit_soil_axis(profile, .x))
})
write_tsv_checked(soil_axis_models, "cross_resolution_soil_axis_models.tsv")

run_procrustes <- function(profile_x, profile_y, subset, coordinates = FALSE) {
  indices <- indices_for(subset)
  x_distance <- vegdist(profile_x$community[indices, , drop = FALSE], method = "bray")
  y_distance <- vegdist(profile_y$community[indices, , drop = FALSE], method = "bray")
  x_ordination <- wcmdscale(x_distance, k = 2, eig = TRUE)
  y_ordination <- wcmdscale(y_distance, k = 2, eig = TRUE)
  x_coordinates <- x_ordination$points[, 1:2, drop = FALSE]
  y_coordinates <- y_ordination$points[, 1:2, drop = FALSE]
  fit <- procrustes(x_coordinates, y_coordinates, symmetric = TRUE)
  statistic <- sqrt(1 - fit$ss)
  p_value <- NA_real_
  permutation <- permutation_for(subset)
  if (!is.null(permutation)) {
    set.seed(42)
    test <- protest(
      x_coordinates,
      y_coordinates,
      permutations = permutation,
      symmetric = TRUE
    )
    statistic <- unname(test$t0)
    p_value <- test$signif
  }

  coordinate_table <- NULL
  if (coordinates) {
    coordinate_table <- tibble(
      sample_id = rownames(fit$X),
      x_axis1 = fit$X[, 1],
      x_axis2 = fit$X[, 2],
      y_axis1 = fit$Yrot[, 1],
      y_axis2 = fit$Yrot[, 2]
    ) %>%
      left_join(
        metadata %>% select(sample_id, site, site_pit, depth_label),
        by = "sample_id"
      )
  }

  list(
    summary = tibble(
      profile_x = profile_x$label,
      profile_y = profile_y$label,
      subset = subset,
      n_samples = length(indices),
      procrustes_r = statistic,
      p = p_value,
      permutation_scheme = permutation_label(subset)
    ),
    coordinates = coordinate_table
  )
}

gene_ko_results <- map(SUBSET_LEVELS, function(subset) {
  run_procrustes(
    whole_gene_profile,
    fractional_ko_profile,
    subset,
    coordinates = subset == "all_17"
  )
})
gene_ko_summary <- map_dfr(gene_ko_results, "summary")
gene_ko_coordinates <- gene_ko_results[[1]]$coordinates
write_tsv_checked(gene_ko_summary, "whole_gene_ko_procrustes.tsv")
write_tsv_checked(gene_ko_coordinates, "whole_gene_ko_procrustes_coordinates.tsv")

message("Constructing recovered-MAG taxonomic distances...")
mag_counts <- fread(PATHS[["mag_counts"]])
setnames(mag_counts, "Sample", "MAG")
mag_quality <- fread(
  PATHS[["mag_quality"]],
  select = c("Bin Id", "Genome_Size")
)
setnames(mag_quality, "Bin Id", "MAG")
mag_table <- merge(mag_counts, mag_quality, by = "MAG", all.x = TRUE)
if (anyNA(mag_table$Genome_Size) || any(mag_table$Genome_Size <= 0)) {
  stop("Missing or invalid MAG genome sizes.")
}
for (sample_id in SAMPLES) {
  set(
    mag_table,
    j = sample_id,
    value = mag_table[[sample_id]] / (mag_table$Genome_Size / 1000)
  )
}
mag_profile <- profile_from_table(mag_table[, c("MAG", SAMPLES), with = FALSE], "Recovered MAG taxonomy")

function_taxonomy <- map_dfr(profiles, function(profile) {
  map_dfr(SUBSET_LEVELS, function(subset) {
    run_procrustes(profile, mag_profile, subset)$summary
  })
})
write_tsv_checked(function_taxonomy, "function_taxonomy_cross_resolution.tsv")

fit_conditioned_axis <- function(profile, subset, n_axes) {
  indices <- indices_for(subset)
  response_distance <- vegdist(profile$community[indices, , drop = FALSE], method = "bray")
  tax_distance <- vegdist(mag_profile$community[indices, , drop = FALSE], method = "bray")
  tax_ordination <- wcmdscale(tax_distance, k = 4, eig = TRUE)
  tax_axes <- as.data.frame(tax_ordination$points[, seq_len(n_axes), drop = FALSE])
  colnames(tax_axes) <- paste0("TaxPCoA", seq_len(n_axes))
  model_data <- cbind(
    soil_model_data[indices, , drop = FALSE],
    tax_axes
  )
  condition_terms <- paste(colnames(tax_axes), collapse = " + ")
  formula <- as.formula(
    paste0(
      "response_distance ~ SoilPC1 + SoilPC2 + Condition(",
      condition_terms,
      ")"
    ),
    env = environment()
  )
  model <- dbrda(formula, data = model_data)
  r2 <- RsquareAdj(model)
  p_value <- NA_real_
  permutation <- permutation_for(subset)
  if (!is.null(permutation)) {
    set.seed(42)
    p_value <- anova.cca(model, permutations = permutation)$`Pr(>F)`[1]
  }
  positive_eigenvalues <- tax_ordination$eig[tax_ordination$eig > 0]
  tibble(
    profile = profile$label,
    subset = subset,
    n_samples = length(indices),
    taxonomic_axes_conditioned = n_axes,
    taxonomic_variance_fraction_conditioned =
      sum(positive_eigenvalues[seq_len(n_axes)]) / sum(positive_eigenvalues),
    unique_environment_R2 = r2$r.squared,
    unique_environment_adjusted_R2 = r2$adj.r.squared,
    p = p_value,
    permutation_scheme = permutation_label(subset)
  )
}

taxonomy_conditioned <- map_dfr(profiles, function(profile) {
  map_dfr(SUBSET_LEVELS, function(subset) {
    map_dfr(1:4, ~ fit_conditioned_axis(profile, subset, .x))
  })
})
write_tsv_checked(
  taxonomy_conditioned,
  "taxonomy_conditioned_axis_sensitivity.tsv"
)

# -----------------------------------------------------------------------------
# Figure 3A: whole-gene versus fractionally assigned KO concordance.
# -----------------------------------------------------------------------------
gene_ko_all <- gene_ko_summary %>% filter(subset == "all_17") %>% slice_head(n = 1)
panel_a <- ggplot(gene_ko_coordinates) +
  geom_segment(
    aes(x = x_axis1, y = x_axis2, xend = y_axis1, yend = y_axis2),
    colour = "grey70",
    linewidth = 0.42,
    lineend = "round"
  ) +
  geom_point(
    aes(x = x_axis1, y = x_axis2, fill = site, shape = depth_label),
    size = 2.75,
    stroke = 0.52,
    colour = "grey12"
  ) +
  geom_point(
    aes(x = y_axis1, y = y_axis2, colour = site, shape = depth_label),
    size = 2.75,
    stroke = 0.78,
    fill = "white"
  ) +
  scale_fill_manual(
    name = "Site",
    values = SITE_COLORS,
    breaks = SITE_LEVELS,
    labels = SITE_LABELS
  ) +
  scale_colour_manual(
    name = "Site",
    values = SITE_COLORS,
    breaks = SITE_LEVELS,
    labels = SITE_LABELS,
    guide = "none"
  ) +
  scale_shape_manual(
    name = "Depth",
    values = DEPTH_SHAPES,
    breaks = DEPTH_LEVELS,
    labels = DEPTH_LABELS
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, colour = "grey12", size = 3.0)
    ),
    shape = guide_legend(
      override.aes = list(fill = "white", colour = "grey12", size = 3.0)
    )
  ) +
  labs(
    title = "Whole gene-KO concordance",
    subtitle = sprintf(
      "Filled: genes; hollow: fractional KOs<br>All 17: r = %.3f",
      gene_ko_all$procrustes_r
    ),
    x = "Procrustes axis 1",
    y = "Procrustes axis 2"
  ) +
  coord_equal() +
  theme_sbb(base_size = 9.2) +
  theme(
    legend.position = "bottom",
    plot.subtitle = element_markdown(size = 8.2, lineheight = 1.05),
    legend.box = "vertical",
    legend.margin = margin(t = -2),
    plot.margin = margin(7, 8, 5, 8)
  )

# -----------------------------------------------------------------------------
# Figure 3B: the same soil-axis association at gene and KO resolution.
# -----------------------------------------------------------------------------
panel_b_data <- soil_axis_models %>%
  mutate(
    profile = factor(profile, levels = PROFILE_LEVELS),
    subset = factor(subset, levels = rev(SUBSET_LEVELS)),
    subset_label = factor(
      unname(SUBSET_LABELS[as.character(subset)]),
      levels = unname(SUBSET_LABELS[rev(SUBSET_LEVELS)])
    ),
    p_label = case_when(
      is.na(p) ~ "effect size",
      p < 0.001 ~ "<i>p</i> < 0.001",
      TRUE ~ sprintf("<i>p</i> = %.3f", p)
    )
  )

panel_b <- ggplot(
  panel_b_data,
  aes(x = adjusted_R2, y = subset_label, colour = profile, shape = profile)
) +
  geom_vline(xintercept = 0, colour = "grey76", linewidth = 0.35) +
  geom_point(
    position = position_dodge(width = 0.42),
    size = 3.0,
    stroke = 0.75
  ) +
  ggtext::geom_richtext(
    aes(label = p_label),
    position = position_dodge(width = 0.42),
    hjust = -0.10,
    size = 2.30,
    family = SBB_FONT,
    colour = "grey25",
    fill = NA,
    label.colour = NA,
    label.padding = unit(rep(0, 4), "pt"),
    show.legend = FALSE
  ) +
  scale_colour_manual(
    name = NULL,
    values = PROFILE_COLORS,
    breaks = PROFILE_LEVELS,
    labels = PROFILE_LABELS
  ) +
  scale_shape_manual(
    name = NULL,
    values = c(`Whole gene catalog` = 16, `Fractional KO catalog` = 17),
    breaks = PROFILE_LEVELS,
    labels = PROFILE_LABELS
  ) +
  scale_x_continuous(
    limits = c(min(-0.02, min(panel_b_data$adjusted_R2) - 0.01), 0.48),
    breaks = seq(0, 0.4, 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Soil-axis association",
    subtitle = "Response-independent soil PC1-PC2 model",
    x = expression("Adjusted " * R^2),
    y = NULL
  ) +
  theme_sbb(base_size = 9.2) +
  theme(
    legend.position = "bottom",
    plot.subtitle = element_text(size = 8.2),
    legend.box = "horizontal",
    legend.text = element_text(size = 7.4),
    legend.key.width = unit(4.0, "mm"),
    legend.spacing.x = unit(1.0, "mm"),
    legend.margin = margin(t = -2),
    axis.text.y = element_text(size = 8.2),
    plot.margin = margin(7, 12, 5, 8)
  )

# -----------------------------------------------------------------------------
# Figure 3C: unique soil-axis fraction after conditioning on MAG taxonomy.
# -----------------------------------------------------------------------------
panel_c_data <- taxonomy_conditioned %>%
  filter(subset != "all_17") %>%
  mutate(
    profile = factor(profile, levels = PROFILE_LEVELS),
    subset_label = factor(
      unname(SUBSET_LABELS[subset]),
      levels = c("Surface", "NA depth")
    ),
    label_vjust = case_when(
      subset_label == "Surface" & profile == "Whole gene catalog" ~ -0.75,
      subset_label == "Surface" ~ 1.55,
      subset_label == "NA depth" & profile == "Whole gene catalog" ~ 1.55,
      TRUE ~ -0.75
    )
  )

two_axis_labels <- panel_c_data %>%
  filter(taxonomic_axes_conditioned == 2) %>%
  mutate(
    label = case_when(
      p < 0.001 ~ "<i>p</i> < 0.001",
      TRUE ~ sprintf("<i>p</i> = %.3f", p)
    )
  )

panel_c <- ggplot(
  panel_c_data,
  aes(
    x = taxonomic_axes_conditioned,
    y = unique_environment_adjusted_R2,
    colour = profile,
    shape = profile,
    group = profile
  )
) +
  geom_hline(yintercept = 0, colour = "grey68", linewidth = 0.38, linetype = "dashed") +
  geom_line(linewidth = 0.68) +
  geom_point(size = 2.65, stroke = 0.72) +
  ggtext::geom_richtext(
    data = two_axis_labels,
    aes(label = label, vjust = label_vjust),
    hjust = -0.08,
    size = 2.35,
    family = SBB_FONT,
    fill = NA,
    label.colour = NA,
    label.padding = unit(rep(0, 4), "pt"),
    show.legend = FALSE
  ) +
  facet_wrap(~ subset_label, nrow = 1) +
  scale_colour_manual(
    name = NULL,
    values = PROFILE_COLORS,
    breaks = PROFILE_LEVELS,
    labels = PROFILE_LABELS
  ) +
  scale_shape_manual(
    name = NULL,
    values = c(`Whole gene catalog` = 16, `Fractional KO catalog` = 17),
    breaks = PROFILE_LEVELS,
    labels = PROFILE_LABELS
  ) +
  scale_x_continuous(breaks = 1:4, minor_breaks = NULL) +
  labs(
    title = "Soil fraction after conditioning on MAG taxonomy",
    subtitle = "Labels show the design-aware test after conditioning on two taxonomic axes",
    x = "MAG taxonomic PCoA axes conditioned",
    y = expression("Unique soil-axis adjusted " * R^2)
  ) +
  theme_sbb(base_size = 9.2) +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(size = 8.2),
    strip.text = element_text(size = 8.5, face = "bold"),
    panel.spacing.x = unit(8, "mm"),
    legend.margin = margin(t = -2),
    plot.margin = margin(7, 9, 6, 8)
  )

ggsave(
  file.path(PANEL_DIR, "Figure_3_A.pdf"),
  panel_a,
  width = 8.45,
  height = 8.25,
  units = "cm",
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(PANEL_DIR, "Figure_3_B.pdf"),
  panel_b,
  width = 8.45,
  height = 8.25,
  units = "cm",
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(PANEL_DIR, "Figure_3_C.pdf"),
  panel_c,
  width = 17.45,
  height = 7.75,
  units = "cm",
  device = cairo_pdf,
  bg = "white"
)

figure_3 <- (panel_a + panel_b) / panel_c +
  plot_layout(heights = c(1.05, 0.95)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(
      family = SBB_FONT,
      face = "bold",
      size = 14,
      colour = "black"
    ),
    plot.tag.position = c(0.01, 0.99)
  )

ggsave(
  file.path(OUT_DIR, "Figure_3_candidate_cross_resolution.pdf"),
  figure_3,
  width = 17.8,
  height = 16.5,
  units = "cm",
  device = cairo_pdf,
  bg = "white"
)

readme <- c(
  "# V05 Figure 3 cross-resolution reconstruction",
  "",
  "- The response-independent soil PC1-PC2 model is applied unchanged to the whole gene catalog and a KO catalog built by fractionally assigning genes with multiple valid KOs.",
  "- All-17 coefficients are effect sizes; no mixed-design permutation p value is reported.",
  "- Surface tests use 999 unrestricted permutations; Nahuelbuta tests restrict permutations within site_pit.",
  "- Panel A shows whole-gene/KO Procrustes concordance.",
  "- Panel B compares soil-axis adjusted R2 at the two functional resolutions.",
  "- Panel C evaluates the unique soil-axis fraction after conditioning on one to four MAG taxonomic PCoA axes.",
  "- The figure uses only vector PDF outputs; no persistent PNG is generated.",
  "- V04 and the accepted V05 Figure 2 source remain untouched."
)
writeLines(readme, file.path(OUT_DIR, "README.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(OUT_DIR, "R_sessionInfo.txt"))

message("V05 Figure 3 reconstruction complete: ", normalizePath(OUT_DIR))
