# =============================================================================
# Candidate V06 ecological-spine audit
#
# This script asks the two direct design-level questions that precede any
# environmental-driver interpretation:
#   1. Do whole-gene, fractional-KO, and recovered-MAG configurations differ
#      among the four sampled surface sites?
#   2. Do the same configurations differ among Nahuelbuta depths after
#      restricting permutations within the three pit profiles?
#
# The analysis is non-destructive. It does not replace the V05 soil-axis model
# or write into any V05 submission folder.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(vegan)
  library(permute)
})

options(scipen = 999)
set.seed(42)

OUT_DIR <- "analysis/V06_ecological_spine_audit_2026-08-18"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PATHS <- c(
  reads = "data/processed/metagenome_output/Genecatalog/Nmapped_reads.tsv",
  lengths = "data/processed/metagenome_output/Genecatalog/seq_len_nucl.tsv",
  annotation = "data/processed/metagenome_output/Genecatalog/eggNog.tsv",
  mag_counts = "data/processed/metagenome_output/MAGs/counts_genomes.tsv",
  mag_quality = "data/processed/metagenome_output/MAGs/genome_quality.tsv",
  metadata = "data/metadata/env_2024_meta.tsv"
)

missing_inputs <- PATHS[!file.exists(PATHS)]
if (length(missing_inputs)) {
  stop("Missing required inputs: ", paste(missing_inputs, collapse = ", "))
}

SAMPLES <- c("C01", sprintf("C%02d", 3:18))
SITE_LEVELS <- c("AZ", "SG", "LC", "NB")
DEPTH_LEVELS <- c("5_10cm", "20_40cm", "60_80cm")
PROFILE_LEVELS <- c(
  "Whole gene catalog",
  "Fractional KO catalog",
  "Genome-size-normalised MAG assemblage"
)

write_tsv_checked <- function(x, filename) {
  path <- file.path(OUT_DIR, filename)
  readr::write_tsv(x, path, na = "NA")
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Failed to write output: ", path)
  }
  path
}

metadata <- read_tsv(PATHS[["metadata"]], show_col_types = FALSE, na = c("", "NA")) %>%
  filter(sample_id %in% SAMPLES) %>%
  arrange(match(sample_id, SAMPLES)) %>%
  mutate(
    site = factor(site, levels = SITE_LEVELS),
    depth_label = factor(depth_label, levels = DEPTH_LEVELS),
    site_pit = factor(site_pit)
  )

if (!identical(metadata$sample_id, SAMPLES)) {
  stop("Metadata and expected metagenome order do not agree.")
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
) %>%
  arrange(site_pit, match(depth_label, DEPTH_LEVELS))

observed_nb_map <- metadata %>%
  filter(site == "NB") %>%
  transmute(
    site_pit = as.character(site_pit),
    depth_label = as.character(depth_label),
    sample_id
  ) %>%
  arrange(site_pit, match(depth_label, DEPTH_LEVELS))

if (!identical(observed_nb_map, expected_nb_map)) {
  stop("Nahuelbuta pit-depth map differs from the authoritative design map.")
}
write_tsv_checked(observed_nb_map, "nahuelbuta_pit_depth_map.tsv")

profile_from_feature_table <- function(table, feature_id, label) {
  feature_matrix <- as.matrix(table[, ..SAMPLES])
  storage.mode(feature_matrix) <- "double"
  keep <- rowSums(feature_matrix, na.rm = TRUE) > 0
  feature_matrix <- feature_matrix[keep, , drop = FALSE]
  community <- decostand(t(feature_matrix), method = "hellinger")
  rownames(community) <- SAMPLES
  list(
    label = label,
    feature_id = feature_id,
    n_features = ncol(community),
    community = community
  )
}

message("Building the whole-gene RPK profile...")
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

whole_gene_profile <- profile_from_feature_table(
  gene_rpk,
  "GeneName",
  "Whole gene catalog"
)

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
  gene_rpk,
  by = "GeneName",
  allow.cartesian = TRUE
)
for (sample_id in SAMPLES) {
  set(
    ko_joined,
    j = sample_id,
    value = ko_joined[[sample_id]] / ko_joined$n_valid_KO
  )
}
ko_rpk <- ko_joined[
  , lapply(.SD, sum, na.rm = TRUE),
  by = KO,
  .SDcols = SAMPLES
]
fractional_ko_profile <- profile_from_feature_table(
  ko_rpk,
  "KO",
  "Fractional KO catalog"
)

ko_summary <- tibble(
  metric = c(
    "nonzero whole-gene features",
    "genes with at least one valid KO",
    "genes with multiple valid KOs",
    "detected fractional-KO features"
  ),
  value = c(
    whole_gene_profile$n_features,
    uniqueN(ko_map$GeneName),
    uniqueN(ko_map[n_valid_KO > 1]$GeneName),
    fractional_ko_profile$n_features
  )
)
write_tsv_checked(ko_summary, "ko_annotation_summary.tsv")

rm(annotation, gene_rpk, ko_joined, ko_map, ko_rpk)
invisible(gc())

message("Building the genome-size-normalised MAG profile...")
mag_counts <- fread(PATHS[["mag_counts"]])
setnames(mag_counts, "Sample", "MAG")
mag_quality <- fread(
  PATHS[["mag_quality"]],
  select = c("Bin Id", "Genome_Size")
)
setnames(mag_quality, "Bin Id", "MAG")
mag_rpk <- merge(mag_counts, mag_quality, by = "MAG", all.x = TRUE)
if (anyNA(mag_rpk$Genome_Size) || any(mag_rpk$Genome_Size <= 0)) {
  stop("Missing or invalid MAG genome sizes.")
}
for (sample_id in SAMPLES) {
  set(
    mag_rpk,
    j = sample_id,
    value = mag_rpk[[sample_id]] / (mag_rpk$Genome_Size / 1000)
  )
}
mag_rpk[, Genome_Size := NULL]
mag_profile <- profile_from_feature_table(
  mag_rpk,
  "MAG",
  "Genome-size-normalised MAG assemblage"
)
rm(mag_counts, mag_quality, mag_rpk)
invisible(gc())

profiles <- list(whole_gene_profile, fractional_ko_profile, mag_profile)

index_surface <- which(metadata$depth_label == "5_10cm")
index_depth <- which(metadata$site == "NB")
depth_control <- how(
  nperm = 999,
  blocks = droplevels(metadata$site_pit[index_depth])
)

run_location_test <- function(profile, subset) {
  profile_label <- profile$label
  profile_features <- profile$n_features

  if (subset == "surface_site") {
    indices <- index_surface
    group <- droplevels(metadata$site[indices])
    permutations <- 999
    scheme <- "999 unrestricted permutations among 11 surface samples"
    design_note <- "Site-associated contrast; climate category is confounded with named site"
  } else if (subset == "nahuelbuta_depth") {
    indices <- index_depth
    group <- droplevels(metadata$depth_label[indices])
    permutations <- depth_control
    scheme <- "permutations restricted within three Nahuelbuta pit profiles"
    design_note <- "Depth contrast within one site, blocked by site_pit"
  } else {
    stop("Unknown subset: ", subset)
  }

  community <- profile$community[indices, , drop = FALSE]
  distance <- vegdist(community, method = "bray")
  model_data <- tibble(group = group)

  set.seed(42)
  location_fit <- adonis2(
    distance ~ group,
    data = model_data,
    permutations = permutations,
    by = NULL
  )

  dispersion_fit <- betadisper(distance, group, type = "median", bias.adjust = TRUE)
  set.seed(42)
  dispersion_test <- permutest(
    dispersion_fit,
    permutations = permutations,
    pairwise = FALSE
  )

  list(
    location = tibble(
      profile = profile_label,
      subset = subset,
      n_samples = length(indices),
      n_features = profile_features,
      groups = nlevels(group),
      R2 = location_fit$R2[1],
      F = location_fit$F[1],
      p = location_fit$`Pr(>F)`[1],
      permutation_scheme = scheme,
      design_note = design_note
    ),
    dispersion = tibble(
      profile = profile_label,
      subset = subset,
      n_samples = length(indices),
      groups = nlevels(group),
      F = dispersion_test$tab$F[1],
      p = dispersion_test$tab$`Pr(>F)`[1],
      permutation_scheme = scheme
    )
  )
}

test_results <- map(profiles, function(profile) {
  map(
    c("surface_site", "nahuelbuta_depth"),
    function(subset) run_location_test(profile, subset)
  )
}) %>%
  flatten()

location_results <- map_dfr(test_results, "location") %>%
  mutate(
    profile = factor(profile, levels = PROFILE_LEVELS),
    q_within_subset = ave(p, subset, FUN = function(x) p.adjust(x, method = "BH"))
  ) %>%
  arrange(subset, profile) %>%
  mutate(profile = as.character(profile))

dispersion_results <- map_dfr(test_results, "dispersion") %>%
  mutate(
    profile = factor(profile, levels = PROFILE_LEVELS),
    q_within_subset = ave(p, subset, FUN = function(x) p.adjust(x, method = "BH"))
  ) %>%
  arrange(subset, profile) %>%
  mutate(profile = as.character(profile))

write_tsv_checked(location_results, "direct_location_permanova.tsv")
write_tsv_checked(dispersion_results, "direct_location_dispersion.tsv")

ordination_coordinates <- map_dfr(profiles, function(profile) {
  map_dfr(c("surface_site", "nahuelbuta_depth"), function(subset) {
    indices <- if (subset == "surface_site") index_surface else index_depth
    distance <- vegdist(profile$community[indices, , drop = FALSE], method = "bray")
    ordination <- wcmdscale(distance, k = 2, eig = TRUE, add = "lingoes")
    positive <- ordination$eig[ordination$eig > 0]
    variance <- positive / sum(positive)
    tibble(
      profile = profile$label,
      subset = subset,
      sample_id = rownames(ordination$points),
      axis1 = ordination$points[, 1],
      axis2 = ordination$points[, 2],
      axis1_variance = variance[1],
      axis2_variance = variance[2]
    ) %>%
      left_join(
        metadata %>%
          transmute(
            sample_id,
            site = as.character(site),
            site_pit = as.character(site_pit),
            depth_label = as.character(depth_label)
          ),
        by = "sample_id"
      )
  })
})
write_tsv_checked(ordination_coordinates, "direct_location_pcoa_coordinates.tsv")

readme_lines <- c(
  "# V06 ecological-spine audit",
  "",
  "This candidate analysis tests direct surface-site and Nahuelbuta-depth turnover before interpreting environmental predictors. It does not replace V05.",
  "",
  "## Design",
  "",
  "- Surface: 11 independent biological samples across AZ, SG, LC, and NA; site and climate position are confounded.",
  "- Nahuelbuta depth: three pit-resolved profiles; permutations are restricted within site_pit.",
  "- Whole-gene and fractional-KO profiles use length-normalised abundance followed by Hellinger transformation and Bray-Curtis distance.",
  "- The MAG profile is genome-size normalised before the same Hellinger/Bray-Curtis workflow.",
  "- PERMANOVA location and betadisper tests are reported together. BH q values are calculated across the three resolution layers within each design subset.",
  "",
  "## Output interpretation",
  "",
  "- A supported location test with unsupported dispersion supports configuration turnover at that resolution.",
  "- A supported dispersion test prevents a simple centroid-shift interpretation.",
  "- Agreement across layers is robustness across representations derived from the same metagenomes, not independent experimental replication.",
  "",
  "## Files",
  "",
  "- direct_location_permanova.tsv",
  "- direct_location_dispersion.tsv",
  "- direct_location_pcoa_coordinates.tsv",
  "- ko_annotation_summary.tsv",
  "- nahuelbuta_pit_depth_map.tsv",
  "- R_sessionInfo.txt"
)
writeLines(readme_lines, file.path(OUT_DIR, "README.md"))
capture.output(sessionInfo(), file = file.path(OUT_DIR, "R_sessionInfo.txt"))

message("Completed V06 ecological-spine audit: ", OUT_DIR)
print(location_results)
print(dispersion_results)
