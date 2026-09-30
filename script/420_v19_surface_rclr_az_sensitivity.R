#!/usr/bin/env Rscript
# V19 targeted surface-site sensitivity analysis.
#
# The prespecified Bray-Curtis and rCLR/Aitchison analyses remain the primary
# and compositional-sensitivity endpoints, respectively. This script asks two
# bounded influence questions without changing those endpoints:
#   1. Does the site pattern persist among SG, LC, and NA after excluding AZ?
#   2. Is it stable when all four sites contribute exactly two samples?

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(vegan)
})

options(scipen = 999)
set.seed(42)

ROOT <- normalizePath(Sys.getenv("CHILE_PROJECT_ROOT", getwd()), mustWork = TRUE)
setwd(ROOT)
OUT <- file.path(ROOT, "analysis/V19_peer_review_targeted_revision_2026-09-01/robustness")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SAMPLES <- c("C01", sprintf("C%02d", 3:18))
PATHS <- c(
  counts = "data/processed/metagenome_output/Genecatalog/Nmapped_reads.tsv",
  lengths = "data/processed/metagenome_output/Genecatalog/seq_len_nucl.tsv",
  annotation = "data/processed/metagenome_output/Genecatalog/eggNog.tsv",
  metadata = "data/metadata/env_2024_meta.tsv"
)
stopifnot(all(file.exists(PATHS)))

write_tsv_checked <- function(x, filename) {
  path <- file.path(OUT, filename)
  readr::write_tsv(x, path, na = "NA")
  stopifnot(file.exists(path), file.info(path)$size > 0)
  invisible(path)
}

metadata <- read_tsv(PATHS[["metadata"]], show_col_types = FALSE, na = "") %>%
  filter(sample_id %in% SAMPLES) %>%
  mutate(
    site = recode(as.character(site), NB = "NA"),
    site = factor(site, levels = c("AZ", "SG", "LC", "NA"))
  )

surface <- metadata %>%
  filter(depth_label == "5_10cm") %>%
  arrange(match(sample_id, SAMPLES))
stopifnot(nrow(surface) == 11L, identical(as.integer(table(surface$site)), c(2L, 3L, 3L, 3L)))

message("Reconstructing the fractional KO abundance matrix...")
counts <- fread(PATHS[["counts"]], select = c("GeneName", SAMPLES))
lengths <- fread(PATHS[["lengths"]], header = FALSE, col.names = c("GeneName", "Length"))
annotation <- fread(PATHS[["annotation"]], select = c("GeneName", "KO"))

gene_rpk <- merge(counts, lengths, by = "GeneName", all.x = TRUE)
gene_rpk[, total_reads := rowSums(.SD, na.rm = TRUE), .SDcols = SAMPLES]
gene_rpk <- gene_rpk[total_reads > 0 & is.finite(Length) & Length > 0]
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

run_test <- function(sample_ids, transformation, scenario, combination = NA_character_) {
  meta <- surface %>% filter(sample_id %in% sample_ids) %>% arrange(match(sample_id, sample_ids))
  group <- droplevels(meta$site)
  x <- t(as.matrix(ko_rpk[, ..sample_ids]))
  storage.mode(x) <- "double"
  rownames(x) <- sample_ids
  x <- x[, colSums(x, na.rm = TRUE) > 0, drop = FALSE]

  if (transformation == "Hellinger-Bray") {
    transformed <- decostand(x, method = "hellinger")
    distance <- vegdist(transformed, method = "bray")
    rank <- NA_integer_
  } else {
    rank <- as.integer(min(3L, nrow(x) - 1L, ncol(x) - 1L))
    transformed <- decostand(x, method = "rclr", impute = TRUE, ropt = rank)
    stopifnot(all(is.finite(transformed)))
    distance <- dist(transformed, method = "euclidean")
  }

  model_data <- data.frame(group = group)
  set.seed(42)
  location <- adonis2(distance ~ group, data = model_data, permutations = 999)
  # PERMDISP is not interpretable when every group has only two observations:
  # each pair is equidistant from its group median and the within-group term is
  # degenerate. Retain dispersion only for the primary and AZ-exclusion runs.
  if (!all(table(group) == 2L)) {
    dispersion <- betadisper(distance, group, type = "median", bias.adjust = TRUE)
    set.seed(42)
    dispersion_test <- permutest(dispersion, permutations = 999)
    dispersion_F <- dispersion_test$tab[1, "F"]
    dispersion_p <- dispersion_test$tab[1, "Pr(>F)"]
  } else {
    dispersion_F <- NA_real_
    dispersion_p <- NA_real_
  }

  tibble(
    scenario = scenario,
    combination = combination,
    transformation = transformation,
    n_samples = nrow(meta),
    n_sites = nlevels(group),
    samples = paste(sample_ids, collapse = ","),
    n_features = ncol(x),
    rclr_rank = rank,
    location_R2 = location$R2[1],
    location_F = location$F[1],
    location_p = location$`Pr(>F)`[1],
    dispersion_F = dispersion_F,
    dispersion_p = dispersion_p
  )
}

primary_ids <- surface$sample_id
non_az_ids <- surface %>% filter(site != "AZ") %>% pull(sample_id)

influence <- map_dfr(c("Hellinger-Bray", "rCLR-Aitchison"), function(transformation) {
  bind_rows(
    run_test(primary_ids, transformation, "primary_four_sites"),
    run_test(non_az_ids, transformation, "exclude_AZ_three_sites")
  )
})
write_tsv_checked(influence, "surface_az_exclusion_sensitivity.tsv")

site_choices <- lapply(levels(surface$site), function(site_name) {
  ids <- surface %>% filter(site == site_name) %>% pull(sample_id)
  if (length(ids) == 2L) list(ids) else combn(ids, 2L, simplify = FALSE)
})
names(site_choices) <- levels(surface$site)
grid <- expand.grid(lapply(site_choices, seq_along), KEEP.OUT.ATTRS = FALSE)

balanced <- map_dfr(seq_len(nrow(grid)), function(i) {
  selected <- unlist(map2(site_choices, as.list(grid[i, ]), ~ .x[[.y]]), use.names = FALSE)
  combination <- sprintf("B%02d", i)
  map_dfr(c("Hellinger-Bray", "rCLR-Aitchison"), function(transformation) {
    run_test(selected, transformation, "balanced_two_per_site", combination)
  })
})
stopifnot(nrow(balanced) == 54L)
write_tsv_checked(balanced, "surface_balanced_two_per_site_sensitivity.tsv")

balanced_summary <- balanced %>%
  group_by(transformation) %>%
  summarise(
    n_balanced_combinations = n(),
    location_R2_min = min(location_R2),
    location_R2_median = median(location_R2),
    location_R2_max = max(location_R2),
    location_p_lt_0_05_n = sum(location_p < 0.05),
    location_p_lt_0_05_fraction = mean(location_p < 0.05),
    dispersion_test_status = "not estimated because each balanced site has n = 2",
    .groups = "drop"
  )
write_tsv_checked(balanced_summary, "surface_balanced_two_per_site_summary.tsv")

capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))
message("V19 surface rCLR/AZ sensitivity completed: ", OUT)
