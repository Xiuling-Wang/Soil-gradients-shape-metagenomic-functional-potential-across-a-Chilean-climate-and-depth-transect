# =============================================================================
# Exploratory screening of focal bacterial MAGs for the Chile metagenome paper.
#
# This script deliberately separates candidate selection from interpretation:
# four candidates are chosen by taxonomy, genome quality, and curated marker
# breadth; a fifth is chosen as a cross-sample Dormibacterota representative.
# Group trajectories are interpreted only after those roles are fixed. All
# abundance values are relative to the recovered, genome-size-normalised
# 228-MAG fraction; they are not estimates of absolute soil-community abundance.
#
# Outputs are exploratory source tables and standalone PDF panels. No PNG files
# are generated and no manuscript or submission-package file is modified.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)))
set.seed(42)

PROJECT <- normalizePath(Sys.getenv("CHILE_PROJECT_ROOT", getwd()), mustWork = TRUE)
OUT_ROOT <- file.path(PROJECT, "analysis/V06_bacterial_MAG_exploration_2026-08-18")
TABLE_DIR <- file.path(OUT_ROOT, "source_data")
FIGURE_DIR <- file.path(OUT_ROOT, "figures")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(PROJECT, "script/_theme_sbb.R"))

MAG_ROOT <- file.path(PROJECT, "data/processed/metagenome_output/MAGs")
MARKER_FILE <- file.path(
  PROJECT,
  "analysis/V20_8_candidate_Figure5_amoA_validation_2026-08-10/source_data",
  "Figure5_MAG_marker_assignments.tsv"
)

samples <- c(
  "C01", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10",
  "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18"
)
group_levels <- c("AZ", "SG", "LC", "NA-S", "NA-M", "NA-D")
group_labels <- c(
  AZ = "AZ", SG = "SG", LC = "LC", `NA-S` = "NA 5-10 cm",
  `NA-M` = "NA 20-40 cm", `NA-D` = "NA 60-80 cm"
)
module_labels <- c(
  nitrogen_transformations = "Nitrogen transformations",
  pedogenic_resource_acquisition = "Pedogenic resource acquisition",
  alternative_energy_chemolithotrophy = "Alternative-energy chemolithotrophy"
)
module_colours <- c(
  nitrogen_transformations = "#4F7898",
  pedogenic_resource_acquisition = "#C79A48",
  alternative_energy_chemolithotrophy = "#5C9482"
)

taxonomy <- read_tsv(file.path(MAG_ROOT, "gtdb_taxonomy.tsv"), show_col_types = FALSE) %>%
  rename(MAG = user_genome)
quality <- read_tsv(file.path(MAG_ROOT, "genome_quality.tsv"), show_col_types = FALSE) %>%
  rename(MAG = `Bin Id`, assembly_sample = Sample)
counts <- read_tsv(file.path(MAG_ROOT, "counts_genomes.tsv"), show_col_types = FALSE) %>%
  rename(MAG = Sample)
markers <- read_tsv(MARKER_FILE, show_col_types = FALSE)
metadata <- read_tsv(file.path(PROJECT, "data/metadata/env_2024_meta.tsv"), show_col_types = FALSE) %>%
  filter(sample_id %in% samples) %>%
  mutate(
    group = case_when(
      site != "NB" ~ site,
      depth_label == "5_10cm" ~ "NA-S",
      depth_label == "20_40cm" ~ "NA-M",
      depth_label == "60_80cm" ~ "NA-D",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = group_levels),
    sample_id = factor(sample_id, levels = samples)
  ) %>%
  select(sample_id, site, depth_label, group)

stopifnot(
  nrow(taxonomy) == 228,
  nrow(quality) == 228,
  nrow(counts) == 228,
  nrow(metadata) == 17,
  setequal(taxonomy$MAG, quality$MAG),
  setequal(taxonomy$MAG, counts$MAG),
  setequal(samples, as.character(metadata$sample_id)),
  all(markers$MAG %in% taxonomy$MAG)
)

# Genome-size normalisation followed by closure within the recovered MAG layer.
mag_abundance <- counts %>%
  pivot_longer(all_of(samples), names_to = "sample_id", values_to = "mapped_reads") %>%
  left_join(quality %>% select(MAG, Genome_Size, assembly_sample), by = "MAG") %>%
  mutate(RPK = mapped_reads / (Genome_Size / 1000)) %>%
  group_by(sample_id) %>%
  mutate(relative_abundance_pct = 100 * RPK / sum(RPK, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(metadata %>% mutate(sample_id = as.character(sample_id)), by = "sample_id") %>%
  mutate(
    sample_id = factor(sample_id, levels = samples),
    group = factor(group, levels = group_levels)
  )

sample_sums <- mag_abundance %>%
  group_by(sample_id) %>%
  summarise(total = sum(relative_abundance_pct), .groups = "drop")
stopifnot(max(abs(sample_sums$total - 100)) < 1e-7)

# Only high/medium-quality bacterial MAGs enter the exploratory universe.
eligible <- taxonomy %>%
  inner_join(
    quality %>%
      select(
        MAG, Completeness, Contamination, Genome_Size, Total_Coding_Sequences,
        Contig_N50, assembly_sample
      ),
    by = "MAG"
  ) %>%
  filter(Domain == "Bacteria", Completeness >= 80, Contamination < 5)

marker_summary <- markers %>%
  filter(MAG %in% eligible$MAG) %>%
  group_by(MAG) %>%
  summarise(
    marker_hits = n(),
    marker_KOs = n_distinct(KO),
    marker_categories = n_distinct(category),
    marker_modules = n_distinct(module),
    .groups = "drop"
  )

abundance_summary <- mag_abundance %>%
  filter(MAG %in% eligible$MAG) %>%
  group_by(MAG, assembly_sample) %>%
  summarise(
    mean_relative_abundance_pct = mean(relative_abundance_pct),
    median_relative_abundance_pct = median(relative_abundance_pct),
    max_relative_abundance_pct = max(relative_abundance_pct),
    assembly_sample_abundance_pct = relative_abundance_pct[sample_id == assembly_sample][1],
    assembly_sample_share_of_cross_sample_sum =
      assembly_sample_abundance_pct / sum(relative_abundance_pct),
    assembly_sample_fold_vs_other_median = assembly_sample_abundance_pct /
      median(relative_abundance_pct[sample_id != assembly_sample]),
    n_samples_ge_0_05_pct = sum(relative_abundance_pct >= 0.05),
    n_nonassembly_samples_ge_0_05_pct =
      sum(relative_abundance_pct >= 0.05 & sample_id != assembly_sample),
    n_samples_ge_10pct_of_max =
      sum(relative_abundance_pct >= 0.1 * max_relative_abundance_pct),
    .groups = "drop"
  )

group_means <- mag_abundance %>%
  filter(MAG %in% eligible$MAG) %>%
  group_by(MAG, group) %>%
  summarise(group_mean_pct = mean(relative_abundance_pct), .groups = "drop") %>%
  mutate(group = as.character(group)) %>%
  pivot_wider(names_from = group, values_from = group_mean_pct, names_prefix = "mean_")

surface_gradient <- mag_abundance %>%
  filter(MAG %in% eligible$MAG, group %in% c("AZ", "SG", "LC", "NA-S")) %>%
  group_by(MAG, group) %>%
  summarise(group_mean_pct = mean(relative_abundance_pct), .groups = "drop") %>%
  mutate(surface_order = match(as.character(group), c("AZ", "SG", "LC", "NA-S"))) %>%
  group_by(MAG) %>%
  summarise(
    surface_gradient_spearman_rho = cor(
      surface_order, group_mean_pct, method = "spearman", use = "complete.obs"
    ),
    .groups = "drop"
  )

candidate_universe <- eligible %>%
  left_join(marker_summary, by = "MAG") %>%
  left_join(abundance_summary, by = c("MAG", "assembly_sample")) %>%
  left_join(group_means, by = "MAG") %>%
  left_join(surface_gradient, by = "MAG") %>%
  mutate(
    across(
      c(marker_hits, marker_KOs, marker_categories, marker_modules),
      ~ replace_na(.x, 0L)
    ),
    quality_score = Completeness - 5 * Contamination,
    assembly_source_flag = case_when(
      assembly_sample_share_of_cross_sample_sum >= 0.35 ~ "high",
      assembly_sample_share_of_cross_sample_sum >= 0.20 ~ "moderate",
      TRUE ~ "low"
    ),
    taxon_label = coalesce(genus, family, order, class, phylum)
  )

# Candidate selection rules are role-based and do not use abundance gradients.
pick_one <- function(data) {
  if (nrow(data) == 0) stop("A predefined candidate role returned no eligible MAG.")
  data %>% slice(1) %>% pull(MAG)
}

mag_bridge <- candidate_universe %>%
  filter(genus == "Bradyrhizobium") %>%
  arrange(desc(marker_categories), desc(quality_score), MAG) %>%
  pick_one()

mag_actino <- candidate_universe %>%
  filter(phylum == "Actinomycetota", MAG != mag_bridge) %>%
  arrange(desc(marker_categories), desc(marker_modules), desc(quality_score), MAG) %>%
  pick_one()

mag_pseudo <- candidate_universe %>%
  filter(phylum == "Pseudomonadota", MAG != mag_bridge) %>%
  arrange(desc(marker_categories), desc(marker_modules), desc(quality_score), MAG) %>%
  pick_one()

mag_dormi_breadth <- candidate_universe %>%
  filter(phylum == "Dormibacterota") %>%
  arrange(desc(marker_categories), desc(marker_modules), desc(quality_score), MAG) %>%
  pick_one()

mag_dormi_comparator <- candidate_universe %>%
  filter(
    phylum == "Dormibacterota", MAG != mag_dormi_breadth,
    assembly_sample_share_of_cross_sample_sum < 0.20,
    n_nonassembly_samples_ge_0_05_pct >= 12
  ) %>%
  arrange(desc(mean_relative_abundance_pct), desc(quality_score), MAG) %>%
  pick_one()

selection_rules <- tribble(
  ~MAG, ~selection_role, ~selection_rule,
  mag_bridge,
  "Cross-project taxonomic bridge",
  "Best eligible Bradyrhizobium MAG by curated marker breadth, then genome quality",
  mag_actino,
  "Marker-rich Actinomycetota",
  "Highest curated marker-category breadth among eligible Actinomycetota",
  mag_pseudo,
  "High-quality functional comparator",
  "Highest curated marker-category breadth among non-Bradyrhizobium Pseudomonadota",
  mag_dormi_breadth,
  "Dormibacterota functional candidate",
  "Highest marker breadth among eligible Dormibacterota, with quality as tie-breaker",
  mag_dormi_comparator,
  "Dormibacterota cross-sample representative",
  "Most abundant remaining Dormibacterota with low assembly-source share and support in at least 12 non-source samples"
)

shortlist <- selection_rules %>%
  left_join(candidate_universe, by = "MAG") %>%
  mutate(
    candidate_order = row_number(),
    candidate_label = paste0(MAG, " | ", taxon_label),
    candidate_label = factor(candidate_label, levels = candidate_label)
  )

stopifnot(
  nrow(shortlist) == 5,
  n_distinct(shortlist$MAG) == 5,
  all(shortlist$MAG %in% candidate_universe$MAG)
)

shortlist_abundance <- mag_abundance %>%
  filter(MAG %in% shortlist$MAG) %>%
  left_join(
    shortlist %>% select(MAG, candidate_order, candidate_label, selection_role),
    by = "MAG"
  ) %>%
  arrange(candidate_order, sample_id)

shortlist_markers <- markers %>%
  filter(MAG %in% shortlist$MAG) %>%
  left_join(shortlist %>% select(MAG, candidate_order), by = "MAG") %>%
  mutate(
    contig = str_remove(gene_id, "_[0-9]+$"),
    orf_index = as.integer(str_extract(gene_id, "[0-9]+$"))
  ) %>%
  arrange(candidate_order, module, category, contig, orf_index)

neighbourhoods <- shortlist_markers %>%
  group_by(MAG, module, category, contig) %>%
  summarise(
    n_marker_genes = n(),
    n_distinct_KOs = n_distinct(KO),
    KOs = paste(sort(unique(KO)), collapse = ";"),
    symbols = paste(sort(unique(symbol)), collapse = ";"),
    first_orf = min(orf_index),
    last_orf = max(orf_index),
    orf_span = last_orf - first_orf,
    compact_multigene_neighbourhood = n_distinct_KOs >= 2 & orf_span <= 15,
    .groups = "drop"
  ) %>%
  arrange(MAG, desc(compact_multigene_neighbourhood), module, category, contig)

# These are explicit marker sets, not claims of complete metabolic pathways.
marker_sets <- tribble(
  ~marker_set, ~KO,
  "CBB core markers", "K00855",
  "CBB core markers", "K01601",
  "CBB core markers", "K01602",
  "Urease structural genes", "K01428",
  "Urease structural genes", "K01429",
  "Urease structural genes", "K01430",
  "PQQ biosynthesis B-E", "K06136",
  "PQQ biosynthesis B-E", "K06137",
  "PQQ biosynthesis B-E", "K06138",
  "PQQ biosynthesis B-E", "K06139",
  "Sox carrier AXYZ", "K17222",
  "Sox carrier AXYZ", "K17223",
  "Sox carrier AXYZ", "K17226",
  "Sox carrier AXYZ", "K17227",
  "Cytosolic DNRA pair", "K00362",
  "Cytosolic DNRA pair", "K00363",
  "Hydrogenase homolog pair", "K06281",
  "Hydrogenase homolog pair", "K06282",
  "APS-reductase homolog pair", "K00394",
  "APS-reductase homolog pair", "K00395",
  "Catecholate synthesis ABCE", "K00216",
  "Catecholate synthesis ABCE", "K01252",
  "Catecholate synthesis ABCE", "K02361",
  "Catecholate synthesis ABCE", "K02363"
)

marker_set_summary <- expand_grid(
  MAG = shortlist$MAG,
  marker_set = unique(marker_sets$marker_set)
) %>%
  left_join(marker_sets %>% count(marker_set, name = "expected_KOs"), by = "marker_set") %>%
  left_join(
    shortlist_markers %>%
      inner_join(marker_sets, by = "KO") %>%
      group_by(MAG, marker_set) %>%
      summarise(
        detected_KOs = n_distinct(KO),
        detected_symbols = paste(sort(unique(symbol)), collapse = ";"),
        n_contigs = n_distinct(contig),
        all_detected_on_one_contig = n_contigs == 1,
        .groups = "drop"
      ),
    by = c("MAG", "marker_set")
  ) %>%
  mutate(
    detected_KOs = replace_na(detected_KOs, 0L),
    detected_symbols = replace_na(detected_symbols, ""),
    marker_set_fraction = detected_KOs / expected_KOs,
    complete_marker_set = detected_KOs == expected_KOs,
    all_detected_on_one_contig = replace_na(all_detected_on_one_contig, FALSE)
  ) %>%
  arrange(MAG, desc(marker_set_fraction), marker_set)

marker_matrix <- shortlist_markers %>%
  count(MAG, module, category, name = "marker_hits") %>%
  left_join(shortlist %>% select(MAG, candidate_order, candidate_label), by = "MAG") %>%
  mutate(
    module = factor(module, levels = names(module_labels)),
    category = fct_reorder(category, as.integer(module), .fun = min),
    candidate_label = factor(
      as.character(candidate_label),
      levels = rev(as.character(shortlist$candidate_label))
    )
  )

# Context needed to decide whether a Dormibacterota observation is best framed
# at phylum level or as a focal-genome example.
dormibacterota_context <- candidate_universe %>%
  filter(phylum == "Dormibacterota") %>%
  arrange(desc(mean_relative_abundance_pct), desc(quality_score), MAG)

dormibacterota_phylum_abundance <- mag_abundance %>%
  left_join(taxonomy %>% select(MAG, phylum), by = "MAG") %>%
  filter(phylum == "Dormibacterota") %>%
  group_by(sample_id, group) %>%
  summarise(
    relative_abundance_pct = sum(relative_abundance_pct),
    n_recovered_MAGs = n_distinct(MAG),
    .groups = "drop"
  ) %>%
  arrange(sample_id)

write_tsv(candidate_universe %>% arrange(desc(marker_categories), desc(quality_score), MAG),
          file.path(TABLE_DIR, "bacterial_MAG_candidate_universe.tsv"))
write_tsv(shortlist %>% mutate(candidate_label = as.character(candidate_label)),
          file.path(TABLE_DIR, "bacterial_MAG_candidate_shortlist.tsv"))
write_tsv(selection_rules, file.path(TABLE_DIR, "bacterial_MAG_selection_rules.tsv"))
write_tsv(shortlist_abundance %>% mutate(across(where(is.factor), as.character)),
          file.path(TABLE_DIR, "bacterial_MAG_shortlist_sample_abundance.tsv"))
write_tsv(shortlist_markers, file.path(TABLE_DIR, "bacterial_MAG_shortlist_marker_assignments.tsv"))
write_tsv(neighbourhoods, file.path(TABLE_DIR, "bacterial_MAG_shortlist_marker_neighbourhoods.tsv"))
write_tsv(marker_set_summary, file.path(TABLE_DIR, "bacterial_MAG_shortlist_marker_sets.tsv"))
write_tsv(dormibacterota_context,
          file.path(TABLE_DIR, "Dormibacterota_eligible_MAG_context.tsv"))
write_tsv(
  dormibacterota_phylum_abundance %>% mutate(across(where(is.factor), as.character)),
  file.path(TABLE_DIR, "Dormibacterota_recovered_fraction_by_sample.tsv")
)

# ---- Standalone exploratory figures ---------------------------------------
abundance_plot <- ggplot(
  shortlist_abundance,
  aes(x = group, y = relative_abundance_pct, colour = group)
) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 42),
    size = 2.0, alpha = 0.82
  ) +
  stat_summary(
    aes(group = group), fun = median, geom = "crossbar",
    width = 0.42, linewidth = 0.35, colour = "black"
  ) +
  facet_wrap(vars(candidate_label), ncol = 1, scales = "free_y") +
  scale_colour_manual(values = pal_group, guide = "none") +
  scale_x_discrete(labels = group_labels, drop = FALSE) +
  labs(
    x = NULL,
    y = "Relative abundance within recovered MAG fraction (%)",
    title = "Exploratory sample distributions of five bacterial MAG candidates",
    subtitle = "Points are independent samples; black bars are group medians. Free y-scales are used for readability."
  ) +
  theme_sbb(base_size = 8.2) +
  theme(
    strip.text = element_text(hjust = 0, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.spacing.y = unit(4, "mm")
  )

ggsave(
  file.path(FIGURE_DIR, "bacterial_MAG_shortlist_abundance_profiles.pdf"),
  abundance_plot, width = 7.2, height = 9.2, units = "in", device = cairo_pdf
)

quality_long <- shortlist %>%
  mutate(candidate_label = fct_reorder(candidate_label, candidate_order, .desc = TRUE)) %>%
  transmute(
    candidate_label,
    phylum,
    `Completeness (%)` = Completeness,
    `Contamination (%)` = Contamination,
    `Assembly-source share (%)` = 100 * assembly_sample_share_of_cross_sample_sum,
    `Non-source samples >= 0.05% (n)` = n_nonassembly_samples_ge_0_05_pct
  ) %>%
  pivot_longer(
    -c(candidate_label, phylum),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c(
        "Completeness (%)", "Contamination (%)", "Assembly-source share (%)",
        "Non-source samples >= 0.05% (n)"
      )
    )
  )

quality_plot <- ggplot(quality_long, aes(x = value, y = candidate_label, colour = phylum)) +
  geom_segment(aes(x = 0, xend = value, yend = candidate_label), colour = "grey82", linewidth = 0.45) +
  geom_point(size = 2.6) +
  geom_text(aes(label = number(value, accuracy = 0.1)), nudge_x = 0.8, size = 2.4, colour = "black") +
  facet_wrap(vars(metric), ncol = 2, scales = "free_x") +
  scale_colour_manual(values = phylum_palette, labels = taxon_plotmath_labels) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    x = NULL, y = NULL, colour = "Phylum",
    title = "Genome quality and cross-sample support of bacterial MAG candidates",
    subtitle = "Assembly-source share flags possible recovery-source leverage; it is not a formal exclusion test."
  ) +
  theme_sbb(base_size = 8.2) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(hjust = 0),
    legend.position = "bottom"
  )

ggsave(
  file.path(FIGURE_DIR, "bacterial_MAG_shortlist_quality_and_support.pdf"),
  quality_plot, width = 8.0, height = 5.7, units = "in", device = cairo_pdf
)

marker_plot <- ggplot(
  marker_matrix,
  aes(x = category, y = candidate_label, size = marker_hits, fill = module)
) +
  geom_point(shape = 21, colour = "grey35", stroke = 0.35, alpha = 0.95) +
  scale_size_area(max_size = 8, breaks = c(1, 3, 6, 10), name = "Marker hits") +
  scale_fill_manual(values = module_colours, labels = module_labels, name = "Module") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 16)) +
  labs(
    x = NULL, y = NULL,
    title = "Curated functional-marker profiles of bacterial MAG candidates",
    subtitle = "Marker presence denotes encoded potential; bubble size is the number of retained gene annotations."
  ) +
  theme_sbb(base_size = 8.2) +
  theme(
    axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.30),
    legend.position = "right"
  )

ggsave(
  file.path(FIGURE_DIR, "bacterial_MAG_shortlist_marker_profiles.pdf"),
  marker_plot, width = 10.2, height = 5.6, units = "in", device = cairo_pdf
)

dormibacterota_plot <- ggplot(
  dormibacterota_phylum_abundance,
  aes(x = group, y = relative_abundance_pct, colour = group)
) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 42),
    size = 2.4, alpha = 0.84
  ) +
  stat_summary(
    aes(group = group), fun = median, geom = "crossbar",
    width = 0.42, linewidth = 0.38, colour = "black"
  ) +
  scale_colour_manual(values = pal_group, guide = "none") +
  scale_x_discrete(labels = group_labels, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(
    x = NULL,
    y = "Relative abundance within recovered MAG fraction (%)",
    title = "Recovered Dormibacterota fraction across sites and Nahuelbuta depths",
    subtitle = paste0(
      "Sum across ", unique(dormibacterota_phylum_abundance$n_recovered_MAGs),
      " recovered Dormibacterota MAGs; points are independent samples and black bars are medians."
    )
  ) +
  theme_sbb(base_size = 8.2) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  file.path(FIGURE_DIR, "Dormibacterota_recovered_fraction_context.pdf"),
  dormibacterota_plot, width = 7.2, height = 4.6, units = "in", device = cairo_pdf
)

audit_lines <- c(
  "Bacterial MAG exploratory screening audit",
  paste0("Generated: ", Sys.Date()),
  paste0("Eligible bacterial MAGs (completeness >=80%, contamination <5%): ", nrow(candidate_universe)),
  paste0("Shortlisted MAGs: ", paste(shortlist$MAG, collapse = ", ")),
  paste0("Recovered Dormibacterota MAGs in phylum-level context: ",
         unique(dormibacterota_phylum_abundance$n_recovered_MAGs)),
  "Four roles used taxonomy, quality, and curated marker breadth; the Dormibacterota representative additionally used cross-sample abundance/support, not a site or depth contrast.",
  "Abundance denominator: genome-size-normalised relative abundance within the recovered 228-MAG layer.",
  "No paired-depth test was run; Nahuelbuta depths are displayed as independent biological samples.",
  "Marker sets are annotation checks, not claims of complete pathways, activity, or biogeochemical flux.",
  "No PNG files were generated."
)
write_lines(audit_lines, file.path(OUT_ROOT, "RUN_AUDIT.txt"))

message("Completed bacterial MAG exploration: ", OUT_ROOT)
