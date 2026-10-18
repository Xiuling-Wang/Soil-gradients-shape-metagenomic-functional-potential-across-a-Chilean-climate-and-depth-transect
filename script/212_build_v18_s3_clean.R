# =============================================================================
# V18 Figure S3: phylogenetic placement and genome traits. The former DRAM
# panel remains in Figure S4 and is not duplicated here.
#
# This script promotes the previously reviewed Ayala-style candidate into a
# non-destructive submission-version figure. All panels use real project data
# and describe only the recovered MAG fraction of the metagenomes.
#
# Page 1
#   A  Bacterial phylogeny with four genome-level tracks
#   B  Genome size versus GC content, with marginal densities
#   C  Mean DRAM pathway completeness by major bacterial phylum
# Page 2
#   D  Detailed archaeal tree, abundance profile, and amoA/pmoA annotation
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(cowplot)
  library(dplyr)
  library(ggnewscale)
  library(ggplot2)
  library(ggtext)
  library(ggtree)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(tidyr)
  library(treeio)
})

set.seed(42)

PROJECT <- Sys.getenv("CHILE_PROJECT_ROOT", unset=getwd())
MAG_DIR <- file.path(PROJECT, "data/processed/metagenome_output/MAGs")
META_FILE <- file.path(PROJECT, "data/metadata/env_2024_meta.tsv")
OUT_DIR <- Sys.getenv(
  "CHILE_FIGURE_S3_SUPPORT_DIR",
  file.path(PROJECT, "analysis/V18_figure_caption_cleanup_2026-08-19/source_data/Figure_S3_support")
)
OUT_PDF <- Sys.getenv(
  "CHILE_FIGURE_S3_OUT",
  file.path(PROJECT, "analysis/V18_figure_caption_cleanup_2026-08-19/figures/Figure_S3.pdf")
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

source(file.path(PROJECT, "script/_theme_sbb.R"))

# Shared within-figure evidence palettes. Completeness starts near white at the
# retained 50% threshold and gains contrast only toward high completeness.
completeness_colours <- c(
  "#FAFCFB", "#E5F1ED", "#B9D8CF", "#70AA9A", "#286B5E"
)
pmo_palette <- c(absent = "#F1F3F5", present = "#4B5563")

# ---- 1. Input data and contracts --------------------------------------------
bac_tree <- ape::ladderize(ape::read.tree(file.path(MAG_DIR, "tree/gtdbtk.bac120.nwk")))
ar_tree <- ape::ladderize(ape::read.tree(file.path(MAG_DIR, "tree/gtdbtk.ar53.nwk")))

taxonomy <- read_tsv(
  file.path(MAG_DIR, "gtdb_taxonomy.tsv"), show_col_types = FALSE
) %>%
  rename(MAG = user_genome)

quality <- read_tsv(
  file.path(MAG_DIR, "genome_quality.tsv"), show_col_types = FALSE
) %>%
  rename(MAG = `Bin Id`) %>%
  transmute(
    MAG,
    Completeness,
    Contamination,
    Genome_size_Mb = Genome_Size / 1e6,
    Genome_Size,
    GC_percent = if_else(GC_Content <= 1, GC_Content * 100, GC_Content),
    Coding_density_percent = if_else(Coding_Density <= 1, Coding_Density * 100, Coding_Density)
  )

counts <- read_tsv(
  file.path(MAG_DIR, "counts_genomes.tsv"), show_col_types = FALSE
) %>%
  rename(MAG = Sample)

env <- read_tsv(META_FILE, show_col_types = FALSE) %>%
  select(sample_id, site, depth_label)

product <- read_tsv(
  file.path(MAG_DIR, "dram/distil/product.tsv"), show_col_types = FALSE
) %>%
  rename(MAG = genome)

stopifnot(
  nrow(taxonomy) == 228L,
  nrow(quality) == 228L,
  nrow(product) == 228L,
  length(bac_tree$tip.label) == 221L,
  length(ar_tree$tip.label) == 7L,
  setequal(bac_tree$tip.label, taxonomy$MAG[taxonomy$Domain == "Bacteria"]),
  setequal(ar_tree$tip.label, taxonomy$MAG[taxonomy$Domain == "Archaea"]),
  length(names(counts)[-1]) == 17L,
  all(names(counts)[-1] %in% env$sample_id),
  !"C02" %in% names(counts)[-1]
)

# ---- 2. MAG relative abundance ---------------------------------------------
counts_long <- counts %>%
  pivot_longer(-MAG, names_to = "sample_id", values_to = "reads") %>%
  left_join(quality %>% select(MAG, Genome_Size), by = "MAG") %>%
  mutate(RPK = reads / (Genome_Size / 1000)) %>%
  group_by(sample_id) %>%
  mutate(rel_abundance = RPK / sum(RPK, na.rm = TRUE)) %>%
  ungroup()

mean_abundance <- counts_long %>%
  group_by(MAG) %>%
  summarise(
    mean_rel_abundance = mean(rel_abundance, na.rm = TRUE),
    mean_per_mille = mean_rel_abundance * 1000,
    log10_mean_per_mille = log10(mean_per_mille + 0.001),
    .groups = "drop"
  )

site_order <- c("AZ", "SG", "LC", "NB")
depth_order <- c("5_10cm", "20_40cm", "60_80cm")
group_levels <- c("AZ", "SG", "LC", "NB_5_10cm", "NB_20_40cm", "NB_60_80cm")

group_abundance <- counts_long %>%
  left_join(
    env %>%
      mutate(
        site = factor(site, levels = site_order),
        depth_label = factor(depth_label, levels = depth_order)
      ),
    by = "sample_id"
  ) %>%
  filter(!is.na(site)) %>%
  mutate(
    group = if_else(
      as.character(site) != "NB",
      as.character(site),
      paste0("NB_", as.character(depth_label))
    ),
    group = factor(group, levels = group_levels)
  ) %>%
  group_by(MAG, group) %>%
  summarise(mean_rel = mean(rel_abundance, na.rm = TRUE), .groups = "drop") %>%
  complete(MAG, group = factor(group_levels, levels = group_levels), fill = list(mean_rel = 0)) %>%
  mutate(log10_per_mille = log10(mean_rel * 1000 + 0.001))

# ---- 3. Integrated MAG table ------------------------------------------------
mag_metrics <- taxonomy %>%
  left_join(quality, by = "MAG") %>%
  left_join(mean_abundance, by = "MAG") %>%
  mutate(
    quality_tier = if_else(
      Completeness >= 90 & Contamination < 5,
      "HC/LC",
      "Other retained"
    ),
    phylum_plot = if_else(phylum %in% names(phylum_palette), phylum, "Other")
  )

write_tsv(
  mag_metrics %>%
    select(
      MAG, Domain, phylum, class, order, family, genus,
      Completeness, Contamination, Genome_size_Mb, GC_percent,
      Coding_density_percent, quality_tier,
      mean_rel_abundance, mean_per_mille, log10_mean_per_mille
    ),
  file.path(OUT_DIR, "Figure_S3_MAG_tree_tracks.tsv")
)

# ---- 4. DRAM category summary ----------------------------------------------
classify_pathway <- function(x) {
  case_when(
    str_detect(x, "^CAZy:") ~ "CAZy",
    str_detect(x, "^Complex") ~ "Energy / ETC",
    str_detect(x, "^Methanogenesis and methanotrophy") ~ "Methane metabolism",
    str_detect(x, "^Nitrogen metabolism") ~ "N cycle",
    str_detect(x, "^Other Reductases") ~ "Other reductases",
    str_detect(x, "^Photosynthesis") ~ "Photosynthesis",
    str_detect(x, "^SCFA and alcohol") ~ "SCFA / alcohol",
    str_detect(x, "^Sulfur metabolism") ~ "S cycle",
    TRUE ~ "C metabolism"
  )
}

category_levels <- c(
  "C metabolism", "Energy / ETC", "CAZy", "N cycle", "S cycle",
  "Methane metabolism", "SCFA / alcohol", "Other reductases", "Photosynthesis"
)

pathway_meta <- tibble(pathway = names(product)[-1]) %>%
  mutate(
    category = classify_pathway(pathway),
    category = factor(category, levels = category_levels)
  )

major_bacterial_phyla <- taxonomy %>%
  filter(Domain == "Bacteria") %>%
  count(phylum, name = "n_MAG") %>%
  filter(n_MAG >= 4) %>%
  arrange(desc(n_MAG))

dram_mag_category <- product %>%
  mutate(across(-MAG, as.numeric)) %>%
  pivot_longer(-MAG, names_to = "pathway", values_to = "pathway_completeness") %>%
  left_join(pathway_meta, by = "pathway") %>%
  inner_join(
    taxonomy %>%
      filter(Domain == "Bacteria", phylum %in% major_bacterial_phyla$phylum) %>%
      select(MAG, phylum),
    by = "MAG"
  ) %>%
  group_by(MAG, phylum, category) %>%
  summarise(
    mean_category_completeness = mean(pathway_completeness, na.rm = TRUE),
    n_pathways = n_distinct(pathway),
    .groups = "drop"
  ) %>%
  mutate(mean_category_completeness = replace_na(mean_category_completeness, 0))

dram_phylum_category <- dram_mag_category %>%
  group_by(phylum, category) %>%
  summarise(
    mean_completeness = mean(mean_category_completeness),
    n_MAG = n_distinct(MAG),
    n_pathways = unique(n_pathways),
    .groups = "drop"
  ) %>%
  mutate(
    phylum = factor(phylum, levels = major_bacterial_phyla$phylum),
    category = factor(category, levels = rev(category_levels))
  )

write_tsv(
  dram_phylum_category %>% mutate(across(where(is.factor), as.character)),
  file.path(OUT_DIR, "Figure_S3_DRAM_category_by_bacterial_phylum.tsv")
)

# ---- 5. Panel A: bacterial phylogeny and genome tracks ----------------------
tree_track_data <- mag_metrics %>%
  filter(Domain == "Bacteria") %>%
  transmute(
    bin = MAG,
    Completeness,
    Contamination,
    Genome_size_Mb,
    mean_per_mille,
    log10_mean_per_mille
  )

tree_meta <- mag_metrics %>%
  filter(Domain == "Bacteria") %>%
  transmute(label = MAG, phylum_plot, quality_tier)

abundance_limits <- range(tree_track_data$log10_mean_per_mille, finite = TRUE)
genome_size_limits <- range(tree_track_data$Genome_size_Mb, finite = TRUE)
contamination_limits <- c(0, max(6, ceiling(max(tree_track_data$Contamination, na.rm = TRUE))))

make_track_matrix <- function(data, value, label) {
  out <- data %>%
    select(bin, value = all_of(value)) %>%
    tibble::column_to_rownames("bin") %>%
    as.data.frame()
  names(out) <- label
  out
}

abundance_matrix <- make_track_matrix(
  tree_track_data, "log10_mean_per_mille", "Mean abundance"
)
completeness_matrix <- make_track_matrix(
  tree_track_data, "Completeness", "Completeness"
)
contamination_matrix <- make_track_matrix(
  tree_track_data, "Contamination", "Contamination"
)
genome_size_matrix <- make_track_matrix(
  tree_track_data, "Genome_size_Mb", "Genome size"
)

p_tree_base <- ggtree(
  bac_tree, layout = "fan", open.angle = 18, size = 0.16, color = "grey55"
) %<+% tree_meta +
  geom_tippoint(
    aes(fill = phylum_plot, alpha = quality_tier),
    shape = 21, size = 1.45, stroke = 0.12, color = "white"
  ) +
  scale_fill_manual(
    name = "Phylum",
    values = phylum_palette,
    breaks = c(major_bacterial_phyla$phylum, "Other"),
    labels = function(x) taxon_plotmath_labels(x),
    guide = guide_legend(ncol = 4, override.aes = list(size = 2.4, alpha = 1))
  ) +
  scale_alpha_manual(
    name = "Operational stratum",
    values = c("HC/LC" = 1, "Other retained" = 0.42),
    breaks = c("HC/LC", "Other retained"),
    labels = c(
      "HC/LC" = "HC/LC (\u226590%/<5%)",
      "Other retained" = "Other retained (\u226550%/<10%)"
    ),
    guide = guide_legend(nrow = 1, override.aes = list(size = 2.6))
  )

p_tree_abundance <- p_tree_base + ggnewscale::new_scale_fill()
p_tree_abundance <- gheatmap(
  p_tree_abundance, abundance_matrix,
  offset = 0.035, width = 0.105, color = SBB_HEATMAP_BORDER, colnames = FALSE
) +
  scale_fill_gradientn(
    colours = pal_seq,
    name = expression(atop("Mean abundance", log[10] * " per-mille")),
    limits = abundance_limits,
    breaks = pretty(abundance_limits, n = 3),
    guide = "none"
  )

p_tree_completeness <- p_tree_abundance + ggnewscale::new_scale_fill()
p_tree_completeness <- gheatmap(
  p_tree_completeness, completeness_matrix,
  offset = 0.155, width = 0.105, color = SBB_HEATMAP_BORDER, colnames = FALSE
) +
  scale_fill_gradientn(
    colours = completeness_colours,
    values = scales::rescale(c(50, 65, 80, 90, 100), from = c(50, 100)),
    name = "Completeness (%)", limits = c(50, 100), oob = squish,
    guide = "none"
  )

p_tree_contamination <- p_tree_completeness + ggnewscale::new_scale_fill()
p_tree_contamination <- gheatmap(
  p_tree_contamination, contamination_matrix,
  offset = 0.275, width = 0.105, color = SBB_HEATMAP_BORDER, colnames = FALSE
) +
  scale_fill_gradient(
    low = "white", high = "grey38",
    name = "Contamination (%)", limits = contamination_limits, oob = squish,
    guide = "none"
  )

p_tree_genome_size <- p_tree_contamination + ggnewscale::new_scale_fill()
p_tree <- gheatmap(
  p_tree_genome_size, genome_size_matrix,
  offset = 0.395, width = 0.105, color = SBB_HEATMAP_BORDER, colnames = FALSE
) +
  scale_fill_gradient(
    low = "#F2ECF5", high = unname(pal_category[["mauve"]]),
    name = "Genome size (Mb)", limits = genome_size_limits, oob = squish,
    guide = "none"
  ) +
  labs(
    tag = "A",
    title = "Recovered bacterial MAGs",
    subtitle = NULL
  ) +
  theme_void(base_family = SBB_FONT) +
  theme(
    text = element_text(family = SBB_FONT, color = "black"),
    plot.title = element_text(size = 9.2, face = "bold", hjust = 0, margin = margin(l = 22)),
    plot.subtitle = element_text(size = 6.3, color = "grey25", hjust = 0, lineheight = 1.15),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.015, 0.985),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.title = element_text(size = 5.8, face = "bold"),
    legend.text = element_text(size = 5.5),
    legend.spacing.y = unit(0.02, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(7, 4, 5, 4)
  )

# ---- 6. Panel B: observed genome trait space -------------------------------
trait_data <- mag_metrics %>%
  filter(Domain == "Bacteria") %>%
  mutate(
    phylum_plot = factor(
      phylum_plot,
      levels = c(major_bacterial_phyla$phylum, "Other")
    )
  )

p_trait_top <- ggplot(trait_data, aes(x = Genome_size_Mb)) +
  geom_density(fill = "grey78", color = "grey35", linewidth = 0.35, alpha = 0.75) +
  labs(tag = "B", title = "Bacterial MAG genome traits") +
  theme_void(base_family = SBB_FONT) +
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.015, 0.97),
    plot.title = element_text(size = 9, face = "bold", margin = margin(l = 22, b = 2)),
    plot.margin = margin(4, 0, 0, 0)
  )

p_trait <- ggplot(
  trait_data,
  aes(x = Genome_size_Mb, y = GC_percent, fill = phylum_plot)
) +
  geom_point(shape = 21, color = "white", stroke = 0.18, alpha = 0.80, size = 2.2) +
  scale_fill_manual(values = phylum_palette, guide = "none") +
  labs(x = "Genome size (Mb)", y = "GC content (%)") +
  theme_sbb(base_size = 7.2) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 2, 2)
  )

p_trait_right <- ggplot(trait_data, aes(y = GC_percent)) +
  geom_density(fill = "grey78", color = "grey35", linewidth = 0.35, alpha = 0.75) +
  coord_flip() +
  theme_void(base_family = SBB_FONT) +
  theme(plot.margin = margin(0, 2, 2, 0))

trait_design <- c(
  patchwork::area(t = 1, l = 1, b = 1, r = 1),
  patchwork::area(t = 1, l = 2, b = 1, r = 2),
  patchwork::area(t = 2, l = 1, b = 2, r = 1),
  patchwork::area(t = 2, l = 2, b = 2, r = 2)
)
p_trait_block <- wrap_plots(
  p_trait_top, plot_spacer(), p_trait, p_trait_right,
  design = trait_design,
  widths = c(4.4, 1), heights = c(1, 4.2)
)

# ---- 7. Panel C: broad DRAM category summary -------------------------------
phylum_axis_labels <- major_bacterial_phyla %>%
  transmute(phylum, label = paste0(phylum, "\n(n=", n_MAG, ")")) %>%
  tibble::deframe()

p_dram <- ggplot(
  dram_phylum_category %>%
    mutate(category = factor(as.character(category), levels = category_levels)),
  aes(x = category, y = phylum, fill = mean_completeness)
) +
  geom_tile(color = SBB_HEATMAP_BORDER, linewidth = 0.38) +
  scale_fill_gradientn(
    colours = pal_seq,
    limits = c(0, max(dram_phylum_category$mean_completeness, na.rm = TRUE)),
    oob = squish,
    name = "Mean pathway\ncompleteness",
    guide = guide_colorbar(barheight = 2.8, barwidth = 0.35)
  ) +
  scale_x_discrete(
    labels = c(
      "C metabolism" = "C metabolism",
      "Energy / ETC" = "Energy / ETC",
      "CAZy" = "CAZy",
      "N cycle" = "N cycle",
      "S cycle" = "S cycle",
      "Methane metabolism" = "Methane",
      "SCFA / alcohol" = "SCFA / alcohol",
      "Other reductases" = "Other reductases",
      "Photosynthesis" = "Photosynthesis"
    ),
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = rev(major_bacterial_phyla$phylum),
    labels = phylum_axis_labels,
    drop = FALSE
  ) +
  labs(
    tag = "C",
    title = "Phylum-level DRAM context",
    subtitle = "Within-MAG pathway mean; then phylum mean",
    x = NULL, y = NULL
  ) +
  theme_sbb(base_size = 7.0) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 50, hjust = 1, vjust = 1, size = 5.0
    ),
    axis.text.y = element_text(size = 5.9, face = "italic"),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.01, 0.985),
    plot.title = element_text(size = 9, margin = margin(l = 22, b = 2)),
    plot.subtitle = element_text(size = 6.4, margin = margin(b = 3)),
    legend.position = "right",
    plot.margin = margin(4, 4, 7, 2)
  )

# ---- 8. Panel D: detailed archaeal view ------------------------------------
pmo_col <- "Methanogenesis and methanotrophy: methane => methanol, with oxygen (pmo)"
stopifnot(pmo_col %in% names(product))

pmo_data <- product %>%
  transmute(
    MAG,
    pmo_raw = .data[[pmo_col]],
    pmo_present = case_when(
      str_to_lower(as.character(pmo_raw)) == "true" ~ "present",
      str_to_lower(as.character(pmo_raw)) == "false" ~ "absent",
      is.na(pmo_raw) ~ "absent",
      suppressWarnings(as.numeric(as.character(pmo_raw))) > 0 ~ "present",
      TRUE ~ "absent"
    )
  )

ar_tips <- ar_tree$tip.label
ar_meta <- mag_metrics %>%
  filter(MAG %in% ar_tips) %>%
  left_join(pmo_data, by = "MAG") %>%
  mutate(
    lower_taxon = case_when(
      !is.na(genus) & genus != "" &
        !str_detect(genus, "^[A-Z0-9]+(?:-[A-Z0-9]+){2,}$") ~ genus,
      !is.na(class) & class != "" ~ class,
      !is.na(order) & order != "" ~ order,
      TRUE ~ "Unclassified"
    ),
    tip_display = paste0("<b>", MAG, "</b>")
  ) %>%
  transmute(
    label = MAG, MAG, tip_display, phylum, class, lower_taxon, phylum_plot,
    Completeness, pmo_present
  )

archaeal_phylum_levels <- c("Thermoproteota", "Thermoplasmatota")
archaeal_class_levels <- c("Bathyarchaeia", "Nitrososphaeria", "Thermoplasmata")
archaeal_class_palette <- c(
  Bathyarchaeia = "#377EB8",
  Nitrososphaeria = "#4DAF4A",
  Thermoplasmata = "#984EA3"
)
archaeal_phylum_counts <- ar_meta %>%
  count(phylum, name = "n") %>%
  arrange(match(phylum, archaeal_phylum_levels))
archaeal_class_counts <- ar_meta %>%
  count(class, name = "n") %>%
  arrange(match(class, archaeal_class_levels))
stopifnot(
  identical(archaeal_phylum_counts$phylum, archaeal_phylum_levels),
  identical(archaeal_phylum_counts$n, c(6L, 1L)),
  identical(archaeal_class_counts$class, archaeal_class_levels),
  identical(archaeal_class_counts$n, c(4L, 2L, 1L)),
  all(ar_meta$phylum_plot == ar_meta$phylum)
)

write_tsv(
  bind_rows(
    tibble(
      role = "Phylum", level = archaeal_phylum_levels,
      hex = unname(phylum_palette[archaeal_phylum_levels])
    ),
    tibble(
      role = "Class", level = archaeal_class_levels,
      hex = unname(archaeal_class_palette[archaeal_class_levels])
    ),
    tibble(role = "amoA/pmoA", level = names(pmo_palette), hex = unname(pmo_palette)),
    tibble(
      role = "Completeness", level = c("50", "65", "80", "90", "100"),
      hex = completeness_colours
    )
  ),
  file.path(OUT_DIR, "Figure_S3_palette.tsv")
)

write_tsv(
  ar_meta %>%
    mutate(tree_input_order = match(MAG, ar_tips)) %>%
    arrange(tree_input_order) %>%
    select(
      tree_input_order, MAG, phylum, class, lower_taxon,
      Completeness, pmo_present
    ),
  file.path(OUT_DIR, "Figure_S3_archaeal_taxonomy.tsv")
)

ar_abundance_matrix <- group_abundance %>%
  filter(MAG %in% ar_tips) %>%
  select(MAG, group, log10_per_mille) %>%
  pivot_wider(names_from = group, values_from = log10_per_mille) %>%
  tibble::column_to_rownames("MAG") %>%
  as.data.frame()
ar_abundance_matrix <- ar_abundance_matrix[
  ar_tips, group_levels, drop = FALSE
]
colnames(ar_abundance_matrix) <- c("AZ", "SG", "LC", "NA\n5-10", "NA\n20-40", "NA\n60-80")

ar_class_matrix <- ar_meta %>%
  select(MAG, class) %>%
  mutate(class = factor(class, levels = archaeal_class_levels)) %>%
  tibble::column_to_rownames("MAG") %>%
  as.data.frame()
ar_class_matrix <- ar_class_matrix[ar_tips, , drop = FALSE]
colnames(ar_class_matrix) <- "Class"

ar_completeness_matrix <- ar_meta %>%
  select(MAG, Completeness) %>%
  tibble::column_to_rownames("MAG") %>%
  as.data.frame()
ar_completeness_matrix <- ar_completeness_matrix[ar_tips, , drop = FALSE]
colnames(ar_completeness_matrix) <- "Comp."

ar_pmo_matrix <- ar_meta %>%
  select(MAG, pmo_present) %>%
  tibble::column_to_rownames("MAG") %>%
  as.data.frame()
ar_pmo_matrix <- ar_pmo_matrix[ar_tips, , drop = FALSE]
colnames(ar_pmo_matrix) <- "amoA/pmoA"

tree_theme <- theme(
  text = element_text(family = SBB_FONT, color = "black"),
  legend.background = element_blank(),
  legend.key = element_blank(),
  legend.text = element_text(family = SBB_FONT, size = 6.8),
  legend.title = element_text(family = SBB_FONT, size = 7.2, face = "bold"),
  plot.title = element_text(
    family = SBB_FONT, size = 9.5, face = "bold", hjust = 0,
    margin = margin(l = 22)
  ),
  plot.subtitle = element_text(
    family = SBB_FONT, size = 7, color = "grey25",
    margin = margin(t = 4, l = 22)
  ),
  plot.margin = margin(8, 10, 12, 8)
)

p_ar <- ggtree(ar_tree, layout = "rectangular", size = 0.62, color = "grey28") %<+% ar_meta +
  ggtext::geom_richtext(
    data = function(x) dplyr::filter(x, !is.na(tip_display)),
    aes(label = tip_display), hjust = 0, vjust = 0.5, nudge_x = 0.18,
    size = 3.35, color = "grey10", family = SBB_FONT,
    fill = NA, label.color = NA, label.padding = unit(0, "pt"),
    label.r = unit(0, "pt")
  ) +
  geom_tippoint(aes(color = phylum_plot), size = 4.1, alpha = 0.92) +
  scale_color_manual(
    name = "Phylum",
    values = phylum_palette,
    breaks = archaeal_phylum_levels,
    labels = function(x) taxon_plotmath_labels(x),
    guide = guide_legend(
      order = 1,
      override.aes = list(size = 4.1, alpha = 1)
    )
  ) +
  labs(
    tag = "C",
    title = "Recovered archaeal MAGs",
    subtitle = NULL
  ) +
  xlim(0, 8.4) +
  ylim(0.25, 8.75) +
  tree_theme +
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.015, 0.985),
    plot.subtitle = ggtext::element_markdown(
      family = SBB_FONT, size = 7, color = "grey25",
      margin = margin(t = 4, l = 22)
    )
  )

p_ar_c <- gheatmap(
  p_ar, ar_class_matrix,
  offset = 1.45, width = 0.34, color = SBB_HEATMAP_BORDER,
  colnames = TRUE, colnames_position = "top", colnames_angle = 0,
  hjust = 0, colnames_offset_y = 0.30, font.size = 2.7
) +
  scale_fill_manual(
    name = "Class",
    values = archaeal_class_palette,
    breaks = archaeal_class_levels,
    labels = function(x) taxon_plotmath_labels(x),
    guide = guide_legend(order = 2, override.aes = list(size = 3))
  )

p_ar_q <- p_ar_c + ggnewscale::new_scale_fill()
p_ar_q <- gheatmap(
  p_ar_q, ar_completeness_matrix,
  offset = 2.25, width = 0.24, color = SBB_HEATMAP_BORDER,
  colnames = TRUE, colnames_position = "top", colnames_angle = 0,
  hjust = 0, colnames_offset_y = 0.30, font.size = 2.7
) +
  scale_fill_gradientn(
    colours = completeness_colours,
    values = scales::rescale(c(50, 65, 80, 90, 100), from = c(50, 100)),
    name = "Completeness (%)", limits = c(50, 100), oob = squish,
    guide = guide_colorbar(order = 5, barwidth = 0.5, barheight = 4.2)
  )

p_ar_h <- p_ar_q + ggnewscale::new_scale_fill()
p_ar_h <- gheatmap(
  p_ar_h, ar_abundance_matrix,
  offset = 2.80, width = 3.00, color = SBB_HEATMAP_BORDER,
  colnames = TRUE, colnames_position = "top", colnames_angle = 0,
  hjust = 0.5, colnames_offset_y = 0.30, font.size = 2.5
) +
  scale_fill_gradientn(
    colours = pal_seq,
    name = expression(atop("Relative abundance", log[10] * " per-mille")),
    limits = c(-3, 1.5), oob = squish,
    breaks = c(-3, -1.5, 0, 1.5),
    guide = guide_colorbar(order = 4, barwidth = 0.5, barheight = 4.2)
  )

p_ar_f <- p_ar_h + ggnewscale::new_scale_fill()
p_ar_f <- gheatmap(
  p_ar_f, ar_pmo_matrix,
  offset = 5.95, width = 0.30, color = SBB_HEATMAP_BORDER,
  colnames = TRUE, colnames_position = "top", colnames_angle = 0,
  hjust = 0, colnames_offset_y = 0.30, font.size = 2.7
) +
  scale_fill_manual(
    name = expression(italic(amoA) * "/" * italic(pmoA)),
    values = pmo_palette,
    breaks = c("absent", "present"),
    guide = guide_legend(order = 3, override.aes = list(size = 3))
  ) +
  theme(legend.position = "right")

# ---- 9. Assemble vector PDF -------------------------------------------------
# Keep the accepted source-panel proportions. The empty lower-right area is
# deliberate: it prevents the genome-trait panel from being stretched after
# the duplicate DRAM panel was removed.
right_column <- plot_spacer() / p_trait_block / plot_spacer() +
  plot_layout(heights = c(0.10, 0.42, 0.48))
page_1 <- cowplot::plot_grid(
  p_tree, right_column,
  nrow = 1, rel_widths = c(1.42, 1), align = "none"
)

page_1_pdf <- tempfile(fileext = ".pdf")
page_2_pdf <- tempfile(fileext = ".pdf")

ggsave(
  page_1_pdf, page_1,
  width = 20.0, height = 18.0, units = "cm",
  device = cairo_pdf, bg = "white"
)
ggsave(
  page_2_pdf, p_ar_f,
  width = 20.0, height = 17.0, units = "cm",
  device = cairo_pdf, bg = "white"
)

if (file.exists(OUT_PDF)) unlink(OUT_PDF)
merge_status <- system2(
  "pdfunite",
  c(shQuote(page_1_pdf), shQuote(page_2_pdf), shQuote(OUT_PDF))
)
unlink(c(page_1_pdf, page_2_pdf))
if (!identical(merge_status, 0L)) stop("Could not merge integrated PDF pages")

qa_summary <- tibble(
  check = c(
    "total_MAGs", "bacterial_MAGs", "archaeal_MAGs", "operational_HC_LC_MAGs",
    "archaeal_phyla", "archaeal_phylum_counts",
    "archaeal_classes", "archaeal_class_counts",
    "samples", "major_bacterial_phyla_n_ge_4", "DRAM_pathways", "output_pdf"
  ),
  value = c(
    nrow(mag_metrics),
    sum(mag_metrics$Domain == "Bacteria"),
    sum(mag_metrics$Domain == "Archaea"),
    sum(mag_metrics$quality_tier == "HC/LC"),
    nrow(archaeal_phylum_counts),
    paste0(
      archaeal_phylum_counts$phylum, "=", archaeal_phylum_counts$n,
      collapse = ";"
    ),
    nrow(archaeal_class_counts),
    paste0(
      archaeal_class_counts$class, "=", archaeal_class_counts$n,
      collapse = ";"
    ),
    n_distinct(counts_long$sample_id),
    nrow(major_bacterial_phyla),
    ncol(product) - 1L,
    OUT_PDF
  )
)
write_tsv(qa_summary, file.path(OUT_DIR, "Figure_S3_QA_summary.tsv"))

message("Integrated Figure S3 written to: ", OUT_PDF)
