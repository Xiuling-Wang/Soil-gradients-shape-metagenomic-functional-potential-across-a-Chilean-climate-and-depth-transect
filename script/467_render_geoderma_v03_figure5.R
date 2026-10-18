#!/usr/bin/env Rscript
# Geoderma V03 derivative of 458_render_v28_figure4.R : former Fig. 4 minus panel A; B->A, C->B; page 17.2 x 23.8 cm (7.60 + 0.50 + 15.70); phylum names roman in legend and x-axis; output Figure_5.pdf.
# V26 derivative of _tmp_455_narrow.R (443_build_v26_figures_no_panel_titles.py): page width 19.0 -> 17.2 cm; panel-C tile 4.05 (label font 2.85); frame slack 0.40; tag margin 26 pt; compact legends.
# V26 derivative of 455_render_v27_figure4.R (443_build_v26_figures_no_panel_titles.py): panel C: top margin so the tag clears the frame; category tile widened 3.30 -> 3.70 so labels sit clear of its left edge.
# V26 derivative of 447_render_v26_figure4.R (443_build_v26_figures_no_panel_titles.py): A/B facet (site/depth) labels 8.2 -> 10.5 pt.
# V26 derivative of 447_render_v26_figure4.R (443_build_v26_figures_no_panel_titles.py): panel C title (redrawn by the wrapper) removed as well.
# V26 derivative of 435_render_v24_figure4.R (443_build_v26_figures_no_panel_titles.py): sources the no-title Fig. 4 base script.
# V24 derivative of 390_render_readability_figure4.R (built by 438_build_v24_figure_refinement.py). Shorter y-axis titles (panel titles already name the 228-MAG reference set).
# Readability-first Figure 4 refinement.
#
# Reuse the frozen V18 scientific construction, then enlarge only the legend,
# bubble, and supporting-label geometry. No data, marker registry, denominator,
# sample order, or statistic is changed.

suppressPackageStartupMessages({
  library(jsonlite)
})

PROJECT <- Sys.getenv("CHILE_PROJECT_ROOT", unset=getwd())
BUILD_ROOT <- Sys.getenv(
  "CHILE_READABILITY_FIGURE4_ROOT",
  file.path(PROJECT, "analysis/readability_first_figure_refinement_2026-08-30")
)
dir.create(file.path(BUILD_ROOT, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BUILD_ROOT, "source_data"), recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  CHILE_V20_8_FIGURE_ROOT = BUILD_ROOT,
  CHILE_FIGURE4_MAPPING_AXIS_TITLE = "Recruited reads (%)",
  CHILE_FIGURE4_MAPPING_PANEL_TITLE = "Read recruitment to the recovered 228-MAG reference set",
  CHILE_FIGURE4_COMPOSITION_AXIS_TITLE = "Relative abundance (%)",
  CHILE_FIGURE4_COMPOSITION_PANEL_TITLE = "Phylum composition within the recovered MAG fraction",
  CHILE_FIGURE4_CARRIER_PANEL_TITLE = "MAGs carrying ≥1 accepted marker ORF"
)

# Evaluate the canonical renderer only through construction of the three plot
# objects. Its top-level quit() separates the retained figure from retired DRAM
# code, so it is a stable stopping point for this style-only wrapper.
source_path <- file.path(PROJECT, "script/448_build_v26_figure4_base.R")
expressions <- parse(source_path)
for (expression in expressions) {
  if (is.call(expression) && identical(expression[[1]], as.name("quit"))) break
  eval(expression, envir = .GlobalEnv)
}

FONT <- "Arial"

# Geoderma taxon typography: phylum names roman (genus/species only are italic).
# Same phyla, same order, same n; only the <i></i> markup of the source labels is dropped.
phylum_axis_labels <- setNames(
  paste0(phylum_denominators$phylum, "<br>(n = ", phylum_denominators$n_MAG, ")"),
  phylum_denominators$phylum
)

p_mapping_readable <- p_mapping +
  geom_text(
    data = mapping %>% filter(as.character(sample_id) == "C08"),
    aes(
      x = sample_id,
      y = mapping_median + 3.0,
      label = sprintf("Median = %.1f%%", mapping_median)
    ),
    inherit.aes = FALSE,
    family = FONT,
    size = 3.65,
    colour = "grey25"
  ) +
  theme(
    text = element_text(family = FONT),
    axis.text.x = element_text(size = 10.5, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10.5),
    strip.text = element_text(size = 10.5, face = "bold"),
    axis.title.y = element_text(size = 11.2, margin = margin(r = 6)),
    plot.title = element_text(size = 11.7, face = "bold", margin = margin(l = 28, b = 3)),
    plot.tag = element_text(size = 15.5, face = "bold"),
    plot.tag.position = c(0.005, 0.990),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.42)
  )

p_composition_readable <- p_composition +
  theme(
    text = element_text(family = FONT),
    axis.text.x = element_text(size = 10.5, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10.5),
    strip.text = element_text(size = 10.5, face = "bold"),
    axis.title.y = element_text(size = 11.4),
    plot.title = element_text(size = 11.7, face = "bold", margin = margin(l = 28, b = 3)),
    plot.tag = element_text(size = 15.5, face = "bold"),
    plot.tag.position = c(0.005, 0.990),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9.8),
    legend.key.width = unit(7.4, "mm"),
    legend.key.height = unit(4.8, "mm"),
    legend.spacing.x = unit(1.6, "mm"),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.42)
  ) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(alpha = 1))) +
  # Geoderma V03: former panel B becomes panel A; phylum legend labels roman (same values/limits/order).
  labs(tag = "A") +
  scale_fill_manual(
    values = phylum_palette, limits = phylum_levels, drop = FALSE,
    labels = function(x) x
  )

recode_carrier_module <- function(x) {
  recode(
    as.character(x),
    `Nitrogen transformations` = "Nitrogen\ntransformations",
    `Catecholate siderophore synthesis` = "Catecholate siderophore\nsynthesis",
    `Pedogenic resource acquisition` = "Mineral-resource\nacquisition",
    `Alternative-energy chemolithotrophy` = "Alternative-energy\nchemolithotrophy"
  )
}
recode_carrier_category <- function(x) {
  display_labels <- recode(
    as.character(x),
    `Hydroxylamine oxidation` = "Hydroxylamine oxidation",
    `Nitrate/nitrite redox` = "Nitrate/nitrite redox",
    `Assimilatory nitrate reduction` = "Assimilatory nitrate reduction",
    `Organic-acid weathering` = "Organic-acid weathering",
    `Catecholate siderophore synthesis` = "Catecholate siderophore synthesis",
    `CBB carbon-fixation markers` = "CBB-associated carbon fixation",
    `Sulfur redox, direction unresolved` = "Sulfur redox (direction unresolved)",
    `Sulfur redox, direction\nunresolved` = "Sulfur redox (direction unresolved)",
    `Hydrogenase homologs` = "Hydrogenase homologues"
  )
  gsub("[[:space:]]*\\n[[:space:]]*", " ", display_labels)
}
carrier_module_levels <- unique(recode_carrier_module(levels(p_carriers$data$module_label)))
carrier_category_levels <- unique(recode_carrier_category(levels(p_carriers$data$category_label)))

p_carriers$data <- p_carriers$data %>%
  mutate(
    module_label = factor(recode_carrier_module(module_label), levels = carrier_module_levels),
    category_label = factor(recode_carrier_category(category_label), levels = carrier_category_levels)
  )

carrier_plot_data <- p_carriers$data %>%
  mutate(phylum_index = as.numeric(phylum))
carrier_category_cells <- carrier_plot_data %>%
  distinct(module_label, category_label)
bubble_frame_cells <- carrier_plot_data %>%
  distinct(module_label) %>%
  mutate(
    frame_xmin = 0.88,
    frame_xmax = length(carrier_phyla) + 0.40
  )
carrier_category_right <- 0.70
# Compact the category tile so labels stay attached to the bubble matrix rather
# than carrying a broad, empty left margin.
carrier_category_width <- 4.05
carrier_category_x <- carrier_category_right - carrier_category_width / 2
carrier_category_text_x <- carrier_category_right - 0.14
carrier_x_min <- carrier_category_right - carrier_category_width

p_carriers_readable <- ggplot(
  carrier_plot_data,
  aes(phylum_index, category_label)
) +
  geom_tile(
    data = carrier_category_cells,
    aes(x = carrier_category_x, y = category_label),
    inherit.aes = FALSE,
    width = carrier_category_width, height = 0.94,
    fill = "#F1F3F5", colour = "#CCD1D5", linewidth = 0.42
  ) +
  geom_text(
    data = carrier_category_cells,
    aes(x = carrier_category_text_x, y = category_label, label = category_label),
    inherit.aes = FALSE,
    family = FONT, size = 2.85, hjust = 1, lineheight = 0.96,
    colour = "black"
  ) +
  # The two grey text columns form one annotation block. The bubble matrix is
  # a separate data block with its own frame and a positive white gutter.
  geom_rect(
    data = bubble_frame_cells,
    aes(
      xmin = frame_xmin, xmax = frame_xmax,
      ymin = -Inf, ymax = Inf
    ),
    inherit.aes = FALSE,
    fill = NA, colour = "grey72", linewidth = 0.32
  ) +
  geom_point(
    data = ~ filter(.x, carrier_n > 0),
    aes(size = carrier_pct, fill = mean_marker_hits_per_1000_ORFs),
    shape = 21, colour = "grey30", stroke = 0.28
  ) +
  facet_grid(
    rows = vars(module_label), scales = "free_y", space = "free_y",
    switch = "y", drop = FALSE
  ) +
  scale_x_continuous(
    breaks = seq_along(carrier_phyla),
    labels = unname(phylum_axis_labels[carrier_phyla]),
    limits = c(carrier_x_min, length(carrier_phyla) + 0.40),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_size_area(
    name = "Carrier prevalence\nwithin phylum (%)",
    max_size = 10.0,
    limits = c(0, 100),
    breaks = c(25, 50, 75, 100),
    guide = guide_legend(order = 1, nrow = 1, byrow = TRUE)
  ) +
  scale_fill_gradientn(
    colours = pal_seq,
    limits = c(0, marker_fill_max),
    oob = scales::squish,
    name = "Mean accepted marker hits\nper 1,000 ORFs",
    guide = guide_colourbar(
      order = 2,
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(26, "mm"),
      barheight = unit(4.0, "mm"),
      ticks = FALSE
    )
  ) +
  labs(
    tag = "B", x = NULL, y = NULL,
    title = NULL,
    subtitle = NULL
  ) +
  theme_sbb(base_size = 7.5) +
  theme(
    text = element_text(family = FONT),
    axis.text.x = ggtext::element_markdown(
      size = 9.6, angle = 55, hjust = 1, lineheight = 0.98
    ),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_line(colour = "grey91", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.placement = "outside",
    strip.switch.pad.grid = unit(0, "mm"),
    strip.background = element_rect(fill = "#F1F3F5", colour = "#CCD1D5", linewidth = 0.42),
    strip.text.y.left = element_text(
      angle = 0, face = "bold", size = 9.8, hjust = 1, lineheight = 0.98
    ),
    plot.title = element_text(size = 11.7, face = "bold", margin = margin(l = 28, b = 3)),
    plot.tag = element_text(size = 15.5, face = "bold"),
    plot.tag.position = "topleft",
    plot.margin = margin(t = 26, r = 2, b = 2, l = 2),
    plot.tag.location = "margin",
    legend.title = element_text(size = 9.8, face = "bold"),
    legend.text = element_text(size = 9.2),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.key.height = unit(5.5, "mm"),
    legend.key.width = unit(5.5, "mm"),
    legend.spacing.x = unit(2.5, "mm"),
    legend.margin = margin(t = 3, r = 0, b = 0, l = 0),
    panel.spacing.y = unit(1.2, "mm")
  )

# Geoderma V03: two-panel page (former B -> A, former C -> B). Width kept at 17.2 cm (<= 190 mm);
# panel heights 7.60 + 0.50 gap + 15.70 = 23.80 cm (<= 240 mm). Font sizes are those of the V28 source.
panel_heights_cm <- c(A = 7.60, gap_AB = 0.50, B = 15.70)
native_panel_dir <- file.path(BUILD_ROOT, "native_panels")
dir.create(native_panel_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path(native_panel_dir, "Figure_5_A.pdf"), p_composition_readable,
  width = 17.2, height = panel_heights_cm[["A"]], units = "cm",
  bg = "white", device = cairo_pdf
)
ggsave(
  file.path(native_panel_dir, "Figure_5_B.pdf"), p_carriers_readable,
  width = 17.2, height = panel_heights_cm[["B"]], units = "cm",
  bg = "white", device = cairo_pdf
)

figure5_geoderma <- wrap_elements(full = p_composition_readable) /
  plot_spacer() /
  wrap_elements(full = p_carriers_readable) +
  plot_layout(heights = unname(panel_heights_cm))

output <- file.path(BUILD_ROOT, "figures/Figure_5.pdf")
ggsave(
  output,
  figure5_geoderma,
  width = 17.2,
  height = sum(panel_heights_cm),
  units = "cm",
  bg = "white",
  device = cairo_pdf
)

pt_per_cm <- 72 / 2.54
page_height_cm <- sum(panel_heights_cm)
panel_b_bottom <- 0
panel_b_top <- panel_heights_cm[["B"]]
panel_a_bottom <- panel_b_top + panel_heights_cm[["gap_AB"]]
panel_a_top <- page_height_cm
panel_record <- function(id, bottom_cm, top_cm) {
  list(
    id = id,
    role = "panel",
    source_width_pt = 17.2 * pt_per_cm,
    source_height_pt = (top_cm - bottom_cm) * pt_per_cm,
    placed_box_pt = c(0, bottom_cm, 17.2, top_cm) * pt_per_cm,
    ink_box_pt = c(0.15, bottom_cm + 0.12, 17.05, top_cm - 0.12) * pt_per_cm
  )
}

write_json(
  list(
    figure_id = "geoderma_v03_Figure_5_composition_and_carriers",
    scientific_source = "script/448_build_v26_figure4_base.R",
    style_source = "script/458_render_v28_figure4.R",
    data_or_statistical_change = FALSE,
    page_size_cm = c(width = 17.2, height = sum(panel_heights_cm)),
    source_panels_locked_before_composition = TRUE,
    uses_outer_boxes_for_alignment = FALSE,
    aspect_tolerance_pct = 0.05,
    native_panel_sizes_cm = list(
      A = c(width = 17.2, height = panel_heights_cm[["A"]]),
      B = c(width = 17.2, height = panel_heights_cm[["B"]])
    ),
    panel_gaps_cm = c(A_B = panel_heights_cm[["gap_AB"]]),
    panels = list(
      panel_record("A", panel_a_bottom, panel_a_top),
      panel_record("B", panel_b_bottom, panel_b_top)
    ),
    gaps = list(
      list(panel_a = "A", panel_b = "B", axis = "y", kind = "vertical", min_pt = 10, max_pt = 35)
    ),
    changes = c(
      "former Fig. 4A (recruited reads %) moved to Fig. S3",
      "former Fig. 4B relabelled A, former Fig. 4C relabelled B",
      "panel heights 8.00/16.00 -> 7.60/15.70 cm, gap 0.55 -> 0.50 cm; page 17.2 x 23.8 cm",
      "phylum names in the composition legend and carrier x-axis set roman (Geoderma taxon typography)"
    )
  ),
  file.path(BUILD_ROOT, "geometry_ledger_Figure_5.json"),
  pretty = TRUE, auto_unbox = TRUE
)

stopifnot(file.exists(output), file.info(output)$size > 0)
message("Geoderma V03 Figure 5 saved -> ", output)
