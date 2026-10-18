#!/usr/bin/env Rscript
# Geoderma V03 derivative of 445_render_v26_figure3.R : split into two pages via CHILE_GEODERMA_SPLIT_MODE=figure3 (nitrogen block A only -> Figure_3.pdf) | figure4 (mineral-resource + alternative-energy blocks relabelled A/B -> Figure_4.pdf); heatmap body 12.96 -> 11.4 cm so the page is 18.9 cm wide (<= 190 mm); data, colours, order, z-scores, fonts unchanged; single-panel Fig. 3 drawn without a panel letter.
# V26 derivative of 389_render_readability_figure3.R (443_build_v26_figures_no_panel_titles.py): module titles next to A/B/C removed; title band 1.20 -> 0.60 cm (heatmap viewport unchanged).
# Readability-first Figure 3 refinement.
#
# The three heatmap bodies use one fixed width and one shared right edge. Row
# labels and process labels may extend to the left, while the shared legend is
# drawn independently on the right so it cannot displace panel A.

suppressPackageStartupMessages({
  library(circlize)
  library(ComplexHeatmap)
  library(dplyr)
  library(grid)
  library(jsonlite)
  library(readr)
  library(tidyr)
})
if (!requireNamespace("gridtext", quietly = TRUE)) {
  stop("Package 'gridtext' is required for Arial-only italic gene labels.")
}

options(scipen = 999)
set.seed(42)

PROJECT <- Sys.getenv("CHILE_PROJECT_ROOT", unset=getwd())
BUILD_ROOT <- Sys.getenv(
  "CHILE_READABILITY_FIGURE3_ROOT",
  file.path(PROJECT, "analysis/readability_first_figure_refinement_2026-08-30")
)
FIGURE_DIR <- file.path(BUILD_ROOT, "figures")
SOURCE_DIR <- file.path(BUILD_ROOT, "source_data")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

FONT <- "Arial"
pdfFonts(Arial = pdfFonts("ArialMT")[[1]])
if (capabilities("aqua")) quartzFonts(
  Arial = quartzFont(c("Arial", "Arial Bold", "Arial Italic", "Arial Bold Italic"))
)
SAMPLES <- c("C01", sprintf("C%02d", 3:18))
GROUP_LEVELS <- c("AZ", "SG", "LC", "NA-S", "NA-M", "NA-D")
HEATMAP_BODY_WIDTH_CM <- 11.4
GROUP_COLOURS <- c(
  AZ = "#C25B3F", SG = "#D6A141", LC = "#6FA055",
  `NA-S` = "#8DB5CD", `NA-M` = "#5F88AA", `NA-D` = "#2E5577"
)
PRIMARY_MODULES <- c(
  "nitrogen_transformations",
  "pedogenic_resource_acquisition",
  "alternative_energy_chemolithotrophy"
)
PANEL_LETTERS <- c(
  nitrogen_transformations = "A",
  pedogenic_resource_acquisition = "B",
  alternative_energy_chemolithotrophy = "C"
)
MODULE_TITLES <- c(
  nitrogen_transformations = "Nitrogen-related functional genes",
  pedogenic_resource_acquisition = "Mineral-resource acquisition genes",
  alternative_energy_chemolithotrophy = "Alternative-energy chemolithotrophy genes"
)
# Geoderma V03 split switch. The full three-module data set is still read, filtered and
# frozen to source_data exactly as in 445; only the modules drawn on THIS page differ.
SPLIT_MODE <- Sys.getenv("CHILE_GEODERMA_SPLIT_MODE", "figure3")
if (SPLIT_MODE == "figure3") {
  RENDER_MODULES <- "nitrogen_transformations"
  PANEL_LETTERS <- c(nitrogen_transformations = "A")  # ledger bookkeeping only; not drawn
  DRAW_PANEL_LETTERS <- FALSE  # single-panel figure: no (A) tag on the page
  OUTPUT_NAME <- "Figure_3.pdf"
} else if (SPLIT_MODE == "figure4") {
  RENDER_MODULES <- c("pedogenic_resource_acquisition", "alternative_energy_chemolithotrophy")
  PANEL_LETTERS <- c(pedogenic_resource_acquisition = "A", alternative_energy_chemolithotrophy = "B")
  DRAW_PANEL_LETTERS <- TRUE
  OUTPUT_NAME <- "Figure_4.pdf"
} else {
  stop("CHILE_GEODERMA_SPLIT_MODE must be 'figure3' or 'figure4'")
}

marker_path <- file.path(
  PROJECT,
  "analysis/V42_ecological_story_revision_2026-08-26/source_data",
  "Figure_S4_full_sample_marker_values.tsv"
)
if (!file.exists(marker_path)) stop("Missing frozen sample-level Figure 3 values")

sample_values <- read_tsv(marker_path, show_col_types = FALSE) %>%
  filter(module %in% PRIMARY_MODULES, sample_id %in% SAMPLES) %>%
  mutate(
    category = recode(
      category,
      `Assimilatory nitrate reduction` = "Assimilatory nitrate\nreduction",
      `Hydroxylamine oxidation` = "Hydroxylamine oxidation",
      `Nitrate/nitrite redox` = "Nitrate/nitrite\nredox",
      `Catecholate siderophore synthesis` = "Catecholate siderophore\nsynthesis",
      `EPS and aggregation` = "EPS and aggregation",
      `Organic-acid weathering` = "Organic-acid weathering",
      `CBB carbon-fixation markers` = "CBB carbon fixation",
      `Hydrogenase homologs` = "Hydrogenase homologues",
      `Putative CO oxidation` = "Putative CO oxidation",
      `Sulfur redox, direction unresolved` = "Sulfur redox\n(direction unresolved)"
    ),
    module = factor(module, levels = PRIMARY_MODULES),
    sample_id = factor(sample_id, levels = SAMPLES),
    group = factor(group, levels = GROUP_LEVELS)
  ) %>%
  arrange(module, category, marker, sample_id)

stopifnot(
  n_distinct(sample_values$sample_id) == 17L,
  !anyNA(sample_values$display_z),
  all(table(sample_values$module, sample_values$sample_id) > 0)
)

marker_order <- sample_values %>%
  distinct(module, category, marker, marker_label) %>%
  mutate(marker_order = row_number())
sample_group_key <- sample_values %>%
  distinct(sample_id, group) %>%
  arrange(sample_id)
stopifnot(identical(as.character(sample_group_key$sample_id), SAMPLES))

write_tsv(
  sample_values %>%
    select(module, category, marker, marker_label, sample_id, group, TPM, display_z),
  file.path(SOURCE_DIR, "Figure_3_sample_resolved_marker_values.tsv")
)

z_col_fun <- colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "#F7F7F7", "#B2182B"))

build_heatmap <- function(module_name) {
  marker_info <- marker_order %>%
    filter(module == module_name) %>%
    arrange(marker_order) %>%
    mutate(row_key = make.unique(paste0(marker_label, " [", marker, "]")))

  matrix_data <- sample_values %>%
    filter(module == module_name) %>%
    inner_join(
      marker_info %>% select(category, marker, marker_label, row_key, marker_order),
      by = c("category", "marker", "marker_label")
    ) %>%
    select(row_key, marker_order, sample_id, display_z) %>%
    pivot_wider(names_from = sample_id, values_from = display_z) %>%
    arrange(marker_order) %>%
    select(row_key, all_of(SAMPLES))

  matrix_value <- as.matrix(matrix_data[, SAMPLES])
  rownames(matrix_value) <- matrix_data$row_key
  category_split <- factor(marker_info$category, levels = unique(marker_info$category))
  column_split <- factor(as.character(sample_group_key$group), levels = GROUP_LEVELS)
  row_labels <- sprintf(
    "<i>%s</i> [%s]",
    marker_info$marker_label,
    marker_info$marker
  )

  group_annotation <- HeatmapAnnotation(
    Group = sample_group_key$group,
    col = list(Group = GROUP_COLOURS),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = unit(3.3, "mm"),
    gp = gpar(col = "white", lwd = 0.40)
  )

  # Separate the functional-group and marker identifiers into two visual table
  # columns. The marker column is drawn per heatmap slice so its fine grid stays
  # aligned with the corresponding KO rows.
  gene_ko_annotation <- AnnotationFunction(
    fun = function(index, k = 1, n = 1) {
      n_rows <- length(index)
      y_positions <- unit((n_rows - seq_len(n_rows) + 0.5) / n_rows, "npc")
      grid.rect(
        x = unit(0.5, "npc"), y = unit(0.5, "npc"),
        width = unit(1, "npc"), height = unit(1, "npc"),
        gp = gpar(fill = "#F1F3F5", col = "#CCD1D5", lwd = 0.55)
      )
      grid.rect(
        x = unit(0.5, "npc"), y = y_positions,
        width = unit(1, "npc"), height = unit(1 / n_rows, "npc"),
        gp = gpar(fill = "#F1F3F5", col = "#D8DDE1", lwd = 0.35)
      )
      grid.draw(gridtext::richtext_grob(
        row_labels[index], x = unit(0.97, "npc"), y = y_positions,
        hjust = 1, vjust = 0.5,
        gp = gpar(fontfamily = FONT, fontsize = 9.0, col = "black"),
        box_gp = gpar(col = NA, fill = NA)
      ))
    },
    which = "row", width = unit(30, "mm")
  )

  functional_group_annotation <- AnnotationFunction(
    fun = function(index, k = 1, n = 1) {
      group_label <- as.character(category_split[index[[1]]])
      grid.rect(
        x = unit(0.5, "npc"), y = unit(0.5, "npc"),
        width = unit(1, "npc"), height = unit(1, "npc"),
        gp = gpar(fill = "#F1F3F5", col = "#CCD1D5", lwd = 0.55)
      )
      grid.text(
        group_label, x = unit(0.97, "npc"), y = unit(0.5, "npc"),
        just = "right", gp = gpar(fontfamily = FONT, fontsize = 8.8, fontface = "bold", col = "black")
      )
    },
    which = "row", width = unit(40, "mm")
  )

  label_annotation <- rowAnnotation(
    `Functional group` = functional_group_annotation,
    `Gene and KO` = gene_ko_annotation,
    show_annotation_name = FALSE,
    gap = unit(1.2, "mm")
  )

  Heatmap(
    matrix_value,
    name = paste0("z_", module_name),
    col = z_col_fun,
    top_annotation = group_annotation,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = category_split,
    column_split = column_split,
    cluster_column_slices = FALSE,
    row_gap = unit(0.45, "mm"),
    column_gap = unit(0.75, "mm"),
    row_title = NULL,
    left_annotation = label_annotation,
    show_row_names = FALSE,
    column_names_side = "bottom",
    column_names_rot = 0,
    column_names_centered = TRUE,
    column_names_gp = gpar(fontfamily = FONT, fontsize = 9.0, col = "black"),
    column_title = GROUP_LEVELS,
    column_title_gp = gpar(fontfamily = FONT, fontsize = 10.3, fontface = "bold"),
    rect_gp = gpar(col = "grey70", lwd = 0.35),
    border = TRUE,
    width = unit(HEATMAP_BODY_WIDTH_CM, "cm"),
    height = unit(max(6.0, nrow(matrix_value) * 0.42), "cm"),
    show_heatmap_legend = FALSE
  )
}

page_width_cm <- 18.3 + (HEATMAP_BODY_WIDTH_CM - 10.8)
panel_viewport_left_cm <- 0.05
panel_viewport_width_cm <- page_width_cm - 0.15
title_height_cm <- 0.60
panel_gap_cm <- 0.35
legend_band_height_cm <- 1.60

panel_metrics <- lapply(RENDER_MODULES, function(module_name) {
  n_rows <- marker_order %>% filter(module == module_name) %>% nrow()
  heatmap_body_height <- max(6.0, n_rows * 0.42)
  list(
    module = module_name,
    n_rows = n_rows,
    body_height_cm = heatmap_body_height,
    # ComplexHeatmap adds column-slice titles, the site strip and sample labels
    # outside the declared heatmap body. Reserve their physical height so the
    # panel viewport can clip safely without cutting those labels.
    panel_height_cm = heatmap_body_height + 2.70
  )
})
names(panel_metrics) <- RENDER_MODULES
page_height_cm <- sum(vapply(panel_metrics, `[[`, numeric(1), "panel_height_cm")) +
  panel_gap_cm * (length(RENDER_MODULES) - 1) + 0.40 + legend_band_height_cm

output <- file.path(FIGURE_DIR, OUTPUT_NAME)
quartz(
  type = "pdf",
  file = output,
  width = page_width_cm / 2.54,
  height = page_height_cm / 2.54,
  family = FONT
)
grid.newpage()
grid.rect(gp = gpar(fill = "white", col = NA))

cursor_top_cm <- page_height_cm - 0.18
geometry_rows <- list()
for (module_name in RENDER_MODULES) {
  metric <- panel_metrics[[module_name]]
  panel_top_cm <- cursor_top_cm
  panel_bottom_cm <- panel_top_cm - metric$panel_height_cm

  if (DRAW_PANEL_LETTERS) grid.text(
    PANEL_LETTERS[[module_name]],
    x = unit(0.20, "cm"), y = unit(panel_top_cm - 0.08, "cm"),
    just = c("left", "top"),
    gp = gpar(fontfamily = FONT, fontsize = 15.5, fontface = "bold")
  )

  heatmap_viewport_height_cm <- metric$panel_height_cm - title_height_cm
  pushViewport(viewport(
    x = unit(panel_viewport_left_cm, "cm"),
    y = unit(panel_bottom_cm, "cm"),
    width = unit(panel_viewport_width_cm, "cm"),
    height = unit(heatmap_viewport_height_cm, "cm"),
    just = c("left", "bottom"),
    # Prevent any ComplexHeatmap child grob from painting over the neighbouring
    # panel title or the shared legend band.
    clip = "on",
    name = paste0("panel_", module_name)
  ))
  draw(
    build_heatmap(module_name),
    newpage = FALSE,
    padding = unit(c(1.4, 1.0, 1.6, 0.8), "mm")
  )
  popViewport()

  geometry_rows[[module_name]] <- tibble(
    panel = PANEL_LETTERS[[module_name]],
    module = module_name,
    displayed_samples = 17L,
    displayed_markers = metric$n_rows,
    displayed_cells = 17L * metric$n_rows,
    panel_top_cm = panel_top_cm,
    panel_bottom_cm = panel_bottom_cm,
    heatmap_body_width_cm = HEATMAP_BODY_WIDTH_CM,
    heatmap_body_alignment = "fixed left annotation and row-name widths",
    panel_A_optical_translation_cm = 0,
    statistical_unit = "individual biological sample"
  )
  cursor_top_cm <- panel_bottom_cm - panel_gap_cm
}

shared_legend <- Legend(
  title = "Row z-score [log10(TPM + 1)]",
  title_position = "topcenter",
  col_fun = z_col_fun,
  at = c(-2, 0, 2),
  direction = "horizontal",
  legend_width = unit(31, "mm"),
  legend_height = unit(4.5, "mm"),
  title_gp = gpar(fontfamily = FONT, fontsize = 10.0, fontface = "bold"),
  labels_gp = gpar(fontfamily = FONT, fontsize = 9.2)
)
draw(
  shared_legend,
  x = unit(page_width_cm / 2, "cm"),
  y = unit(0.12, "cm"),
  just = c("center", "bottom")
)
dev.off()

geometry <- bind_rows(geometry_rows)
write_tsv(geometry, file.path(SOURCE_DIR, paste0(sub("\\.pdf$", "", OUTPUT_NAME), "_display_geometry.tsv")))
pt_per_cm <- 72 / 2.54
panel_records <- lapply(seq_len(nrow(geometry)), function(i) {
  left <- 0.12 * pt_per_cm
  right <- (page_width_cm - 0.32) * pt_per_cm
  bottom <- geometry$panel_bottom_cm[[i]] * pt_per_cm
  top <- geometry$panel_top_cm[[i]] * pt_per_cm
  list(
    id = geometry$panel[[i]], role = "panel",
    source_width_pt = right - left,
    source_height_pt = top - bottom,
    placed_box_pt = c(left, bottom, right, top),
    ink_box_pt = c(left, bottom + 0.08 * pt_per_cm, right, top - 0.04 * pt_per_cm)
  )
})
write_json(
  list(
    figure_id = paste0("geoderma_v03_split_", SPLIT_MODE, "_final_size"),
    aspect_tolerance_pct = 0.05,
    source_panels_locked_before_composition = TRUE,
    uses_outer_boxes_for_alignment = FALSE,
    page_width_cm = page_width_cm,
    page_height_cm = page_height_cm,
    shared_heatmap_body_width_cm = HEATMAP_BODY_WIDTH_CM,
    shared_heatmap_body_alignment = "fixed functional-group and gene-KO annotation widths",
    legend_is_independent_of_panel_layout = TRUE,
    panels = panel_records,
    gaps = if (length(RENDER_MODULES) == 2L) list(
      list(panel_a = "A", panel_b = "B", axis = "y", kind = "vertical", min_pt = 5, max_pt = 30)
    ) else list()
  ),
  file.path(BUILD_ROOT, paste0("geometry_ledger_", sub("\\.pdf$", "", OUTPUT_NAME), ".json")),
  pretty = TRUE, auto_unbox = TRUE, digits = 10
)

stopifnot(file.exists(output), file.info(output)$size > 0)
if (length(list.files(BUILD_ROOT, pattern = "\\.png$", recursive = TRUE))) {
  stop("Unexpected persistent PNG output")
}
message("Geoderma V03 split (", SPLIT_MODE, ") saved -> ", output)
