#!/usr/bin/env Rscript
# Geoderma V03 derivative of 463_render_v29_figure_s2.R : abundance-legend titles 10.7 -> 11.5 pt so the log10 subscript prints >= 8 pt (both pages).
# Geoderma V03 derivative of 463_render_v29_figure_s2.R : phylum/class legend labels roman (S2A phylum legend, S2B phylum + class legends); page-2 width 21.0 -> 19.0 cm (Elsevier <= 190 mm).
# V26 derivative of 450_render_v26_figure_s2.R (443_build_v26_figures_no_panel_titles.py): legend 'not detected'; page-2 continuation tag.
# V26 derivative of 437_render_v24_figure_s2.R (443_build_v26_figures_no_panel_titles.py): panel titles removed.
# V24 derivative of _tmp_392_fontpatched.R (built by 438_build_v24_figure_refinement.py). gheatmap column names in Arial; italic taxon legends; shorter K10944 header.
# Readability-first Figure S2 refinement.
#
# Panels A and B are retained as one numbered figure but exported on separate
# full-width pages. This avoids making the 221-tip bacterial tree and the
# archaeal metadata tracks unreadable when a 26-cm landscape source is reduced
# to the manuscript's 16.38-cm body width.

suppressPackageStartupMessages({
  library(ape)
  library(circlize)
  library(ComplexHeatmap)
  library(cowplot)
  library(dplyr)
  library(ggnewscale)
  library(ggplot2)
  library(ggplotify)
  library(ggtext)
  library(ggtree)
  library(grid)
  library(jsonlite)
  library(png)
  library(readr)
  library(scales)
})

set.seed(42)
PROJECT <- Sys.getenv("CHILE_PROJECT_ROOT", unset=getwd())
BUILD_ROOT <- Sys.getenv(
  "CHILE_READABILITY_FIGURE_S2_ROOT",
  file.path(PROJECT, "analysis/readability_first_figure_refinement_2026-08-30")
)
FIGURE_DIR <- file.path(BUILD_ROOT, "figures")
SOURCE_DIR <- file.path(BUILD_ROOT, "source_data")
PANEL_DIR <- file.path(BUILD_ROOT, "native_panels")
STAGE_DIR <- file.path(BUILD_ROOT, "_figure_s4_source_stage")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAGE_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(PROJECT, "script/_theme_sbb.R"))
FONT <- "Arial"
pdfFonts(Arial = pdfFonts("ArialMT")[[1]])
if (capabilities("aqua")) quartzFonts(
  Arial = quartzFont(c("Arial", "Arial Bold", "Arial Italic", "Arial Bold Italic"))
)

# Rebuild the accepted tree objects from their canonical data pipeline.
tree_env <- new.env(parent = globalenv())
old_support <- Sys.getenv("CHILE_FIGURE_S3_SUPPORT_DIR", unset = NA_character_)
old_output <- Sys.getenv("CHILE_FIGURE_S3_OUT", unset = NA_character_)
Sys.setenv(
  CHILE_FIGURE_S3_SUPPORT_DIR = file.path(STAGE_DIR, "source_data"),
  CHILE_FIGURE_S3_OUT = file.path(STAGE_DIR, "source.pdf")
)
sys.source(file.path(PROJECT, "script/212_build_v18_s3_clean.R"), envir = tree_env)
if (is.na(old_support)) Sys.unsetenv("CHILE_FIGURE_S3_SUPPORT_DIR") else Sys.setenv(CHILE_FIGURE_S3_SUPPORT_DIR = old_support)
if (is.na(old_output)) Sys.unsetenv("CHILE_FIGURE_S3_OUT") else Sys.setenv(CHILE_FIGURE_S3_OUT = old_output)

focal_mags <- c("MAG122", "MAG138", "MAG156")
tree_meta <- tree_env$mag_metrics %>%
  filter(Domain == "Bacteria") %>%
  transmute(label = MAG, phylum_plot)
phylum_counts <- tree_env$mag_metrics %>%
  filter(Domain == "Bacteria") %>%
  mutate(phylum_plot = if_else(phylum %in% names(phylum_palette), phylum, "Other")) %>%
  count(phylum_plot, name = "n_MAG")
phylum_breaks <- c(
  "Actinomycetota", "Dormibacterota", "Acidobacteriota", "Pseudomonadota",
  "Gemmatimonadota", "Chloroflexota", "Myxococcota", "Other"
)
phylum_labels <- setNames(
  paste0(phylum_counts$phylum_plot, " (", phylum_counts$n_MAG, ")"),
  phylum_counts$phylum_plot
)
# ---- Panel A: circular bacterial tree with one abundance ring ---------------
p_a_base <- ggtree(
  tree_env$bac_tree, layout = "fan", open.angle = 16,
  size = 0.20, colour = "grey52"
) %<+% tree_meta +
  geom_tippoint(
    aes(fill = phylum_plot), shape = 21, size = 1.32,
    stroke = 0.12, colour = "white", alpha = 0.96
  ) +
  geom_tippoint(
    data = function(x) filter(x, isTip, label %in% focal_mags),
    shape = 23, size = 2.05, stroke = 0.45,
    fill = "#F2C14E", colour = "black"
  ) +
  scale_fill_manual(
    values = phylum_palette,
    guide = "none"
  )

p_a <- p_a_base + ggnewscale::new_scale_fill()
p_a <- gheatmap(
  p_a, tree_env$abundance_matrix,
  offset = 0.040, width = 0.130,
  color = SBB_HEATMAP_BORDER, colnames = FALSE
) +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#9ECAE1", "#4A6D8C"),
    limits = tree_env$abundance_limits,
    breaks = pretty(tree_env$abundance_limits, n = 3),
    guide = "none"
  )

tree_radius <- max(p_a_base$data$x, na.rm = TRUE)
focal_positions <- p_a_base$data %>%
  filter(isTip, label %in% focal_mags) %>%
  mutate(
    label_x = case_when(
      label == "MAG122" ~ tree_radius + 0.025,
      label == "MAG156" ~ tree_radius + 0.025,
      TRUE ~ tree_radius + 0.025
    ),
    tangent_angle = (angle + 90) %% 360,
    tangent_angle = if_else(
      tangent_angle > 90 & tangent_angle < 270,
      (tangent_angle + 180) %% 360,
      tangent_angle
    ),
    label_angle = if_else(label == "MAG156", 0, tangent_angle)
  )

p_a <- p_a +
  geom_text(
    data = focal_positions,
    aes(x = label_x, y = y, label = label, angle = label_angle),
    inherit.aes = FALSE,
    family = FONT, fontface = "bold", size = 3.60, colour = "black",
    hjust = 0.5, vjust = 0.5
  ) +
  labs(title = NULL, subtitle = NULL) +
  theme_void(base_family = FONT) +
  theme(
    text = element_text(family = FONT, colour = "black"),
    plot.title = element_text(size = 11.7, face = "bold", hjust = 0, margin = margin(l = 28, b = 3)),
    legend.position = "none",
    plot.margin = margin(4, 1, 1, 3),
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA)
  )

# ---- Panel B: archaeal tree and aligned metadata tracks ----------------------
tree_theme <- theme(
  text = element_text(family = FONT, color = "black"),
  plot.title = element_text(
    family = FONT, size = 11.7, face = "bold", hjust = 0,
    margin = margin(l = 28)
  ),
    plot.margin = margin(3, 3, 2, 2)
)

tree_env$ar_meta <- tree_env$ar_meta %>%
  mutate(
    tip_display = if_else(
      MAG %in% c("MAG005", "MAG070"),
      paste0("<b>", tip_display, "</b>"),
      tip_display
    )
  )
ar_tree_visual <- tree_env$ar_tree
ar_tree_visual$edge.length <- ar_tree_visual$edge.length * 1.20
p_b <- ggtree(
  ar_tree_visual, layout = "rectangular", size = 0.40, color = "grey28"
) %<+% tree_env$ar_meta +
  ggtext::geom_richtext(
    data = function(x) filter(x, !is.na(tip_display)),
    aes(label = tip_display), hjust = 0, vjust = 0.5, nudge_x = 0.27,
    size = 3.55, color = "grey10", family = FONT,
    fill = NA, label.color = NA, label.padding = unit(0, "pt"),
    label.r = unit(0, "pt")
  ) +
  geom_tippoint(aes(color = phylum_plot), size = 2.65, alpha = 0.92) +
  scale_color_manual(values = phylum_palette, guide = "none") +
  labs(title = NULL, subtitle = NULL) +
  xlim(0, 21.50) +
  ylim(-0.60, 9.60) +
  tree_theme +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA)
  )

p_b_class <- tree_env$ar_class_matrix
colnames(p_b_class) <- "Class"
p_b <- gheatmap(
  p_b, p_b_class,
  offset = 2.55, width = 0.55, color = SBB_HEATMAP_BORDER,
  colnames = TRUE, colnames_position = "top", colnames_angle = 0,
  hjust = 0.5, colnames_offset_y = 0.08, font.size = 3.25, family = FONT
) +
  scale_fill_manual(values = tree_env$archaeal_class_palette, guide = "none")

p_b <- p_b + ggnewscale::new_scale_fill()
p_b_completeness <- tree_env$ar_completeness_matrix
  colnames(p_b_completeness) <- "Completeness\n(%)"
p_b <- gheatmap(
  p_b, p_b_completeness,
  offset = 4.15, width = 0.85, color = SBB_HEATMAP_BORDER,
  colnames = TRUE, colnames_position = "top", colnames_angle = 0,
  hjust = 0.5, colnames_offset_y = 0.22, font.size = 3.30, family = FONT
) +
  scale_fill_gradientn(
    colours = tree_env$completeness_colours,
    values = scales::rescale(c(50, 65, 80, 90, 100), from = c(50, 100)),
    limits = c(50, 100), oob = squish, guide = "none"
  )

p_b <- p_b + ggnewscale::new_scale_fill()
ar_abundance_display <- tree_env$ar_abundance_matrix
stopifnot(ncol(ar_abundance_display) == 6)
colnames(ar_abundance_display) <- c("AZ", "SG", "LC", "NA-S", "NA-M", "NA-D")
p_b <- gheatmap(
  p_b, ar_abundance_display,
  offset = 5.85, width = 9.40, color = SBB_HEATMAP_BORDER,
  colnames = TRUE, colnames_position = "top", colnames_angle = 0,
  hjust = 0.5, colnames_offset_y = 0.08, font.size = 3.30, family = FONT
) +
  scale_fill_gradientn(
    colours = pal_seq,
    limits = c(-3, 1.5), oob = squish, guide = "none"
  )

p_b <- p_b + ggnewscale::new_scale_fill()
p_b_k10944 <- tree_env$ar_pmo_matrix
colnames(p_b_k10944) <- "K10944\nAmoA/PmoA"
p_b <- gheatmap(
  p_b, p_b_k10944,
  offset = 16.55, width = 1.30, color = SBB_HEATMAP_BORDER,
  colnames = TRUE, colnames_position = "top", colnames_angle = 0,
  hjust = 0.5, colnames_offset_y = 0.08, font.size = 3.30, family = FONT
) +
  scale_fill_manual(values = tree_env$pmo_palette, guide = "none") +
  theme(legend.position = "none")

# ---- Dedicated legends ------------------------------------------------------
a_phylum_labels <- unname(phylum_labels[phylum_breaks])
a_phylum_labels[is.na(a_phylum_labels)] <- phylum_breaks[is.na(a_phylum_labels)]
a_phylum_legend <- Legend(
  title = "Phylum",
  # Geoderma: phylum names roman (no italics); "Other" unchanged.
  labels = a_phylum_labels,
  type = "points", pch = 16,
  legend_gp = gpar(col = unname(phylum_palette[phylum_breaks])),
  background = "white", border = NA, size = unit(4.0, "mm"),
  ncol = 2,
  title_gp = gpar(fontfamily = FONT, fontsize = 10.7, fontface = "bold"),
  labels_gp = gpar(fontfamily = FONT, fontsize = 9.6),
  grid_width = unit(5.2, "mm"), grid_height = unit(5.2, "mm"),
  gap = unit(1.2, "mm")
)
a_abundance_fun <- colorRamp2(
  seq(tree_env$abundance_limits[[1]], tree_env$abundance_limits[[2]], length.out = length(pal_seq)),
  pal_seq
)
a_abundance_legend <- Legend(
  title = expression(atop("Mean relative abundance", log[10] * " per-mille")),
  col_fun = a_abundance_fun,
  at = pretty(tree_env$abundance_limits, n = 3),
  direction = "horizontal",
  legend_width = unit(34, "mm"),
  # Geoderma: 11.5 pt so the plotmath subscript "10" (0.7x) prints >= 8 pt.
  title_gp = gpar(fontfamily = FONT, fontsize = 11.5, fontface = "bold"),
  labels_gp = gpar(fontfamily = FONT, fontsize = 9.6)
)
a_focal_legend <- Legend(
  title = "Selected local bacterial MAGs",
  labels = "MAG122, MAG138, MAG156",
  type = "points", pch = 23,
  legend_gp = gpar(fill = "#F2C14E", col = "black"),
  background = "white", border = NA, size = unit(4.0, "mm"),
  title_gp = gpar(fontfamily = FONT, fontsize = 10.7, fontface = "bold"),
  labels_gp = gpar(fontfamily = FONT, fontsize = 9.6),
  grid_width = unit(5.2, "mm"), grid_height = unit(5.2, "mm")
)
a_legend_grob <- grid.grabExpr(draw(packLegend(
  a_abundance_legend, a_focal_legend, a_phylum_legend,
  direction = "horizontal", gap = unit(4, "mm"), max_width = unit(17.4, "cm")
)))

b_phylum_levels <- tree_env$archaeal_phylum_levels
b_class_levels <- tree_env$archaeal_class_levels
b_phylum_legend <- Legend(
  title = "Phylum", labels = b_phylum_levels,
  type = "points", pch = 16,
  legend_gp = gpar(col = unname(phylum_palette[b_phylum_levels])),
  background = "white", border = NA, size = unit(3.0, "mm"),
  title_gp = gpar(fontfamily = FONT, fontsize = 10.7, fontface = "bold"),
  labels_gp = gpar(fontfamily = FONT, fontsize = 9.6, fontface = "plain"),
  grid_width = unit(3.3, "mm"), grid_height = unit(3.3, "mm")
)
b_class_legend <- Legend(
  title = "Class", labels = b_class_levels,
  legend_gp = gpar(fill = unname(tree_env$archaeal_class_palette[b_class_levels]), col = NA),
  title_gp = gpar(fontfamily = FONT, fontsize = 10.7, fontface = "bold"),
  labels_gp = gpar(fontfamily = FONT, fontsize = 9.6, fontface = "plain"),
  grid_width = unit(3.3, "mm"), grid_height = unit(3.3, "mm")
)
b_pmo_legend <- Legend(
  title = "K10944 (AmoA/PmoA family)",
  labels = c("not detected", "present"),
  legend_gp = gpar(fill = unname(tree_env$pmo_palette[c("absent", "present")]), col = NA),
  title_gp = gpar(fontfamily = FONT, fontsize = 10.5, fontface = "bold"),
  labels_gp = gpar(fontfamily = FONT, fontsize = 9.6),
  grid_width = unit(3.3, "mm"), grid_height = unit(3.3, "mm")
)
b_categorical_grob <- grid.grabExpr(draw(packLegend(
  b_phylum_legend, b_class_legend, b_pmo_legend,
  direction = "vertical", gap = unit(1.5, "mm")
)))

b_abundance_fun <- colorRamp2(c(-3, -1.5, 0, 1.5), colorRampPalette(pal_seq)(4))
b_completeness_fun <- colorRamp2(
  c(50, 65, 80, 90, 100), tree_env$completeness_colours
)
b_abundance_legend <- Legend(
  title = expression(atop("Relative abundance", log[10] * " per-mille")),
  col_fun = b_abundance_fun, at = c(-3, -1.5, 0, 1.5),
  direction = "horizontal", legend_width = unit(34, "mm"),
  # Geoderma: 11.5 pt so the plotmath subscript "10" (0.7x) prints >= 8 pt.
  title_gp = gpar(fontfamily = FONT, fontsize = 11.5, fontface = "bold"),
  labels_gp = gpar(fontfamily = FONT, fontsize = 9.6)
)
b_completeness_legend <- Legend(
  title = "Completeness (%)",
  col_fun = b_completeness_fun, at = c(50, 75, 100),
  direction = "horizontal", legend_width = unit(34, "mm"),
  title_gp = gpar(fontfamily = FONT, fontsize = 10.7, fontface = "bold"),
  labels_gp = gpar(fontfamily = FONT, fontsize = 9.6)
)
b_continuous_grob <- grid.grabExpr(draw(packLegend(
  b_abundance_legend, b_completeness_legend,
  direction = "vertical", gap = unit(2.5, "mm")
)))

a_legend_plot <- ggdraw() +
  theme(plot.background = element_rect(fill = NA, colour = NA)) +
  draw_plot(as.ggplot(a_legend_grob), x = 0.025, y = 0.03, width = 0.95, height = 0.94)
b_legend_plot <- plot_grid(
  ggdraw() + theme(plot.background = element_rect(fill = NA, colour = NA)) +
    draw_plot(as.ggplot(b_categorical_grob), x = 0.035, y = 0.00, width = 0.93, height = 0.90),
  ggdraw() + theme(plot.background = element_rect(fill = NA, colour = NA)) +
    draw_plot(as.ggplot(b_continuous_grob), x = 0.035, y = 0.03, width = 0.93, height = 0.94),
  nrow = 1, rel_widths = c(1.04, 0.96), align = "h"
) + theme(plot.background = element_rect(fill = NA, colour = NA))

# Lock each source panel at full manuscript width. Figure S2 remains one
# numbered figure, but uses two pages so neither phylogeny is downscaled to an
# unreadable landscape thumbnail.
page_width_cm <- 18.3
page_b_width_cm <- 19.0
page_height_cm <- 19.8
row_gap_cm <- 0.20
legend_height_cm <- 5.40
top_height_cm <- page_height_cm - row_gap_cm - legend_height_cm
# User-requested 70% reduction relative to the current Figure S2B data body.
top_height_b_cm <- top_height_cm * 0.49
page_b_height_cm <- top_height_b_cm + row_gap_cm + legend_height_cm
top_bottom_cm <- legend_height_cm + row_gap_cm

panel_a_pdf <- file.path(PANEL_DIR, "Figure_S2_A.pdf")
panel_b_pdf <- file.path(PANEL_DIR, "Figure_S2_B.pdf")
legend_a_pdf <- file.path(PANEL_DIR, "Figure_S2_A_legend.pdf")
legend_b_pdf <- file.path(PANEL_DIR, "Figure_S2_B_legend.pdf")
ggsave(panel_a_pdf, p_a, width = page_width_cm, height = top_height_cm,
       units = "cm", device = cairo_pdf, bg = "transparent")
ggsave(panel_b_pdf, p_b, width = page_b_width_cm, height = top_height_b_cm,
       units = "cm", device = cairo_pdf, bg = "transparent")
ggsave(legend_a_pdf, a_legend_plot, width = page_width_cm, height = legend_height_cm,
       units = "cm", device = cairo_pdf, bg = "transparent")
ggsave(legend_b_pdf, b_legend_plot, width = page_b_width_cm, height = legend_height_cm,
       units = "cm", device = cairo_pdf, bg = "transparent")

page_a <- ggdraw() +
  theme(plot.background = element_rect(fill = "white", colour = NA)) +
  draw_grob(
    grid::rectGrob(gp = grid::gpar(fill = "white", col = NA)),
    x = 0, y = 0, width = 1, height = 1
  ) +
  draw_plot(
    p_a,
    x = 0, y = top_bottom_cm / page_height_cm,
    width = 1,
    height = top_height_cm / page_height_cm
  ) +
  draw_plot(
    a_legend_plot,
    x = 0, y = 0,
    width = 1,
    height = legend_height_cm / page_height_cm
  ) +
  draw_label(
    "A", x = 0.010, y = 0.982, hjust = 0, vjust = 1,
    fontfamily = FONT, fontface = "bold", size = 15.5
  )

page_b <- ggdraw() +
  theme(plot.background = element_rect(fill = "white", colour = NA)) +
  draw_grob(
    grid::rectGrob(gp = grid::gpar(fill = "white", col = NA)),
    x = 0, y = 0, width = 1, height = 1
  ) +
  draw_plot(
    p_b,
    x = 0, y = top_bottom_cm / page_b_height_cm,
    width = 1,
    height = top_height_b_cm / page_b_height_cm
  ) +
  draw_plot(
    b_legend_plot,
    x = 0, y = 0,
    width = 1,
    height = legend_height_cm / page_b_height_cm
  ) +
  draw_label(
    "B", x = 0.010, y = 0.982,
    hjust = 0, vjust = 1, fontfamily = FONT, fontface = "bold", size = 15.5
  ) +
  draw_label(
    "Fig. S2 (continued)", x = 0.990, y = 0.982,
    hjust = 1, vjust = 1, fontfamily = FONT, size = 9.5, colour = "grey30"
  )

output <- file.path(FIGURE_DIR, "Figure_S2.pdf")
page_a_pdf <- file.path(PANEL_DIR, "Figure_S2_page_A.pdf")
page_b_pdf <- file.path(PANEL_DIR, "Figure_S2_page_B.pdf")
ggsave(
  page_a_pdf, page_a,
  width = page_width_cm, height = page_height_cm, units = "cm",
  device = cairo_pdf, bg = "white"
)
ggsave(
  page_b_pdf, page_b,
  width = page_b_width_cm, height = page_b_height_cm, units = "cm",
  device = cairo_pdf, bg = "white"
)
merge_status <- system2("pdfunite", c(shQuote(page_a_pdf), shQuote(page_b_pdf), shQuote(output)))
if (!identical(merge_status, 0L)) stop("Could not combine the two Figure S2 pages.")
write_tsv(
  tree_env$mag_metrics %>%
    filter(Domain == "Bacteria") %>%
    transmute(
      MAG, phylum, focal = MAG %in% focal_mags,
      mean_relative_abundance = mean_rel_abundance,
      mean_per_mille, log10_mean_per_mille
    ),
  file.path(SOURCE_DIR, "Figure_S2A_bacterial_tree_abundance.tsv")
)
write_tsv(
  tree_env$ar_meta %>%
    select(MAG, phylum, class, Completeness, pmo_present),
  file.path(SOURCE_DIR, "Figure_S2B_archaeal_metadata.tsv")
)
# ---- Geometry ledger -------------------------------------------------------
pt_per_cm <- 72 / 2.54
tmp_dir <- tempfile("v45_figure_s4_geometry_")
dir.create(tmp_dir)
on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

measure_ink <- function(pdf, width_pt, height_pt) {
  prefix <- file.path(tmp_dir, tools::file_path_sans_ext(basename(pdf)))
  status <- system2(
    "pdftoppm",
    c("-png", "-singlefile", "-r", "300", shQuote(pdf), shQuote(prefix)),
    stdout = FALSE, stderr = FALSE
  )
  if (!identical(status, 0L)) stop("Could not rasterize geometry source: ", pdf)
  image <- readPNG(paste0(prefix, ".png"))
  rgb <- image[, , seq_len(min(3, dim(image)[3])), drop = FALSE]
  ink <- rgb[, , 1] < 0.985 | rgb[, , 2] < 0.985 | rgb[, , 3] < 0.985
  rows <- which(rowSums(ink) > 0)
  cols <- which(colSums(ink) > 0)
  c(
    left = (min(cols) - 1) / ncol(ink) * width_pt,
    bottom = (nrow(ink) - max(rows)) / nrow(ink) * height_pt,
    right = max(cols) / ncol(ink) * width_pt,
    top = (nrow(ink) - min(rows) + 1) / nrow(ink) * height_pt
  )
}

panel_record <- function(id, pdf, width_cm, height_cm, placed_box_cm, role = "panel", page = 1L) {
  source_width_pt <- width_cm * pt_per_cm
  source_height_pt <- height_cm * pt_per_cm
  placed_box_pt <- placed_box_cm * pt_per_cm
  source_ink <- measure_ink(pdf, source_width_pt, source_height_pt)
  list(
    id = id,
    role = role,
    page = page,
    source_width_pt = source_width_pt,
    source_height_pt = source_height_pt,
    placed_box_pt = unname(placed_box_pt),
    ink_box_pt = unname(c(
      placed_box_pt[[1]] + source_ink[[1]],
      placed_box_pt[[2]] + source_ink[[2]],
      placed_box_pt[[1]] + source_ink[[3]],
      placed_box_pt[[2]] + source_ink[[4]]
    ))
  )
}

ledger <- list(
  figure_id = "readability_first_Figure_S2_MAG_phylogenetic_context",
  aspect_tolerance_pct = 0.05,
  source_panels_locked_before_composition = TRUE,
  uses_outer_boxes_for_alignment = FALSE,
  page_count = 2L,
  page_size_cm = list(
    page_1 = c(width = page_width_cm, height = page_height_cm),
    page_2 = c(width = page_b_width_cm, height = page_b_height_cm)
  ),
  alignment = list(
    data_bodies = "A and B each use a full-width native box on a separate page",
    legends = "detached lower row on each page; excluded from data-body sizing",
    backgrounds = "all child canvases transparent; one page-wide white background per page; panel tags drawn last"
  ),
  panels = list(
    panel_record(
      "A", panel_a_pdf, page_width_cm, top_height_cm,
      c(0, top_bottom_cm, page_width_cm, page_height_cm), page = 1L
    ),
    panel_record(
      "B", panel_b_pdf, page_b_width_cm, top_height_b_cm,
      c(0, top_bottom_cm, page_b_width_cm, page_b_height_cm), page = 2L
    ),
    panel_record(
      "A_legend", legend_a_pdf, page_width_cm, legend_height_cm,
      c(0, 0, page_width_cm, legend_height_cm), role = "legend", page = 1L
    ),
    panel_record(
      "B_legend", legend_b_pdf, page_b_width_cm, legend_height_cm,
      c(0, 0, page_b_width_cm, legend_height_cm), role = "legend", page = 2L
    )
  ),
  gaps = list(
    list(
      panel_a = "A", panel_b = "A_legend", axis = "y",
      kind = "vertical", min_pt = 40, max_pt = 75
    ),
    list(
      panel_a = "B", panel_b = "B_legend", axis = "y",
      kind = "vertical", min_pt = 35, max_pt = 95
    )
  )
)

write_json(
  ledger,
  file.path(BUILD_ROOT, "geometry_ledger_Figure_S2.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = 10
)

unlink(STAGE_DIR, recursive = TRUE, force = TRUE)
stopifnot(file.exists(output), file.info(output)$size > 0)
if (length(list.files(BUILD_ROOT, pattern = "\\.png$", recursive = TRUE))) {
  stop("Unexpected persistent PNG output")
}
message("Readability-first two-page Figure S2 saved -> ", output)
