#!/usr/bin/env Rscript
# V26 derivative of 434_render_v24_figure2_s1.R (443_build_v26_figures_no_panel_titles.py): panel titles removed (Fig. 2A–C, Fig. S1A–D); explained variance moves to the Fig. S1 caption.
# V24 derivative of 422_render_v19_figure_s1_pit_labels.R (built by 438_build_v24_figure_refinement.py). Fig. 2C zero row without stray dot; Fig. S1 unchanged from V20 (author decision 2026-09-03: no pit trajectories).
# Readability-first render of Figure 2 and Figure S1.
#
# A/B/C in Figure 2 answer three different questions:
#   A. Is surface KO composition differentiated among sites?
#   B. How is that composition arranged along the two frozen soil-property axes?
#   C. How broad is FDR-supported KO differentiation in the two study designs?
#
# All ordination coordinates and inferential labels are read from frozen outputs.
# Soil-PC arrows are descriptive orientations derived from frozen site scores and
# frozen PCA scores with permutations = 0; no model or hypothesis test is rerun.

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggrepel)
  library(ggplot2)
  library(grid)
  library(jsonlite)
  library(png)
  library(readr)
  library(tibble)
  library(vegan)
})

PROJECT <- Sys.getenv("CHILE_PROJECT_ROOT", unset=getwd())
OUT <- Sys.getenv(
  "CHILE_READABILITY_FIGURE2_ROOT",
  file.path(PROJECT, "analysis/V19_peer_review_targeted_revision_2026-09-01/figure_revision")
)
FIG_DIR <- file.path(OUT, "figures")
SRC_DIR <- file.path(OUT, "source_data")
PANEL_DIR <- file.path(OUT, "native_panels")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
source(file.path(PROJECT, "script/_theme_sbb.R"))

FONT <- "Arial"
pdfFonts(Arial = pdfFonts("ArialMT")[[1]])

INPUTS <- c(
  pcoa = file.path(PROJECT, "analysis/V53_author_comments_revision_2026-08-28/source_data/direct_location_pcoa_coordinates.tsv"),
  surface_sites = file.path(PROJECT, "analysis/V59_KO_ordination_architecture_2026-08-29/source_data/Figure_2B_surface_KO_dbrda_coordinates.tsv"),
  depth_sites = file.path(PROJECT, "analysis/V59_KO_ordination_architecture_2026-08-29/source_data/Figure_S3B_depth_KO_dbrda_coordinates.tsv"),
  dbrda_stats = file.path(PROJECT, "analysis/V59_KO_ordination_architecture_2026-08-29/source_data/KO_dbrda_frozen_model_statistics.tsv"),
  soil_scores = file.path(PROJECT, "analysis/V05_comment_review_2026-08-17/figure2_reanalysis/soil_pca_scores.tsv"),
  soil_loadings = file.path(PROJECT, "analysis/V05_comment_review_2026-08-17/figure2_reanalysis/soil_pca_loadings.tsv"),
  soil_variance = file.path(PROJECT, "analysis/V05_comment_review_2026-08-17/figure2_reanalysis/soil_pca_variance.tsv"),
  da_readme = file.path(PROJECT, "analysis/exploratory_KO_differential_abundance_2026-08-29/KO_DA_README.md")
)
if (any(!file.exists(INPUTS))) {
  stop("Missing frozen input(s): ", paste(INPUTS[!file.exists(INPUTS)], collapse = "; "))
}
invisible(vapply(INPUTS, function(path) {
  file.copy(path, file.path(SRC_DIR, basename(path)), overwrite = TRUE)
}, logical(1)))

GROUP_LEVELS <- c("AZ", "SG", "LC", "NA-S", "NA-M", "NA-D")
GROUP_COLOURS <- c(
  AZ = pal_site[["AZ"]], SG = pal_site[["SG"]], LC = pal_site[["LC"]],
  `NA-S` = pal_depth[["5-10 cm"]], `NA-M` = pal_depth[["20-40 cm"]],
  `NA-D` = pal_depth[["60-80 cm"]]
)
GROUP_SHAPES <- c(AZ = 21, SG = 21, LC = 21, `NA-S` = 21, `NA-M` = 22, `NA-D` = 24)
SURFACE_GROUPS <- c("AZ", "SG", "LC", "NA-S")
DEPTH_GROUPS <- c("NA-S", "NA-M", "NA-D")

PRETTY_ENV <- c(
  pH = "pH", Conductivity = "Conductivity", moisture = "Moisture",
  C = "Total C", N = "Total N", Feo = "Feₒₓ", Alo = "Alₒₓ",
  Mno = "Mnₒₓ", Sio = "Siₒₓ", NH4 = "NH₄⁺",
  NO3 = "NO₃⁻", Pt = "Total P"
)
ENV_VARS <- names(PRETTY_ENV)

stat_label <- function(r2, adjusted_p) {
  sprintf("plain(R)^2 == %.3f*';'~'adjusted'~italic(p) == %s", r2, adjusted_p)
}

square_limits <- function(x, y, padding = 0.14) {
  xr <- range(x, finite = TRUE)
  yr <- range(y, finite = TRUE)
  span <- max(diff(xr), diff(yr)) * (1 + 2 * padding)
  list(
    x = mean(xr) + c(-0.5, 0.5) * span,
    y = mean(yr) + c(-0.5, 0.5) * span
  )
}

theme_master <- function(base_size = 10.5) {
  theme_sbb(base_size = base_size, base_family = FONT) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(size = 11.7, face = "bold", margin = margin(l = 28, b = 4)),
      plot.tag = element_text(size = 15.5, face = "bold", colour = "black"),
      plot.tag.position = c(0.010, 0.990),
      axis.title = element_text(size = 11.2),
      axis.text = element_text(size = 10.2),
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.36),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.42),
      plot.margin = margin(6, 7, 6, 7),
      plot.background = element_rect(fill = NA, colour = NA),
      panel.background = element_rect(fill = NA, colour = NA)
    )
}

read_pcoa <- function(subset_name) {
  p <- read_tsv(INPUTS[["pcoa"]], show_col_types = FALSE) %>%
    filter(profile == "Fractional KO catalog", subset == subset_name)
  if (subset_name == "surface_site") {
    p <- p %>% mutate(
      display_group = factor(if_else(site == "NB", "NA-S", as.character(site)), levels = GROUP_LEVELS)
    )
  } else {
    p <- p %>% mutate(
      display_group = factor(
        recode(depth_label, `5_10cm` = "NA-S", `20_40cm` = "NA-M", `60_80cm` = "NA-D"),
        levels = GROUP_LEVELS
      )
    )
  }
  p
}

make_pcoa_panel <- function(
  data, title, r2, adjusted_p, tag, padding = 0.14, show_pit_labels = FALSE
) {
  lim <- square_limits(data$axis1, data$axis2, padding)
  ggplot(data, aes(axis1, axis2)) +
    geom_hline(yintercept = 0, colour = "grey86", linewidth = 0.22) +
    geom_vline(xintercept = 0, colour = "grey86", linewidth = 0.22) +
    geom_point(
      aes(fill = display_group, shape = display_group),
      colour = "grey12", size = 3.20, stroke = 0.46, show.legend = FALSE
    ) +
    geom_text_repel(
      data = if (show_pit_labels) data else data[0, ],
      aes(label = site_pit),
      family = FONT, fontface = "bold", size = 3.05, colour = "grey15",
      box.padding = 0.16, point.padding = 0.18, min.segment.length = 0,
      segment.colour = "grey45", segment.size = 0.24,
      seed = 42, max.overlaps = Inf, show.legend = FALSE
    ) +
    annotate(
      "text",
      x = lim$x[1] + 0.038 * diff(lim$x),
      y = lim$y[2] - 0.010 * diff(lim$y),
      label = stat_label(r2, adjusted_p), parse = TRUE,
      family = FONT, size = 3.65, hjust = 0, vjust = 1
    ) +
    scale_fill_manual(values = GROUP_COLOURS, breaks = GROUP_LEVELS, drop = FALSE) +
    scale_shape_manual(values = GROUP_SHAPES, breaks = GROUP_LEVELS, drop = FALSE) +
    coord_equal(xlim = lim$x, ylim = lim$y, expand = FALSE, clip = "off") +
    labs(
      tag = tag, title = title,
      x = sprintf("PCoA1 (%.1f%%)", 100 * unique(data$axis1_variance)),
      y = sprintf("PCoA2 (%.1f%%)", 100 * unique(data$axis2_variance))
    ) +
    theme_master()
}

soil_scores <- read_tsv(INPUTS[["soil_scores"]], show_col_types = FALSE) %>%
  select(sample_id, PC1, PC2)

make_soil_pc_vectors <- function(sites, subset_name) {
  joined <- sites %>%
    left_join(soil_scores, by = "sample_id")
  if (anyNA(joined[, c("PC1", "PC2")])) {
    stop("Missing frozen Soil PC scores for ", subset_name, " dbRDA sites.")
  }
  site_matrix <- as.matrix(joined[, c("dbRDA1", "dbRDA2")])
  rownames(site_matrix) <- joined$sample_id
  pc_matrix <- as.matrix(joined[, c("PC1", "PC2")])
  rownames(pc_matrix) <- joined$sample_id

  fit <- vegan::envfit(site_matrix, pc_matrix, permutations = 0)
  raw <- vegan::scores(fit, display = "vectors")[c("PC1", "PC2"), , drop = FALSE]
  site_radius <- max(sqrt(rowSums(site_matrix^2)))
  vector_radius <- max(sqrt(rowSums(raw^2)))
  scaled <- raw * (0.66 * site_radius / vector_radius)

  vectors <- as_tibble(scaled, rownames = "soil_axis") %>%
    rename(xend = dbRDA1, yend = dbRDA2) %>%
    mutate(
      label = recode(soil_axis, PC1 = "Soil PC1", PC2 = "Soil PC2"),
      fitted_r2 = unname(fit$vectors$r[soil_axis]),
      x_label = xend * 1.16,
      y_label = yend * 1.16,
      label_hjust = if_else(x_label < 0, 1, 0)
    )

  if (subset_name == "surface") {
    vectors <- vectors %>% mutate(
      x_label = case_when(soil_axis == "PC1" ~ -0.53, TRUE ~ -0.53),
      y_label = case_when(soil_axis == "PC1" ~ -0.33, TRUE ~ 0.43),
      label_hjust = 0.5
    )
  } else if (subset_name == "depth") {
    vectors <- vectors %>% mutate(
      x_label = case_when(soil_axis == "PC1" ~ 0.47, TRUE ~ 0.46),
      y_label = case_when(soil_axis == "PC1" ~ -0.36, TRUE ~ 0.28),
      label_hjust = 0
    )
  }
  vectors
}

make_dbrda_panel <- function(sites, vectors, title, r2, adjusted_p, tag, axis_fraction) {
  lim <- square_limits(
    c(sites$dbRDA1, vectors$x_label),
    c(sites$dbRDA2, vectors$y_label),
    0.17
  )
  ggplot(sites, aes(dbRDA1, dbRDA2)) +
    geom_hline(yintercept = 0, colour = "grey86", linewidth = 0.22) +
    geom_vline(xintercept = 0, colour = "grey86", linewidth = 0.22) +
    geom_point(
      aes(fill = display_group, shape = display_group),
      colour = "grey12", size = 3.20, stroke = 0.46, show.legend = FALSE
    ) +
    geom_segment(
      data = vectors,
      aes(x = 0, y = 0, xend = xend, yend = yend),
      inherit.aes = FALSE, colour = "black", linewidth = 0.38,
      lineend = "round", arrow = arrow(length = unit(1.55, "mm"), type = "open")
    ) +
    geom_text(
      data = vectors,
      aes(x = x_label, y = y_label, label = label, hjust = label_hjust),
      inherit.aes = FALSE, family = FONT, size = 3.45, colour = "black"
    ) +
    annotate(
      "text",
      x = lim$x[1] + 0.038 * diff(lim$x),
      y = lim$y[2] - 0.010 * diff(lim$y),
      label = stat_label(r2, adjusted_p), parse = TRUE,
      family = FONT, size = 3.65, hjust = 0, vjust = 1
    ) +
    scale_fill_manual(values = GROUP_COLOURS, breaks = GROUP_LEVELS, drop = FALSE) +
    scale_shape_manual(values = GROUP_SHAPES, breaks = GROUP_LEVELS, drop = FALSE) +
    coord_equal(xlim = lim$x, ylim = lim$y, expand = FALSE, clip = "off") +
    labs(
      tag = tag, title = title,
      x = sprintf("dbRDA1 (%.1f%%)", 100 * axis_fraction[[1]]),
      y = sprintf("dbRDA2 (%.1f%%)", 100 * axis_fraction[[2]])
    ) +
    theme_master()
}

make_site_legend <- function() {
  dat <- tibble(
    group = factor(SURFACE_GROUPS, levels = GROUP_LEVELS),
    display_label = c("AZ", "SG", "LC", "NA"),
    x = c(0.3, 2.0, 3.7, 5.4), y = 0
  )
  ggplot(dat, aes(x, y)) +
    geom_point(
      aes(fill = group, shape = group), colour = "grey12",
      size = 3.15, stroke = 0.44, show.legend = FALSE
    ) +
    geom_text(aes(label = display_label), nudge_x = 0.22, hjust = 0, family = FONT, size = 3.50) +
    annotate(
      "text", x = mean(range(dat$x)), y = 0.72,
      label = "Site", family = FONT, fontface = "bold", size = 3.80
    ) +
    scale_fill_manual(values = GROUP_COLOURS) +
    scale_shape_manual(values = GROUP_SHAPES) +
    coord_cartesian(
      xlim = c(-0.10, 6.25), ylim = c(-0.35, 0.98), clip = "off"
    ) +
    theme_void(base_family = FONT) +
    theme(
      plot.background = element_rect(fill = NA, colour = NA),
      panel.background = element_rect(fill = NA, colour = NA)
    )
}

make_depth_legend <- function() {
  dat <- tibble(
    group = factor(DEPTH_GROUPS, levels = GROUP_LEVELS),
    display_label = c("5–10 cm", "20–40 cm", "60–80 cm"),
    x = 0.28, y = c(2.65, 1.65, 0.65)
  )
  ggplot(dat, aes(x, y)) +
    geom_point(
      aes(fill = group, shape = group), colour = "grey12",
      size = 3.15, stroke = 0.44, show.legend = FALSE
    ) +
    geom_text(
      aes(label = display_label), nudge_x = 0.28, hjust = 0,
      family = FONT, size = 3.45
    ) +
    annotate(
      "text", x = 0.10, y = 3.38, label = "Depths",
      family = FONT, fontface = "bold", size = 3.80, hjust = 0
    ) +
    scale_fill_manual(values = GROUP_COLOURS) +
    scale_shape_manual(values = GROUP_SHAPES) +
    coord_cartesian(xlim = c(0, 2.45), ylim = c(0.25, 3.62), clip = "off") +
    theme_void(base_family = FONT) +
    theme(
      plot.background = element_rect(fill = NA, colour = NA),
      panel.background = element_rect(fill = NA, colour = NA)
    )
}

# ---- Frozen statistics and coordinate audits -------------------------------
dbrda_stats <- read_tsv(INPUTS[["dbrda_stats"]], show_col_types = FALSE)
surface_stats <- dbrda_stats %>% filter(subset == "surface")
depth_stats <- dbrda_stats %>% filter(subset == "depth")
if (
  nrow(surface_stats) != 1L || nrow(depth_stats) != 1L ||
  abs(surface_stats$R2[[1]] - 0.3148246034256495) > 1e-10 ||
  surface_stats$adjusted_p_label[[1]] != "0.064" ||
  abs(depth_stats$R2[[1]] - 0.3510658546852085) > 1e-10 ||
  depth_stats$adjusted_p_label[[1]] != "0.204"
) {
  stop("Frozen dbRDA statistics failed the audit.")
}

surface_pcoa_data <- read_pcoa("surface_site")
depth_pcoa_data <- read_pcoa("nahuelbuta_depth")
if (nrow(surface_pcoa_data) != 11L || nrow(depth_pcoa_data) != 9L) {
  stop("Frozen KO PCoA sample counts changed.")
}

surface_sites <- read_tsv(INPUTS[["surface_sites"]], show_col_types = FALSE) %>%
  mutate(display_group = factor(display_group, levels = GROUP_LEVELS))
depth_sites <- read_tsv(INPUTS[["depth_sites"]], show_col_types = FALSE) %>%
  mutate(display_group = factor(display_group, levels = GROUP_LEVELS))
if (nrow(surface_sites) != 11L || nrow(depth_sites) != 9L) {
  stop("Frozen KO dbRDA sample counts changed.")
}

surface_pc_vectors <- make_soil_pc_vectors(surface_sites, "surface")
depth_pc_vectors <- make_soil_pc_vectors(depth_sites, "depth")
write_tsv(surface_pc_vectors, file.path(SRC_DIR, "Figure_2B_surface_SoilPC_vectors.tsv"))
write_tsv(depth_pc_vectors, file.path(SRC_DIR, "Figure_S1D_depth_SoilPC_vectors.tsv"))

# The two count summaries are frozen in the publication-facing workbook and
# independently recorded in the local KO-DA audit README. Assert the README so
# this renderer cannot silently drift from the established values.
da_lines <- read_lines(INPUTS[["da_readme"]])
if (
  !any(grepl("Surface sites \\| 7,841 \\| 1,403 \\(17.9%\\)", da_lines)) ||
  !any(grepl("Nahuelbuta depth, blocked \\| 7,801 \\| 0 \\(0%\\)", da_lines))
) {
  stop("KO differential-abundance summary does not match the frozen audit record.")
}
da_summary <- tibble(
  comparison = factor(
    c("Surface sites", "Nahuelbuta depths"),
    levels = c("Nahuelbuta depths", "Surface sites")
  ),
  tested = c(7841L, 7801L),
  supported = c(1403L, 0L),
  percent = c(17.9, 0),
  label = c("1,403 / 7,841 (17.9%)", "0 / 7,801 (0%)")
)
write_tsv(da_summary, file.path(SRC_DIR, "Figure_2C_KO_DA_display_summary.tsv"))

# ---- Main Figure 2 ----------------------------------------------------------
surface_pcoa <- make_pcoa_panel(
  surface_pcoa_data, NULL, 0.479, "0.0375", "A", 0.14
)
surface_dbrda <- make_dbrda_panel(
  surface_sites, surface_pc_vectors,
  NULL,
  surface_stats$R2[[1]], surface_stats$adjusted_p_label[[1]], "B",
  c(0.757, 0.243)
)

da_colours <- c("Surface sites" = "#4A6A73", "Nahuelbuta depths" = pal_depth[["5-10 cm"]])
ko_da_panel <- ggplot(da_summary, aes(percent, comparison)) +
  geom_col(aes(fill = comparison), width = 0.54, colour = "grey25", linewidth = 0.28, show.legend = FALSE) +
  geom_text(
    aes(x = if_else(percent == 0, 0.55, percent + 0.55), label = label),
    hjust = 0, family = FONT, size = 3.70, colour = "black"
  ) +
  scale_fill_manual(values = da_colours) +
  scale_x_continuous(
    breaks = c(0, 5, 10, 15, 20, 25), limits = c(0, 27),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = "C", title = NULL,
    x = "FDR-supported KOs (% of tested)", y = NULL
  ) +
  theme_master() +
  theme(
    axis.text.y = element_text(size = 10.5, face = "bold"),
    axis.title.x = element_text(size = 11.0),
    panel.grid = element_blank()
  )

surface_legend <- make_site_legend()

render_pdf <- function(plot, path, width_cm, height_cm, bg = "transparent") {
  grDevices::cairo_pdf(
    path, width = width_cm / 2.54, height = height_cm / 2.54,
    family = FONT, bg = bg
  )
  print(plot)
  dev.off()
}

fig2_width_cm <- 18.3
fig2_height_cm <- 15.70
fig2_gutter_cm <- 0.42
fig2_top_height_cm <- 9.55
fig2_legend_height_cm <- 1.45
fig2_legend_width_cm <- 9.60
fig2_legend_gap_cm <- 0.10
fig2_row_gap_cm <- 0.18
fig2_bottom_height_cm <- 4.25
fig2_panel_width_cm <- (fig2_width_cm - fig2_gutter_cm) / 2
fig2_legend_y_cm <- fig2_bottom_height_cm + fig2_legend_gap_cm
fig2_top_y_cm <- fig2_legend_y_cm + fig2_legend_height_cm + fig2_row_gap_cm
fig2_right_x_cm <- fig2_panel_width_cm + fig2_gutter_cm

render_pdf(surface_pcoa, file.path(PANEL_DIR, "Figure_2_A.pdf"), fig2_panel_width_cm, fig2_top_height_cm)
render_pdf(surface_dbrda, file.path(PANEL_DIR, "Figure_2_B.pdf"), fig2_panel_width_cm, fig2_top_height_cm)
render_pdf(ko_da_panel, file.path(PANEL_DIR, "Figure_2_C.pdf"), fig2_width_cm, fig2_bottom_height_cm)
render_pdf(
  surface_legend, file.path(PANEL_DIR, "Figure_2_shared_legend.pdf"),
  fig2_legend_width_cm, fig2_legend_height_cm
)

figure_2 <- ggdraw() +
  theme(plot.background = element_rect(fill = "white", colour = NA)) +
  draw_grob(rectGrob(gp = gpar(fill = "white", col = NA)), x = 0, y = 0, width = 1, height = 1) +
  draw_plot(
    surface_pcoa, x = 0, y = fig2_top_y_cm / fig2_height_cm,
    width = fig2_panel_width_cm / fig2_width_cm,
    height = fig2_top_height_cm / fig2_height_cm
  ) +
  draw_plot(
    surface_dbrda, x = fig2_right_x_cm / fig2_width_cm,
    y = fig2_top_y_cm / fig2_height_cm,
    width = fig2_panel_width_cm / fig2_width_cm,
    height = fig2_top_height_cm / fig2_height_cm
  ) +
  draw_plot(
    ko_da_panel, x = 0, y = 0, width = 1,
    height = fig2_bottom_height_cm / fig2_height_cm
  ) +
  # The shared legend is composited last on a transparent child canvas so a
  # panel background can never obscure its text or the Panel C title.
  draw_plot(
    surface_legend,
    x = ((fig2_width_cm - fig2_legend_width_cm) / 2) / fig2_width_cm,
    y = fig2_legend_y_cm / fig2_height_cm,
    width = fig2_legend_width_cm / fig2_width_cm,
    height = fig2_legend_height_cm / fig2_height_cm
  )
render_pdf(figure_2, file.path(FIG_DIR, "Figure_2.pdf"), fig2_width_cm, fig2_height_cm, "white")

# ---- Supplementary Figure S1 -----------------------------------------------
loadings <- read_tsv(INPUTS[["soil_loadings"]], show_col_types = FALSE)
variance <- read_tsv(INPUTS[["soil_variance"]], show_col_types = FALSE)
if (!identical(loadings$variable, ENV_VARS)) {
  stop("Frozen soil-loading variable order changed.")
}
variance_lookup <- setNames(100 * variance$variance_fraction[1:2], c("PC1", "PC2"))

make_loading_panel <- function(axis_name, tag) {
  dat <- loadings %>% transmute(
    variable_label = factor(
      unname(PRETTY_ENV[variable]),
      levels = rev(unname(PRETTY_ENV[ENV_VARS]))
    ),
    loading = .data[[axis_name]],
    direction = factor(
      if_else(loading >= 0, "Positive", "Negative"),
      levels = c("Negative", "Positive")
    )
  )
  ggplot(dat, aes(loading, variable_label, fill = direction)) +
    geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.28) +
    geom_col(width = 0.68, colour = "white", linewidth = 0.18) +
    scale_fill_manual(values = c(Negative = pal_div_low, Positive = pal_div_high), guide = "none") +
    scale_x_continuous(
      breaks = c(-0.5, -0.25, 0, 0.25, 0.5), limits = c(-0.58, 0.58),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      tag = tag,
      x = sprintf("Variable loading on Soil %s", axis_name), y = NULL
    ) +
    theme_master(10.2) +
    theme(
      plot.title = element_text(size = 11.3, face = "bold", margin = margin(l = 28, b = 4)),
      axis.text.y = element_text(size = 9.6),
      axis.text.x = element_text(size = 9.6),
      axis.title.x = element_text(size = 10.7)
    )
}

loading_a <- make_loading_panel("PC1", "A")
loading_b <- make_loading_panel("PC2", "B")
depth_pcoa <- make_pcoa_panel(
  depth_pcoa_data, NULL, 0.389, "0.361", "C", 0.20,
  show_pit_labels = TRUE
)
depth_dbrda <- make_dbrda_panel(
  depth_sites, depth_pc_vectors,
  NULL,
  depth_stats$R2[[1]], depth_stats$adjusted_p_label[[1]], "D",
  c(0.809, 0.191)
)

depth_legend <- make_depth_legend()

s1_width_cm <- 18.3
s1_height_cm <- 20.50
s1_gutter_cm <- 0.32
s1_legend_height_cm <- 2.65
s1_legend_width_cm <- 5.40
s1_row_gap_cm <- 0.32
s1_panel_width_cm <- (s1_width_cm - s1_gutter_cm) / 2
s1_panel_height_cm <- (s1_height_cm - s1_legend_height_cm - s1_row_gap_cm) / 2
s1_top_y_cm <- s1_legend_height_cm + s1_row_gap_cm + s1_panel_height_cm
s1_right_x_cm <- s1_panel_width_cm + s1_gutter_cm

s1_panels <- list(A = loading_a, B = loading_b, C = depth_pcoa, D = depth_dbrda)
for (id in names(s1_panels)) {
  render_pdf(
    s1_panels[[id]], file.path(PANEL_DIR, paste0("Figure_S1_", id, ".pdf")),
    s1_panel_width_cm, s1_panel_height_cm
  )
}
render_pdf(
  depth_legend, file.path(PANEL_DIR, "Figure_S1_shared_legend.pdf"),
  s1_legend_width_cm, s1_legend_height_cm
)

figure_s1 <- ggdraw() +
  theme(plot.background = element_rect(fill = "white", colour = NA)) +
  draw_grob(rectGrob(gp = gpar(fill = "white", col = NA)), x = 0, y = 0, width = 1, height = 1) +
  draw_plot(
    loading_a, x = 0, y = s1_top_y_cm / s1_height_cm,
    width = s1_panel_width_cm / s1_width_cm,
    height = s1_panel_height_cm / s1_height_cm
  ) +
  draw_plot(
    loading_b, x = s1_right_x_cm / s1_width_cm,
    y = s1_top_y_cm / s1_height_cm,
    width = s1_panel_width_cm / s1_width_cm,
    height = s1_panel_height_cm / s1_height_cm
  ) +
  draw_plot(
    depth_pcoa, x = 0, y = s1_legend_height_cm / s1_height_cm,
    width = s1_panel_width_cm / s1_width_cm,
    height = s1_panel_height_cm / s1_height_cm
  ) +
  draw_plot(
    depth_dbrda, x = s1_right_x_cm / s1_width_cm,
    y = s1_legend_height_cm / s1_height_cm,
    width = s1_panel_width_cm / s1_width_cm,
    height = s1_panel_height_cm / s1_height_cm
  ) +
  draw_plot(
    depth_legend,
    x = ((s1_width_cm - s1_legend_width_cm) / 2) / s1_width_cm,
    y = 0, width = s1_legend_width_cm / s1_width_cm,
    height = s1_legend_height_cm / s1_height_cm
  )
render_pdf(figure_s1, file.path(FIG_DIR, "Figure_S1.pdf"), s1_width_cm, s1_height_cm, "white")

# ---- Mechanical geometry ledgers ------------------------------------------
pt_per_cm <- 72 / 2.54
tmp_dir <- tempfile("v0_figure_geometry_")
dir.create(tmp_dir)
on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

measure_ink <- function(pdf, width_pt, height_pt) {
  prefix <- file.path(tmp_dir, tools::file_path_sans_ext(basename(pdf)))
  status <- system2(
    "pdftoppm",
    c("-png", "-singlefile", "-r", "300", shQuote(pdf), shQuote(prefix)),
    stdout = FALSE, stderr = FALSE
  )
  if (!identical(status, 0L)) {
    stop("Could not rasterize source panel for geometry audit: ", pdf)
  }
  image <- readPNG(paste0(prefix, ".png"))
  rgb <- image[, , seq_len(min(3, dim(image)[3])), drop = FALSE]
  ink <- rgb[, , 1] < 0.985 | rgb[, , 2] < 0.985 | rgb[, , 3] < 0.985
  rows <- which(rowSums(ink) > 0)
  cols <- which(colSums(ink) > 0)
  c(
    (min(cols) - 1) / ncol(ink) * width_pt,
    (nrow(ink) - max(rows)) / nrow(ink) * height_pt,
    max(cols) / ncol(ink) * width_pt,
    (nrow(ink) - min(rows) + 1) / nrow(ink) * height_pt
  )
}

panel_record <- function(id, source_pdf, source_width_cm, source_height_cm, left_cm, bottom_cm) {
  placed <- c(
    left_cm, bottom_cm,
    left_cm + source_width_cm, bottom_cm + source_height_cm
  ) * pt_per_cm
  ink <- measure_ink(
    source_pdf,
    source_width_cm * pt_per_cm,
    source_height_cm * pt_per_cm
  )
  list(
    id = id, role = "panel",
    source_width_pt = source_width_cm * pt_per_cm,
    source_height_pt = source_height_cm * pt_per_cm,
    placed_box_pt = unname(placed),
    ink_box_pt = unname(c(
      placed[1] + ink[1], placed[2] + ink[2],
      placed[1] + ink[3], placed[2] + ink[4]
    ))
  )
}

fig2_ledger <- list(
  figure_id = "readability_first_Figure_2_three_panel_KO_evidence_chain",
  aspect_tolerance_pct = 0.05,
  source_panels_locked_before_composition = TRUE,
  uses_outer_boxes_for_alignment = FALSE,
  page_size_cm = c(width = fig2_width_cm, height = fig2_height_cm),
  panels = list(
    panel_record("A", file.path(PANEL_DIR, "Figure_2_A.pdf"), fig2_panel_width_cm, fig2_top_height_cm, 0, fig2_top_y_cm),
    panel_record("B", file.path(PANEL_DIR, "Figure_2_B.pdf"), fig2_panel_width_cm, fig2_top_height_cm, fig2_right_x_cm, fig2_top_y_cm),
    panel_record("C", file.path(PANEL_DIR, "Figure_2_C.pdf"), fig2_width_cm, fig2_bottom_height_cm, 0, 0)
  ),
  gaps = list(
    list(panel_a = "A", panel_b = "B", axis = "x", kind = "top-row optical gutter", min_pt = 3, max_pt = 80),
    list(panel_a = "A", panel_b = "C", axis = "y", kind = "left evidence-step gutter", min_pt = 4, max_pt = 110),
    list(panel_a = "B", panel_b = "C", axis = "y", kind = "right evidence-step gutter", min_pt = 4, max_pt = 110)
  )
)
write_json(fig2_ledger, file.path(OUT, "geometry_ledger_Figure_2.json"), pretty = TRUE, auto_unbox = TRUE)

s1_ledger <- list(
  figure_id = "readability_first_Figure_S1_soil_axes_and_depth_KO_context",
  aspect_tolerance_pct = 0.05,
  source_panels_locked_before_composition = TRUE,
  uses_outer_boxes_for_alignment = FALSE,
  page_size_cm = c(width = s1_width_cm, height = s1_height_cm),
  panels = list(
    panel_record("A", file.path(PANEL_DIR, "Figure_S1_A.pdf"), s1_panel_width_cm, s1_panel_height_cm, 0, s1_top_y_cm),
    panel_record("B", file.path(PANEL_DIR, "Figure_S1_B.pdf"), s1_panel_width_cm, s1_panel_height_cm, s1_right_x_cm, s1_top_y_cm),
    panel_record("C", file.path(PANEL_DIR, "Figure_S1_C.pdf"), s1_panel_width_cm, s1_panel_height_cm, 0, s1_legend_height_cm),
    panel_record("D", file.path(PANEL_DIR, "Figure_S1_D.pdf"), s1_panel_width_cm, s1_panel_height_cm, s1_right_x_cm, s1_legend_height_cm)
  ),
  gaps = list(
    list(panel_a = "A", panel_b = "B", axis = "x", kind = "top-row optical gutter", min_pt = 3, max_pt = 90),
    list(panel_a = "C", panel_b = "D", axis = "x", kind = "bottom-row optical gutter", min_pt = 3, max_pt = 90),
    list(panel_a = "A", panel_b = "C", axis = "y", kind = "left-column optical gutter", min_pt = 3, max_pt = 100),
    list(panel_a = "B", panel_b = "D", axis = "y", kind = "right-column optical gutter", min_pt = 3, max_pt = 100)
  )
)
write_json(s1_ledger, file.path(OUT, "geometry_ledger_Figure_S1.json"), pretty = TRUE, auto_unbox = TRUE)

write_json(list(
  figure_only_candidate = TRUE,
  manuscript_modified = FALSE,
  statistical_inference_rerun = FALSE,
  ordination_coordinates_changed = FALSE,
  figure_2 = list(
    A = list(role = "surface KO PCoA", R2 = 0.479, adjusted_p = 0.0375),
    B = list(role = "surface KO dbRDA constrained by Soil PC1 and Soil PC2", R2 = surface_stats$R2[[1]], adjusted_p = surface_stats$adjusted_p[[1]]),
    C = list(role = "FDR-supported KO differentiation breadth", surface = "1403/7841 (17.9%)", depth = "0/7801 (0%)")
  ),
  figure_s1 = list(
    A = "Soil PC1 loadings", B = "Soil PC2 loadings",
    C = list(role = "depth KO PCoA", R2 = 0.389, adjusted_p = 0.361),
    D = list(role = "depth KO dbRDA constrained by Soil PC1 and Soil PC2", R2 = depth_stats$R2[[1]], adjusted_p = depth_stats$adjusted_p[[1]])
  ),
  soil_pc_vector_method = "envfit directions from frozen dbRDA site scores and frozen Soil PC scores; permutations=0",
  generated = format(Sys.time(), tz = "Asia/Tokyo")
), file.path(OUT, "READABILITY_FIGURE2_S1_QA.json"), pretty = TRUE, auto_unbox = TRUE)

write_lines(c(
  "# Readability-first Figure 2 and Figure S1",
  "",
  "- Figure 2A: frozen surface KO PCoA.",
  "- Figure 2B: frozen surface KO dbRDA coordinates with only Soil PC1 and Soil PC2 shown as descriptive predictor orientations.",
  "- Figure 2C: publication-facing breadth of FDR-supported KO differentiation (surface 1,403/7,841; depth 0/7,801).",
  "- Figure S1A-B: frozen Soil PC loadings.",
  "- Figure S1C-D: frozen Nahuelbuta depth KO PCoA/dbRDA; C labels each point by soil pit, and D shows only Soil PC1 and Soil PC2.",
  "- No manuscript, statistic, ordination coordinate, sample assignment, or scientific conclusion was changed.",
  "- Figure S2 is handled by the dedicated readability-first phylogeny renderer."
), file.path(OUT, "README.md"))

if (length(list.files(OUT, pattern = "\\.png$", recursive = TRUE))) {
  stop("Unexpected persistent PNG output.")
}
if (any(file.info(file.path(FIG_DIR, c("Figure_2.pdf", "Figure_S1.pdf")))$size <= 0)) {
  stop("Candidate PDF export failed.")
}

message("Figure-only candidate written to: ", FIG_DIR)
