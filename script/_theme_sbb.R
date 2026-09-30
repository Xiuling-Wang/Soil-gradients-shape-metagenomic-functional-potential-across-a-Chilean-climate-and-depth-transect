# =============================================================================
# Shared Rhizosphere-style design system for the Chile metagenomic figures.
#   source() this at the top of each figure script, then use theme_sbb() and the
#   shared palettes so all main and supplementary figures read as one set.
#   The Chile site/depth palettes are deliberate project-specific deviations.
# =============================================================================
# The non-interactive render environment can start in the C locale even when
# UTF-8 is requested by the shell. Set the character locale before any labels
# are constructed so chemical subscripts/superscripts stay in Helvetica rather
# than triggering device-specific fallback fonts.
invisible(suppressWarnings(try(
  Sys.setlocale("LC_CTYPE", "en_US.UTF-8"),
  silent = TRUE
)))

suppressPackageStartupMessages(library(ggplot2))

SBB_FONT <- "Helvetica"
SBB_HEATMAP_BORDER <- "white"
SBB_BAR_OUTLINE <- "grey70"
SBB_BAR_WIDTH <- 0.70
SBB_GENERIC_BLUE <- "#3D6F97"
SBB_GENERIC_RED <- "#B05246"
SBB_TEXT <- "black"
SBB_TEXT_MUTED <- "grey30"

# ---- Palettes ---------------------------------------------------------------
# Site palette, ordered arid -> humid (warm earth -> cool water), desaturated.
pal_site <- c(
  AZ   = "#C25B3F",  # arid        terracotta
  SG   = "#D6A141",  # semi-arid   ochre
  LC   = "#6FA055",  # Mediterran. olive-green
  NB   = "#3D6F97",  # humid       slate blue
  `NA` = "#3D6F97"   # alias (NB displayed as NA)
)
# Sequential and diverging ramps follow the Rhizosphere implementation.
pal_seq <- RColorBrewer::brewer.pal(9, "Blues")[c(2, 4, 6, 8, 9)]
pal_div <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101)
pal_div_low <- pal_div[[1]]
pal_div_mid <- pal_div[[51]]
pal_div_high <- pal_div[[101]]

# Depth and six-group palettes. Site colour remains the primary encoding;
# Nahuelbuta depths use a sequential version of the same humid-site blue.
pal_depth <- c(
  `5-10 cm` = "#8DB5CD",
  `20-40 cm` = "#5F88AA",
  `60-80 cm` = "#2E5577"
)
pal_group <- c(
  AZ = pal_site[["AZ"]],
  SG = pal_site[["SG"]],
  LC = pal_site[["LC"]],
  `NA-S` = pal_depth[["5-10 cm"]],
  `NA-M` = pal_depth[["20-40 cm"]],
  `NA-D` = pal_depth[["60-80 cm"]],
  `NA-surf` = pal_depth[["5-10 cm"]],
  `NA-mid` = pal_depth[["20-40 cm"]],
  `NA-deep` = pal_depth[["60-80 cm"]]
)

# Stable Rhizosphere/ITS Paired palette used by every MAG bar, point, and tree.
# The pale-yellow eleventh swatch is omitted because it disappears on white;
# the catch-all category stays neutral grey.
phylum_palette <- c(
  setNames(
    RColorBrewer::brewer.pal(12, "Paired")[-11][seq_len(9)],
    c(
      "Actinomycetota", "Dormibacterota", "Pseudomonadota",
      "Acidobacteriota", "Gemmatimonadota", "Chloroflexota",
      "Thermoproteota", "Thermoplasmatota", "Myxococcota"
    )
  ),
  Other = "grey70"
)
# Preserve the high-separation archaeal pair used in the earlier Chile figures.
phylum_palette[c("Thermoproteota", "Thermoplasmatota")] <- c(
  Thermoproteota = "#4DBBD5",
  Thermoplasmatota = "#E64B35"
)

# Journal-facing labels: Environmental Microbiome asks that taxonomic names at
# phylum/class/genus rank are italicised in both text and figures. Keep neutral
# catch-all labels roman while returning plotmath expressions for ggplot scales.
taxon_plotmath_labels <- function(levels, roman = c("Other", "Unassigned", "Unclassified")) {
  as.expression(lapply(as.character(levels), function(label) {
    if (label %in% roman) bquote(.(label)) else bquote(italic(.(label)))
  }))
}

# Secondary categorical accents. These are deliberately restrained and should
# not replace the site, depth, or phylum palettes when those identities apply.
pal_category <- c(
  blue = "#4F7898",
  teal = "#5C9482",
  olive = "#7C9A61",
  ochre = "#C79A48",
  terracotta = "#B66A4A",
  mauve = "#8B719E",
  grey = "#AAA8A1"
)

# Convenience colourbar/scale builders
scale_fill_sbb_div <- function(..., midpoint = 0)
  scale_fill_gradient2(low = pal_div[[1]], mid = pal_div[[51]], high = pal_div[[101]],
                       midpoint = midpoint, ...)
scale_fill_sbb_seq <- function(...)
  scale_fill_gradientn(colours = pal_seq, ...)

# ComplexHeatmap colour functions. General matrix heatmaps throughout the
# manuscript use ComplexHeatmap with these exact ramps; correlation, z-score,
# DRAM, and curated-marker panels therefore no longer drift in hue.
sbb_div_col_fun <- function(limits = c(-2, 2)) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("The circlize package is required for ComplexHeatmap colour functions.")
  }
  limits <- range(limits, finite = TRUE)
  circlize::colorRamp2(seq(limits[[1]], limits[[2]], length.out = 101), pal_div)
}

sbb_seq_col_fun <- function(limits = c(0, 1)) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("The circlize package is required for ComplexHeatmap colour functions.")
  }
  limits <- range(limits, finite = TRUE)
  circlize::colorRamp2(
    seq(limits[[1]], limits[[2]], length.out = 101),
    colorRampPalette(pal_seq)(101)
  )
}

sbb_heatmap_grob <- function(
  heatmap,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = grid::unit(c(12, 3, 4, 3), "mm")
) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("The ComplexHeatmap package is required for unified matrix figures.")
  }
  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      heatmap,
      heatmap_legend_side = heatmap_legend_side,
      annotation_legend_side = annotation_legend_side,
      padding = padding,
      merge_legend = TRUE
    ),
    wrap.grobs = TRUE
  )
}

sbb_group_annotation <- function(groups, name = "Group", show_legend = FALSE) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("The ComplexHeatmap package is required for unified matrix figures.")
  }
  groups <- as.character(groups)
  missing <- setdiff(unique(groups), names(pal_group))
  if (length(missing)) stop("Missing group colours: ", paste(missing, collapse = ", "))
  values <- factor(groups, levels = unique(groups))
  cols <- pal_group[levels(values)]
  ComplexHeatmap::HeatmapAnnotation(
    Group = values,
    col = list(Group = cols),
    annotation_name_side = "left",
    annotation_name_gp = grid::gpar(fontfamily = SBB_FONT, fontsize = 8, fontface = "bold"),
    simple_anno_size = grid::unit(3.0, "mm"),
    show_legend = show_legend,
    show_annotation_name = FALSE
  )
}

# ---- Theme ------------------------------------------------------------------
# Rhizosphere house style: classic axes, white background, no grid or frame.
theme_sbb <- function(base_size = 8, base_family = SBB_FONT) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text              = element_text(family = base_family, colour = "black"),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      panel.grid        = element_blank(),
      panel.border      = element_blank(),
      axis.line         = element_line(colour = "black", linewidth = 0.35),
      axis.ticks        = element_line(colour = "black", linewidth = 0.3),
      axis.text         = element_text(colour = "black", size = base_size),
      axis.title        = element_text(colour = "black", size = base_size),
      strip.background  = element_blank(),
      strip.text        = element_text(
        face = "bold", colour = "black", size = max(8.5, base_size)
      ),
      legend.key        = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.text       = element_text(colour = "black", size = min(6.8, base_size)),
      legend.title      = element_text(
        colour = "black", size = 7.2, face = "bold"
      ),
      legend.key.size   = unit(0.36, "cm"),
      plot.title        = element_text(
        colour = "black", size = max(9, base_size + 1), face = "bold",
        hjust = 0, margin = margin(b = 3)
      ),
      plot.subtitle     = element_text(
        colour = "black", size = min(7.2, base_size), hjust = 0,
        margin = margin(b = 4)
      ),
      plot.margin       = margin(5, 6, 5, 6)
    )
}

theme_sbb_bar <- function(base_size = 8, base_family = SBB_FONT) {
  theme_sbb(base_size = base_size, base_family = base_family) +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}
