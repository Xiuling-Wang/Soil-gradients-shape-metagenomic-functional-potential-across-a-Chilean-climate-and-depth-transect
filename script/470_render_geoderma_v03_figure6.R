# Geoderma V03 derivative of 449_render_v26_figure5.R : renumbered Fig. 5 -> Fig. 6; page 23.5 x 29.25 -> 19.0 x 23.65 cm (x 0.8085) with every text size, legend key and point size scaled by the same factor; data, tree, colours, labels and italic/roman styling unchanged.
# V26 derivative of 436_render_v24_figure5.R (443_build_v26_figures_no_panel_titles.py): panel titles and the B subtitle removed (the 'Other contig' note is in the caption).
# V24 derivative of 391_render_readability_figure5.R (built by 438_build_v24_figure_refinement.py). Support labels moved off branch lines; tree-scale text in Arial.
# =============================================================================
# Readability-first main Figure 5 refinement. Methods and
# biological interpretation are kept in the manuscript legend and Results.
#
# Page 3 is appended to the currently accepted two-page Figure S3. It validates
# DRAM K10944 homolog identities against named UniProt archaeal AmoA, bacterial
# AmoA, and PmoA references, and shows the contig-context caveat for the two
# Nitrososphaeria MAGs. Existing Figure S3 pages are not modified.
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(ggtext)
  library(ggtree)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(tidyr)
})

invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)))

PROJECT <- Sys.getenv("CHILE_PROJECT_ROOT", unset=getwd())
OUT_ROOT <- Sys.getenv(
  "CHILE_READABILITY_FIGURE5_ROOT",
  file.path(PROJECT, "analysis/readability_first_figure_refinement_2026-08-30")
)
INPUT_ROOT <- Sys.getenv("CHILE_FIGURE5_INPUT_ROOT", OUT_ROOT)
PHYLO_DIR <- file.path(INPUT_ROOT, "phylogeny")
SOURCE_INPUT_DIR <- file.path(INPUT_ROOT, "source_data")
SOURCE_DIR <- file.path(OUT_ROOT, "source_data")
FIGURE_DIR <- file.path(OUT_ROOT, "figures")
DOCS_DIR <- file.path(OUT_ROOT, "docs")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DOCS_DIR, recursive = TRUE, showWarnings = FALSE)

CURRENT_S3 <- file.path(
  PROJECT,
  "EM（BMC）/00_CURRENT_MANUSCRIPT_V17/Supplementary_Figures_PDF/Figure_S3.pdf"
)
TREE_FILE <- file.path(PHYLO_DIR, "amoA_pmoA_ml.treefile")
IQTREE_FILE <- file.path(PHYLO_DIR, "amoA_pmoA_ml.iqtree")
METADATA_FILE <- file.path(SOURCE_INPUT_DIR, "amoA_pmoA_sequence_metadata.tsv")
CONTEXT_FILE <- file.path(SOURCE_INPUT_DIR, "amoA_local_gene_context.tsv")

stopifnot(
  file.exists(CURRENT_S3), file.exists(TREE_FILE), file.exists(IQTREE_FILE),
  file.exists(METADATA_FILE), file.exists(CONTEXT_FILE)
)

source(file.path(PROJECT, "script/_theme_sbb.R"))
SBB_FONT <- "Arial"
if (capabilities("aqua")) quartzFonts(
  Arial = quartzFont(c("Arial", "Arial Bold", "Arial Italic", "Arial Bold Italic"))
)
arial_pdf <- function(filename, width, height, bg = "white", ...) {
  quartz(
    type = "pdf", file = filename, width = width, height = height,
    family = "Arial", bg = bg
  )
}

# ---- Data contracts --------------------------------------------------------
tree_unrooted <- ape::read.tree(TREE_FILE)
metadata <- read_tsv(METADATA_FILE, show_col_types = FALSE)
context <- read_tsv(CONTEXT_FILE, show_col_types = FALSE)

stopifnot(
  length(tree_unrooted$tip.label) == 21L,
  nrow(metadata) == 21L,
  setequal(tree_unrooted$tip.label, metadata$sequence_id),
  sum(metadata$source == "local_MAG") == 5L,
  sum(metadata$reference_class == "Archaeal AmoA") == 6L,
  sum(metadata$reference_class == "Bacterial AmoA") == 4L,
  sum(metadata$reference_class == "PmoA") == 6L,
  setequal(metadata$MAG[metadata$source == "local_MAG"], c("MAG005", "MAG048", "MAG070", "MAG071", "MAG090"))
)

archaeal_tips <- metadata %>%
  filter(Domain == "Archaea") %>%
  pull(sequence_id)
bacterial_tips <- metadata %>%
  filter(Domain == "Bacteria") %>%
  pull(sequence_id)

stopifnot(
  ape::is.monophyletic(tree_unrooted, archaeal_tips),
  ape::is.monophyletic(tree_unrooted, bacterial_tips)
)

tree_rooted <- tree_unrooted %>%
  ape::root(outgroup = bacterial_tips, resolve.root = TRUE) %>%
  ape::ladderize(right = FALSE)

iqtree_lines <- readLines(IQTREE_FILE, warn = FALSE)
model_line <- iqtree_lines[str_detect(iqtree_lines, "^Best-fit model according to BIC:")]
stopifnot(length(model_line) == 1L)
best_model <- str_trim(str_remove(model_line, "^Best-fit model according to BIC:"))

# ---- Nearest-reference audit ----------------------------------------------
tree_distances <- ape::cophenetic.phylo(tree_unrooted)
query_ids <- metadata %>% filter(source == "local_MAG") %>% pull(sequence_id)
reference_ids <- metadata %>% filter(source == "UniProt") %>% pull(sequence_id)

nearest_reference <- bind_rows(lapply(query_ids, function(query_id) {
  nearest <- sort(tree_distances[query_id, reference_ids])
  nearest_ids <- names(nearest)[seq_len(min(3L, length(nearest)))]
  query_row <- metadata %>% filter(sequence_id == query_id)
  reference_rows <- metadata %>%
    filter(sequence_id %in% nearest_ids) %>%
    mutate(nearest_rank = match(sequence_id, nearest_ids)) %>%
    arrange(nearest_rank)
  tibble(
    query_id = query_id,
    MAG = query_row$MAG,
    Domain = query_row$Domain,
    query_phylum = query_row$phylum,
    query_class = query_row$class,
    nearest_rank = seq_along(nearest_ids),
    reference_id = nearest_ids,
    reference_class = reference_rows$reference_class,
    reference_accession = reference_rows$accession,
    reference_organism = reference_rows$organism,
    patristic_distance = as.numeric(nearest[nearest_ids])
  )
})) %>%
  group_by(query_id) %>%
  mutate(
    targeted_interpretation = case_when(
      first(Domain) == "Archaea" & first(reference_class) == "Archaeal AmoA" ~
        "Archaeal AmoA supported",
      TRUE ~ "Divergent K10944 homolog; AmoA/PmoA identity unresolved"
    )
  ) %>%
  ungroup()

stopifnot(
  all(
    nearest_reference %>%
      filter(MAG %in% c("MAG005", "MAG070"), nearest_rank == 1L) %>%
      pull(reference_class) == "Archaeal AmoA"
  )
)

write_tsv(
  nearest_reference,
  file.path(SOURCE_DIR, "FigureS3_amoA_nearest_reference_audit.tsv")
)

# ---- Tree plotting metadata ------------------------------------------------
reference_colours <- c(
  "Local archaeal K10944" = "#B2182B",
  "Archaeal AmoA" = "#1B9E77",
  "Bacterial AmoA" = "#2166AC",
  "PmoA" = "#D95F02",
  "Local bacterial K10944" = "#6B7280"
)
reference_levels <- names(reference_colours)

clean_organism <- function(x) {
  x %>%
    str_remove(" \\(strain.*$") %>%
    str_remove(" \\(isolate.*$") %>%
    str_replace("^Candidatus ", "Ca. ")
}

plot_metadata <- metadata %>%
  mutate(
    plot_class = factor(reference_class, levels = reference_levels),
    source_type = factor(
      if_else(source == "local_MAG", "Local MAG protein", "UniProt reference"),
      levels = c("Local MAG protein", "UniProt reference")
    ),
    local_locus = if_else(
      source == "local_MAG",
      str_remove(gene_id, paste0("^", MAG, "_")),
      ""
    ),
    display_label = if_else(
      source == "local_MAG",
      paste0(
        "<b>", MAG, "</b> ",
        "<span style='color:#6B7280'>(", local_locus, ")</span>"
      ),
      paste0(
        "<i>", clean_organism(organism), "</i> ",
        "<span style='color:#6B7280'>(", accession, ")</span>"
      )
    )
  ) %>%
  select(label = sequence_id, plot_class, source_type, display_label)

archaeal_mrca <- ape::getMRCA(tree_rooted, archaeal_tips)
archaeal_support <- tree_rooted$node.label[archaeal_mrca - ape::Ntip(tree_rooted)]
support_parts <- str_split_fixed(archaeal_support, "/", 2)
archaeal_shalrt <- as.numeric(support_parts[1])
archaeal_ufboot <- as.numeric(support_parts[2])
tree_panel_title <- Sys.getenv(
  "CHILE_FIGURE5_TREE_TITLE",
  "Targeted phylogenetic validation of K10944 homologs"
)
context_panel_title <- Sys.getenv(
  "CHILE_FIGURE5_CONTEXT_TITLE",
  "Genome-context evidence for archaeal AmoA assignments"
)
support_label <- sprintf(
  "SH-aLRT/UFBoot %.1f/%.1f", archaeal_shalrt, archaeal_ufboot
)

p_tree <- ggtree(
  tree_rooted, layout = "rectangular", size = 0.52, colour = "grey32"
)
p_tree$data$ufboot <- suppressWarnings(as.numeric(str_extract(p_tree$data$label, "(?<=/)[0-9.]+")))
p_tree <- p_tree %<+% plot_metadata

tree_xmax <- max(p_tree$data$x, na.rm = TRUE)
tree_support_colour <- Sys.getenv("CHILE_TREE_SUPPORT_COLOUR", "grey28")
tree_support_size <- as.numeric(Sys.getenv("CHILE_TREE_SUPPORT_SIZE", "3.31"))
tree_tip_label_size <- as.numeric(Sys.getenv("CHILE_TREE_TIP_LABEL_SIZE", "3.44"))
tree_xlim_multiplier <- as.numeric(Sys.getenv("CHILE_TREE_XLIM_MULTIPLIER", "2.20"))
if (any(!is.finite(c(tree_support_size, tree_tip_label_size, tree_xlim_multiplier))) ||
    any(c(tree_support_size, tree_tip_label_size, tree_xlim_multiplier) <= 0)) {
  stop("Tree rendering controls must be finite positive numbers.")
}

p_tree <- p_tree +
  geom_hilight(
    node = archaeal_mrca, fill = "#E5F3EE", alpha = 0.68,
    extend = 0.08 * tree_xmax
  ) +
  geom_tippoint(
    aes(fill = plot_class, shape = source_type, size = source_type),
    colour = "grey20", stroke = 0.34
  ) +
  ggtext::geom_richtext(
    data = function(x) filter(x, isTip),
    aes(label = display_label),
    hjust = 0, vjust = 0.5, nudge_x = 0.035 * tree_xmax,
    family = SBB_FONT, size = tree_tip_label_size,
    fill = NA, label.color = NA,
    label.padding = unit(0, "pt"), label.r = unit(0, "pt")
  ) +
  geom_text2(
    aes(
      subset = !isTip & node != archaeal_mrca & !is.na(ufboot) & ufboot >= 75,
      label = round(ufboot)
    ),
    nudge_x = -0.030 * tree_xmax, nudge_y = -0.18, hjust = 1,
    colour = tree_support_colour,
    family = SBB_FONT, size = tree_support_size
  ) +
  geom_label2(
    aes(
      subset = !isTip & node == archaeal_mrca,
      label = support_label
    ),
    nudge_x = -0.08 * tree_xmax, nudge_y = 1.08,
    hjust = 1, family = SBB_FONT, size = 3.40,
    colour = "black", fill = "white", label.size = 0,
    label.padding = unit(0.8, "mm")
  ) +
  scale_fill_manual(
    name = "Sequence class", values = reference_colours,
    breaks = reference_levels, drop = FALSE
  ) +
  scale_shape_manual(
    values = c("Local MAG protein" = 23, "UniProt reference" = 21),
    guide = "none"
  ) +
  scale_size_manual(
    values = c("Local MAG protein" = 3.11, "UniProt reference" = 2.43),
    guide = "none"
  ) +
  geom_treescale(
    x = 0.08, y = 0.40, width = 0.50,
    fontsize = 3.40, linesize = 0.58, family = SBB_FONT
  ) +
  xlim(0, tree_xmax * tree_xlim_multiplier) +
  labs(
    tag = "A",
    title = NULL,
    x = "Branch length (amino-acid substitutions per site)",
    subtitle = NULL,
    caption = "Scale: amino-acid substitutions per site; focal support shown as SH-aLRT/UFBoot"
  ) +
  theme_tree2() +
  theme(
    text = element_text(family = SBB_FONT, colour = "black"),
    plot.tag = element_text(size = 16.2, face = "bold"),
    plot.tag.position = c(0.005, 0.995),
    plot.title.position = "plot",
    plot.title = element_text(size = 15.2, face = "bold", margin = margin(l = 28, b = 4)),
    plot.subtitle = element_text(size = 6.8, colour = "grey25", margin = margin(b = 4)),
    plot.caption = element_text(size = 10.1, colour = "black", hjust = 0.5, margin = margin(t = 5)),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.title = element_text(size = 10.35, face = "bold"),
    legend.text = element_text(size = 9.54),
    legend.key.width = unit(5.26, "mm"),
    plot.margin = margin(4, 6, 2, 2)
  ) +
  guides(
    fill = guide_legend(
      order = 1, nrow = 2, byrow = TRUE,
      override.aes = list(
        shape = c(23, 21, 21, 21, 23),
        size = c(3.11, 2.43, 2.43, 2.43, 3.11),
        colour = "grey20"
      )
    )
  )

# ---- MAG005/MAG070 context caveat -----------------------------------------
target_mags <- c("MAG005", "MAG070")
target_metadata <- metadata %>%
  filter(MAG %in% target_mags) %>%
  select(MAG, pfam_hits)

context_matrix <- expand_grid(
  MAG = factor(target_mags, levels = rev(target_mags)),
  evidence = factor(
    c("K10944\nAmoA", "PF12942\narchaeal AmoA", "K10945\nAmoB", "K10946\nAmoC"),
    levels = c("K10944\nAmoA", "PF12942\narchaeal AmoA", "K10945\nAmoB", "K10946\nAmoC")
  )
) %>%
  rowwise() %>%
  mutate(
    state = case_when(
      as.character(evidence) == "K10944\nAmoA" ~ "Target/domain support",
      as.character(evidence) == "PF12942\narchaeal AmoA" &
        str_detect(target_metadata$pfam_hits[target_metadata$MAG == as.character(MAG)], "PF12942") ~
        "Target/domain support",
      as.character(evidence) == "PF12942\narchaeal AmoA" ~ "Not detected",
      as.character(evidence) == "K10945\nAmoB" &
        any(context$MAG == as.character(MAG) & context$KO == "K10945") ~
        "Detected on another contig",
      as.character(evidence) == "K10946\nAmoC" &
        any(context$MAG == as.character(MAG) & context$KO == "K10946") ~
        "Detected on another contig",
      TRUE ~ "Not detected"
    )
  ) %>%
  ungroup() %>%
  mutate(
    state = factor(
      state,
      levels = c("Target/domain support", "Detected on another contig", "Not detected")
    )
  )

write_tsv(
  context_matrix %>%
    mutate(
      across(where(is.factor), as.character),
      evidence = str_replace_all(evidence, "\\n", " / ")
    ),
  file.path(SOURCE_DIR, "FigureS3_amoA_context_matrix.tsv")
)

context_colours <- c(
  "Target/domain support" = "#DCEFE7",
  "Detected on another contig" = "#F5E7B2",
  "Not detected" = "#F2F2F2"
)

context_text_colours <- c(
  "Target/domain support" = "#14532D",
  "Detected on another contig" = "#6B4F00",
  "Not detected" = "#595959"
)

context_display <- context_matrix %>%
  mutate(
    evidence = factor(
      as.character(evidence),
      levels = c(
        "K10944\nAmoA", "PF12942\narchaeal AmoA",
        "K10945\nAmoB", "K10946\nAmoC"
      ),
      labels = c(
        "K10944 target\nphylogeny", "PF12942\nAmoA domain",
        "K10945\nAmoB", "K10946\nAmoC"
      )
    ),
    cell_label = case_when(
      as.character(evidence) == "K10944 target\nphylogeny" ~ "Tree-supported",
      as.character(state) == "Target/domain support" ~ "Detected",
      as.character(state) == "Detected on another contig" ~ "Other contig",
      TRUE ~ "Not detected"
    ),
    text_colour = unname(context_text_colours[as.character(state)])
  )

p_context <- ggplot(context_display, aes(evidence, MAG)) +
  geom_tile(
    aes(fill = state), colour = "grey38", linewidth = 0.48,
    width = 0.98, height = 0.92
  ) +
  geom_text(
    aes(label = cell_label, colour = text_colour),
    family = SBB_FONT, size = 3.64, fontface = "bold", lineheight = 0.95
  ) +
  scale_fill_manual(values = context_colours, drop = FALSE, guide = "none") +
  scale_colour_identity() +
  scale_x_discrete(position = "top") +
  labs(
    tag = "B", x = NULL, y = NULL,
    title = NULL,
    subtitle = NULL
  ) +
  theme_sbb(base_size = 8.2) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x.top = element_text(size = 10.67, face = "bold", lineheight = 0.95, colour = "black"),
    axis.text.y = element_text(size = 11.32, face = "bold", colour = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.tag = element_text(size = 16.2, face = "bold"),
    plot.tag.position = c(0.005, 0.995),
    plot.title.position = "plot",
    plot.title = element_text(size = 15.2, face = "bold", margin = margin(l = 28, b = 3)),
    plot.subtitle = element_text(size = 12.2, colour = "grey25", lineheight = 1.0, margin = margin(l = 28, b = 6)),
    legend.position = "none",
    plot.margin = margin(4, 8, 5, 2)
  )

page_3 <- wrap_elements(full = p_tree) /
  plot_spacer() /
  wrap_elements(full = p_context) +
  plot_layout(heights = c(20.5, 0.85, 7.9))

native_panel_dir <- Sys.getenv("CHILE_NATIVE_PANEL_DIR", "")
if (nzchar(native_panel_dir)) {
  dir.create(native_panel_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(native_panel_dir, "Figure_6_A.pdf"), p_tree,
         width = 19.0, height = 16.57, units = "cm", device = arial_pdf, bg = "white")
  ggsave(file.path(native_panel_dir, "Figure_6_B.pdf"), p_context,
         width = 19.0, height = 6.39, units = "cm", device = arial_pdf, bg = "white")
}

standalone_pdf <- file.path(FIGURE_DIR, "Figure_6.pdf")
combined_pdf <- file.path(
  FIGURE_DIR, "Figure_S3_candidate_legacy_layout_amoA_validated.pdf"
)

ggsave(
  standalone_pdf, page_3,
  width = 19.0, height = 23.65, units = "cm",
  device = arial_pdf, bg = "white"
)

qa_png <- Sys.getenv("CHILE_QA_PNG", "")
if (nzchar(qa_png)) {
  dir.create(dirname(qa_png), recursive = TRUE, showWarnings = FALSE)
  ragg::agg_png(
    qa_png, width = 19.0, height = 23.65, units = "cm",
    res = 180, background = "white"
  )
  print(page_3)
  dev.off()
}

pt_per_cm <- 72 / 2.54
jsonlite::write_json(
  list(
    figure_id = "geoderma_v03_Figure_6_targeted_amoA_phylogeny",
    aspect_tolerance_pct = 0.05,
    source_panels_locked_before_composition = TRUE,
    uses_outer_boxes_for_alignment = FALSE,
    page_size_cm = c(width = 19.0, height = 23.65),
    panels = list(
      list(id = "A", role = "panel", source_width_pt = 19.0 * pt_per_cm,
           source_height_pt = 16.57 * pt_per_cm,
           placed_box_pt = c(0, 7.07 * pt_per_cm, 19.0 * pt_per_cm, 23.65 * pt_per_cm),
           ink_box_pt = c(0.20 * pt_per_cm, 7.20 * pt_per_cm, 18.80 * pt_per_cm, 23.49 * pt_per_cm)),
      list(id = "B", role = "panel", source_width_pt = 19.0 * pt_per_cm,
           source_height_pt = 6.39 * pt_per_cm,
           placed_box_pt = c(0, 0, 19.0 * pt_per_cm, 6.39 * pt_per_cm),
           ink_box_pt = c(0.20 * pt_per_cm, 0.12 * pt_per_cm, 18.80 * pt_per_cm, 6.27 * pt_per_cm))
    ),
    gaps = list(list(panel_a = "A", panel_b = "B", axis = "y", kind = "vertical", min_pt = 14.6, max_pt = 29.1))
  ),
  file.path(OUT_ROOT, "geometry_ledger_Figure_6.json"), pretty = TRUE, auto_unbox = TRUE, digits = 10
)

message("Geoderma V03 Figure 6 saved -> ", standalone_pdf)
quit(save = "no", status = 0)

if (file.exists(combined_pdf)) unlink(combined_pdf)
merge_status <- system2(
  "pdfunite",
  c(shQuote(CURRENT_S3), shQuote(standalone_pdf), shQuote(combined_pdf))
)
if (!identical(merge_status, 0L)) stop("Could not append amoA validation page to Figure S3")

audit <- c(
  "V20.8 candidate Figure S3 amoA validation audit",
  paste0("tree_tips=", ape::Ntip(tree_rooted)),
  paste0("best_model=", best_model),
  paste0("local_K10944_queries=", length(query_ids)),
  paste0("archaeal_AmoA_references=", sum(metadata$reference_class == "Archaeal AmoA")),
  paste0("bacterial_AmoA_references=", sum(metadata$reference_class == "Bacterial AmoA")),
  paste0("PmoA_references=", sum(metadata$reference_class == "PmoA")),
  paste0("archaeal_clade_SHaLRT_UFBoot=", archaeal_support),
  "MAG005_targeted_interpretation=Archaeal AmoA supported",
  "MAG070_targeted_interpretation=Archaeal AmoA supported",
  "local_bacterial_K10944_interpretation=unresolved divergent homologs",
  "activity_inference=none",
  paste0("standalone_pdf=", standalone_pdf),
  paste0("combined_pdf=", combined_pdf)
)
write_lines(audit, file.path(SOURCE_DIR, "FigureS3_amoA_candidate_audit.txt"))

message("Standalone validation page -> ", standalone_pdf)
message("Three-page Figure S3 candidate -> ", combined_pdf)
