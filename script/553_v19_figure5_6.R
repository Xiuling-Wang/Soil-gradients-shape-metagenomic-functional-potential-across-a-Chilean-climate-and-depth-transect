# V19 geometry-only derivatives of the established renderers; no model refit.
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(ggtext)
  library(patchwork); library(jsonlite)
})
PROJECT <- Sys.getenv("CHILE_PROJECT_ROOT", unset=getwd())
Q <- file.path(PROJECT, "analysis/Geoderma_V19_figure_table_revision")
O <- file.path(Q, "figures")
dir.create(O, recursive = TRUE, showWarnings = FALSE)
source(file.path(PROJECT, "script/_theme_sbb.R"))
FONT <- "Arial"
d <- read_tsv(file.path(PROJECT, "Geoderma/V18/Supplementary_Data/MAG_phylum_marker_carriers.tsv"), show_col_types=FALSE)
phylum_denominators <- d %>% distinct(phylum, n_MAG) %>% arrange(desc(n_MAG), phylum)
carrier_phyla <- phylum_denominators$phylum
module_order <- c("Nitrogen transformations", "Pedogenic resource acquisition", "Alternative-energy chemolithotrophy")
category_order <- c("N fixation", "Hydroxylamine oxidation", "Nitrate/nitrite redox", "Denitrification",
  "Nitrite-to-ammonium reduction", "Assimilatory nitrate reduction", "Urea hydrolysis",
  "Organic-acid weathering", "P acquisition", "Catecholate siderophore synthesis", "EPS and aggregation",
  "Sulfur oxidation", "Sulfur redox, direction unresolved", "Hydrogenase homologs", "Putative CO oxidation", "CBB carbon-fixation markers")
d$module_label <- sub("Mineral-resource acquisition", "Pedogenic resource acquisition", d$module_label)
stopifnot(setequal(d$category_label, category_order), all(d$n_MAG >= 5), nrow(d)==112,
  max(abs(d$carrier_pct-100*d$carrier_n/d$n_MAG)) < 1e-8)
d$phylum <- factor(d$phylum, levels=carrier_phyla)
d$module_label <- factor(d$module_label, levels=module_order)
d$category_label <- factor(d$category_label, levels=rev(category_order))
p_carriers <- list(data=d)
phylum_axis_labels <- setNames(paste0(phylum_denominators$phylum, "<br>(n = ", phylum_denominators$n_MAG, ")"), carrier_phyla)
marker_fill_max <- max(d$mean_marker_hits_per_1000_ORFs)
src <- readLines(file.path(PROJECT, "script/467_render_geoderma_v03_figure5.R"))
first <- which(grepl("^recode_carrier_module <-", src))
last <- which(grepl("^# Geoderma V03: two-panel page", src))-1L
eval(parse(text=src[first:last]))
# A full-size bubble now fits both the row pitch and the facet-edge padding.
p_carriers_readable <- p_carriers_readable +
  scale_y_discrete(expand=expansion(add=0.68)) +
  scale_size_area(name="Carrier prevalence\nwithin phylum (%)", max_size=4.8,
    limits=c(0,100), breaks=c(25,50,75,100),
    guide=guide_legend(order=1, nrow=1, byrow=TRUE))
p_carriers_readable$layers[[3]]$data$frame_xmin <- 0.77
ggsave(file.path(O,"Figure_5_B.pdf"), p_carriers_readable,
  width=487/72, height=(674-230)/72, units="in", bg="white", device=cairo_pdf)
write_json(list(max_size_mm=4.8, stroke=0.28, facet_padding_rows=0.68,
  source_rows=nrow(d), carrier_phyla=carrier_phyla, frame_xmin=0.77,
  no_source_values_changed=TRUE, colour_max=marker_fill_max),
  file.path(O,"Figure_5_render_parameters.json"), pretty=TRUE, auto_unbox=TRUE)

# Reuse the full established Figure 6 renderer with bound V18 tree/context inputs.
code <- paste(readLines(file.path(PROJECT,"script/470_render_geoderma_v03_figure6.R")), collapse="\n")
replace_once <- function(text, old, new) {
  stopifnot(lengths(regmatches(text, gregexpr(old,text,fixed=TRUE)))==1L)
  sub(old,new,text,fixed=TRUE)
}
source_dir <- file.path(PROJECT,"Geoderma/V18/Supplementary_Data/AmoA_sequences")
code <- replace_once(code, 'PHYLO_DIR <- file.path(INPUT_ROOT, "phylogeny")', paste0('PHYLO_DIR <- "',source_dir,'"'))
code <- replace_once(code, 'SOURCE_INPUT_DIR <- file.path(INPUT_ROOT, "source_data")', paste0('SOURCE_INPUT_DIR <- "',source_dir,'"'))
code <- replace_once(code, '  file.exists(CURRENT_S3), file.exists(TREE_FILE), file.exists(IQTREE_FILE),', '  file.exists(TREE_FILE), file.exists(IQTREE_FILE),')
binding <- '
focal_node <- p_tree$data[p_tree$data$node == archaeal_mrca, ]
parent_node <- p_tree$data[p_tree$data$node == focal_node$parent, ]
stopifnot(nrow(focal_node)==1, nrow(parent_node)==1, archaeal_support=="84.4/74")
stem_x <- (focal_node$x + parent_node$x)/2
p_tree <- p_tree + annotate("segment", x=stem_x, xend=stem_x,
  y=focal_node$y, yend=focal_node$y+0.79, colour="grey25", linewidth=0.35)
jsonlite::write_json(list(node=archaeal_mrca, parent=focal_node$parent,
  support=archaeal_support, tips=archaeal_tips, node_x=focal_node$x,
  parent_x=parent_node$x, branch_y=focal_node$y, leader_x=stem_x,
  leader_end_y=focal_node$y+0.79, tree_file=TREE_FILE,
  inference="putative archaeal AmoA; no topology or branch-length change"),
  file.path(OUT_ROOT,"support_binding.json"), pretty=TRUE, auto_unbox=TRUE)
'
code <- replace_once(code, '# ---- MAG005/MAG070 context caveat', paste0(binding, '\n# ---- MAG005/MAG070 context caveat'))
Sys.setenv(CHILE_READABILITY_FIGURE5_ROOT=file.path(O,"figure6"),
  CHILE_NATIVE_PANEL_DIR=file.path(O,"figure6/native_panels"))
eval(parse(text=code))
