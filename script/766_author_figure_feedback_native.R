#!/usr/bin/env Rscript
# Native display-only changes requested after the author's first comparison.
invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent=TRUE)))
PROJECT <- Sys.getenv("CHILE_PROJECT_ROOT", unset=getwd())
AUDIT <- file.path(PROJECT,"analysis/figure_table_refinement_2026-10-18_round2")
NATIVE <- file.path(AUDIT,"native")
PARENT <- file.path(PROJECT,"Environmental_Microbiome/Coauthor_Review_2026-10-01_XW")
dir.create(NATIVE,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(dplyr);library(readr);library(ggplot2);library(grid);library(jsonlite);library(ggtext)})
source(file.path(PROJECT,"script/_theme_sbb.R"))
FONT <- "Arial"
pdfFonts(Arial=pdfFonts("ArialMT")[[1]])
if (capabilities("aqua")) quartzFonts(Arial=quartzFont(c("Arial","Arial Bold","Arial Italic","Arial Bold Italic")))

eval_region <- function(path,start,end,env=globalenv()) {
  lines <- readLines(path,warn=FALSE)
  first <- which(startsWith(lines,start));last <- which(startsWith(lines,end))
  stopifnot(length(first)==1,length(last)==1,first<last)
  expr <- parse(path,keep.source=TRUE);refs <- attr(expr,"srcref")
  for(i in seq_along(expr)) {
    if(refs[[i]][1]<first || refs[[i]][3]>=last) next
    x <- expr[[i]]
    if(is.call(x) && identical(x[[1]],as.name("write_tsv"))) next
    eval(x,env)
  }
}

figure5 <- function() {
  # Reproduce the V19 accepted B grammar, not the rejected horizontal headings.
  src <- readLines(file.path(PROJECT,"script/553_v19_figure5_6.R"))
  first <- which(startsWith(src,"d <- read_tsv"))
  last <- which(startsWith(src,"ggsave(file.path(O,"))[1]-1L
  code <- paste(src[first:last],collapse="\n")
  code <- sub("Geoderma/V18/Supplementary_Data/MAG_phylum_marker_carriers.tsv",
    "Environmental_Microbiome/Coauthor_Review_2026-10-01_XW/Supplementary_Data/MAG_phylum_marker_carriers.tsv",code,fixed=TRUE)
  eval(parse(text=code))
  old <- unserialize(serialize(p_carriers_readable,NULL))
  # A single multiplier keeps zero at zero and bubble area proportional to prevalence.
  p_carriers_readable <- p_carriers_readable +
    scale_size_area(name="Carrier prevalence\nwithin phylum (%)",max_size=6.0,
      limits=c(0,100),breaks=c(25,50,75,100),
      guide=guide_legend(order=1,nrow=1,byrow=TRUE)) +
    theme(strip.background=element_rect(fill="#F8F9FA",colour="#D9DEE2",linewidth=.35))
  p_carriers_readable$layers[[1]]$aes_params$fill <- "#F8F9FA"
  p_carriers_readable$layers[[1]]$aes_params$colour <- "#D9DEE2"
  p_carriers_readable$layers[[1]]$aes_params$linewidth <- .35
  for(id in c("baseline","revised")) {
    p <- if(id=="baseline") old else p_carriers_readable
    ggsave(file.path(NATIVE,paste0("Figure_5_B_",id,".pdf")),p,
      width=487/72,height=444/72,units="in",device=cairo_pdf,bg="white")
  }
  built_old <- ggplot_build(old);built_new <- ggplot_build(p_carriers_readable)
  old_bubbles <- built_old$data[[4]];new_bubbles <- built_new$data[[4]]
  stopifnot(identical(old_bubbles$x,new_bubbles$x),identical(old_bubbles$y,new_bubbles$y),
    identical(old_bubbles$fill,new_bubbles$fill),max(abs(new_bubbles$size/old_bubbles$size-1.25))<1e-12,
    identical(old$data,p_carriers_readable$data))
  write_tsv(old_bubbles,file.path(AUDIT,"Figure5_bubbles_baseline.tsv"))
  write_tsv(new_bubbles,file.path(AUDIT,"Figure5_bubbles_revised.tsv"))
  write_tsv(d,file.path(AUDIT,"Figure5_source_values.tsv"),quote="needed")
  write_json(list(rows=nrow(d),observations_unchanged=TRUE,positive_bubbles=nrow(new_bubbles),
    max_size_old_mm=4.8,max_size_new_mm=6,diameter_multiplier=1.25,area_multiplier=1.5625,
    zero_size_offset=0,category_fill="#F8F9FA",original_legend_structure=TRUE,
    native_size_pt=c(487,444),new_analysis=FALSE),file.path(AUDIT,"Figure5_QA.json"),pretty=TRUE,auto_unbox=TRUE)
}

figureS2 <- function() {
  suppressPackageStartupMessages({library(ape);library(circlize);library(ComplexHeatmap);library(cowplot);library(ggnewscale);library(ggplotify);library(ggtree);library(tidyr);library(stringr);library(scales)})
  e <- new.env(parent=globalenv());te <- new.env(parent=e)
  te$MAG_DIR <- file.path(PROJECT,"data/processed/metagenome_output/MAGs")
  te$META_FILE <- file.path(PROJECT,"data/metadata/env_2024_meta.tsv")
  datafile <- file.path(PROJECT,"script/212_build_v18_s3_clean.R")
  plotfile <- file.path(PROJECT,"script/471_render_geoderma_v03_figure_s2.R")
  eval_region(datafile,"# Shared within-figure evidence palettes.","# ---- 4.",te)
  eval_region(datafile,"# ---- 5.","p_tree_base <-",te)
  eval_region(datafile,"# ---- 8.","tree_theme <-",te)
  e$tree_env <- te
  eval_region(plotfile,"focal_mags <-","# ---- Dedicated legends",e)
  base_b <- unserialize(serialize(e$p_b,NULL))
  ids <- which(vapply(e$p_b$layers,function(l) is.data.frame(l$data)&&all(c("lab","value","x","y","width")%in%names(l$data)),logical(1)))
  stopifnot(length(ids)==4,sum(e$p_a$data$isTip)==221,sum(e$p_b$data$isTip)==7)
  cells <- list();bounds <- list()
  for(i in seq_along(ids)) {
    l <- e$p_b$layers[[ids[i]]];cc <- l$data
    cells[[i]] <- cc %>% mutate(value=as.character(value),track=c("Class","Completeness","Abundance","K10944")[i])
    bounds[[i]] <- data.frame(xmin=min(cc$x-cc$width/2),xmax=max(cc$x+cc$width/2),ymin=min(cc$y)-.5,ymax=max(cc$y)+.5)
    l$aes_params$colour <- "#CDD3D7";l$aes_params$linewidth <- .4/ggplot2::.pt
    e$p_b$layers[[ids[i]]] <- l
  }
  e$p_b <- e$p_b+geom_rect(data=bind_rows(bounds),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
    inherit.aes=FALSE,fill=NA,colour="#959FA6",linewidth=.65/ggplot2::.pt)
  header_ids <- which(vapply(e$p_b$layers,function(l) is.data.frame(l$data)&&
    "from"%in%names(l$data)&&any(grepl("K10944",as.character(l$data$from))),logical(1)))
  stopifnot(length(header_ids)==1)
  e$p_b$layers[[header_ids]]$data$y <- e$p_b$layers[[header_ids]]$data$y+.20
  for(id in c("baseline","revised")) {
    ggsave(file.path(NATIVE,paste0("Figure_S2_B_",id,".pdf")),if(id=="baseline") base_b else e$p_b,
      width=19,height=14.2*.49,units="cm",device=cairo_pdf,bg="white")
  }
  revised_b <- e$p_b
  all_lines <- readLines(plotfile,warn=FALSE)
  lo <- which(startsWith(all_lines,"# ---- Dedicated legends"))
  hi <- which(startsWith(all_lines,"# Lock each source panel"))-1L
  legend_code <- paste(all_lines[lo:hi],collapse="\n")
  e$page_width_cm <- 18.3;e$page_b_width_cm <- 19;e$page_height_cm <- 19.8
  e$top_height_cm <- 14.2;e$top_height_b_cm <- 14.2*.49
  e$legend_height_cm <- 5.4;e$top_bottom_cm <- 5.6
  e$page_b_height_cm <- e$top_height_b_cm+5.6
  for(id in c("baseline","revised")) {
    code <- legend_code
    if(id=="revised") {
      code <- sub("length.out = length(pal_seq)),\n  pal_seq",
        'length.out = 3),\n  c("#F7FBFF", "#9ECAE1", "#4A6D8C")',code,fixed=TRUE)
      code <- sub("at = pretty(tree_env$abundance_limits, n = 3),",
        'at = c(tree_env$abundance_limits[1], 0.5, tree_env$abundance_limits[2]),\n  labels = c(sprintf("%.2f", tree_env$abundance_limits[1]), "0.5", sprintf("%.2f", tree_env$abundance_limits[2])),',code,fixed=TRUE)
      code <- sub("colorRamp2(c(-3, -1.5, 0, 1.5), colorRampPalette(pal_seq)(4))",
        "colorRamp2(seq(-3, 1.5, length.out=length(pal_seq)), pal_seq)",code,fixed=TRUE)
    }
    eval(parse(text=code),e)
    e$p_b <- if(id=="baseline") base_b else revised_b
    eval_region(plotfile,"page_a <-","output <-",e)
    ggsave(file.path(NATIVE,paste0("Figure_S2_A_page_",id,".pdf")),e$page_a,
      width=18.3,height=19.8,units="cm",device=cairo_pdf,bg="white")
    ggsave(file.path(NATIVE,paste0("Figure_S2_B_page_",id,".pdf")),e$page_b,
      width=19,height=e$page_b_height_cm,units="cm",device=cairo_pdf,bg="white")
  }
  # Export native interpolation samples for a full-size colour key. Do not alter
  # the original large legend layout or the bacterial tree to repair a key.
  ring_limits <- te$abundance_limits
  key <- function(name,colours,values=NULL,limits) {
    t <- seq(0,1,length.out=512)
    data.frame(scale=name,t=t,value=limits[1]+t*diff(limits),
      fill=scales::gradient_n_pal(colours,values=values)(t))
  }
  keys <- bind_rows(
    key("A_abundance",c("#F7FBFF","#9ECAE1","#4A6D8C"),limits=ring_limits),
    key("B_abundance",pal_seq,limits=c(-3,1.5)),
    key("B_completeness",te$completeness_colours,
      scales::rescale(c(50,65,80,90,100),from=c(50,100)),c(50,100)))
  write_tsv(keys,file.path(AUDIT,"S2_native_colour_keys.tsv"))
  write_tsv(bind_rows(cells),file.path(AUDIT,"Figure_S2_B_plotted_cells.tsv"),quote="needed")
  write_tsv(te$mag_metrics %>% filter(Domain=="Bacteria") %>% transmute(MAG,phylum,focal=MAG%in%e$focal_mags,
    mean_relative_abundance=mean_rel_abundance,mean_per_mille,log10_mean_per_mille),
    file.path(AUDIT,"Figure_S2A_bacterial_tree_abundance.tsv"))
  write_tsv(te$ar_meta %>% select(MAG,phylum,class,Completeness,pmo_present),file.path(AUDIT,"Figure_S2B_archaeal_metadata.tsv"))
  write_tsv(as.data.frame(te$ar_abundance_matrix) %>% tibble::rownames_to_column("MAG"),file.path(AUDIT,"Figure_S2B_abundance_matrix.tsv"))
  stopifnot(identical(base_b$data,e$p_b$data))
  write_json(list(bacteria=221,archaea=7,tree_coordinates_unchanged=TRUE,
    A_abundance_limits=ring_limits,B_native_size_cm=c(19,14.2*.49),
    original_legend_geometry_preserved=TRUE,K10944_header_offset_y=.20,
    new_analysis=FALSE),file.path(AUDIT,"FigureS2_QA.json"),pretty=TRUE,auto_unbox=TRUE,digits=17)
}

mode <- commandArgs(trailingOnly=TRUE)
stopifnot(length(mode)==1)
switch(mode,figure5=figure5(),figureS2=figureS2(),stop("Unknown mode"))
