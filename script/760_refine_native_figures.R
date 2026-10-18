#!/usr/bin/env Rscript
# Display-only derivatives. Existing observations, scales, topology and tests stay fixed.
invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent=TRUE)))
PROJECT <- Sys.getenv("CHILE_PROJECT_ROOT", unset=getwd())
AUDIT <- file.path(PROJECT, "analysis/figure_table_refinement_2026-10-18")
DELIVERY <- file.path(PROJECT, "Environmental_Microbiome/Figure_table_review_2026-10-18")
PARENT <- file.path(PROJECT, "Environmental_Microbiome/Coauthor_Review_2026-10-01_XW")
NATIVE <- file.path(AUDIT, "native")
dir.create(NATIVE, recursive=TRUE, showWarnings=FALSE)
source(file.path(PROJECT, "script/_theme_sbb.R"))
FONT <- "Arial"
suppressPackageStartupMessages({library(ggplot2); library(grid); library(dplyr); library(readr); library(jsonlite)})
pdfFonts(Arial=pdfFonts("ArialMT")[[1]])
if (capabilities("aqua")) quartzFonts(Arial=quartzFont(c("Arial", "Arial Bold", "Arial Italic", "Arial Bold Italic")))
# Historical renderers register the same font repeatedly; retain the first registration.
pdfFonts <- function(...) {
  args <- list(...)
  if(length(args) && !is.null(names(args)) && "Arial" %in% names(args)) return(invisible(NULL))
  grDevices::pdfFonts(...)
}
quartzFonts <- function(...) {
  args <- list(...)
  if(length(args) && !is.null(names(args)) && "Arial" %in% names(args)) return(invisible(NULL))
  grDevices::quartzFonts(...)
}

replace_exact <- function(code, old, new, count=1L) {
  found <- lengths(regmatches(code, gregexpr(old, code, fixed=TRUE)))
  stopifnot(found == count)
  gsub(old, new, code, fixed=TRUE)
}

eval_region <- function(path, start, end, env=globalenv()) {
  lines <- readLines(path, warn=FALSE)
  first <- which(startsWith(lines,start)); last <- which(startsWith(lines,end))
  stopifnot(length(first)==1, length(last)==1, first<last)
  expr <- parse(path, keep.source=TRUE); refs <- attr(expr,"srcref")
  for (i in seq_along(expr)) {
    if (refs[[i]][1] < first || refs[[i]][3] >= last) next
    x <- expr[[i]]
    if (is.call(x) && identical(x[[1]],as.name("write_tsv"))) next
    eval(x,env)
  }
}

heatmaps <- function() {
  src <- paste(readLines(file.path(PROJECT,"script/465_render_geoderma_v03_figure3_split.R")),collapse="\n")
  src <- replace_exact(src,'HEATMAP_BODY_WIDTH_CM <- 11.4','HEATMAP_BODY_WIDTH_CM <- 11.55')
  src <- replace_exact(src,'width = unit(40, "mm")','width = unit(38.5, "mm")')
  src <- replace_exact(src,'page_width_cm <- 18.3 + (HEATMAP_BODY_WIDTH_CM - 10.8)','page_width_cm <- 18.9')
  src <- replace_exact(src,'rect_gp = gpar(col = "grey70", lwd = 0.35),\n    border = TRUE,',
    'rect_gp = gpar(col = "#D4D8DB", lwd = 0.25),\n    border = "#AAB2B7",')
  src <- gsub('fill = "#F1F3F5", col = "#CCD1D5", lwd = 0.55','fill = "#F8F9FA", col = "#D9DEE2", lwd = 0.35',src,fixed=TRUE)
  src <- gsub('fill = "#F1F3F5", col = "#D8DDE1", lwd = 0.35','fill = "#FAFBFC", col = "#E4E7E9", lwd = 0.25',src,fixed=TRUE)
  src <- replace_exact(src,'group_label <- as.character(category_split[index[[1]]])',
    'group_label <- as.character(category_split[index[[1]]])\n      group_label <- dplyr::recode(group_label, DNRA="Nitrite-to-ammonium\\nreduction", `Organic-acid weathering`="Organic-acid\\nweathering", `Hydrogenase homologues`="Hydrogenase\\nhomologues")')
  for (mode in c("figure3","figure4")) {
    root <- file.path(AUDIT,mode)
    Sys.setenv(CHILE_READABILITY_FIGURE3_ROOT=root,CHILE_GEODERMA_SPLIT_MODE=mode)
    env <- new.env(parent=globalenv())
    eval(parse(text=src),env)
    write_json(env$sample_values,file.path(root,"source_data/plotted_values.json"),dataframe="rows",pretty=TRUE,digits=17)
    write_tsv(env$sample_values,file.path(root,"source_data/Figure_3_sample_resolved_marker_values.tsv"),quote="needed")
    num <- if(mode=="figure3") 3 else 4
    file.copy(file.path(root,"figures",paste0("Figure_",num,".pdf")),DELIVERY,overwrite=TRUE)
  }
}

figure5 <- function() {
  # Same 112 frozen phylum/category observations; only physical layout changes.
  d <- read_tsv(file.path(PARENT,"Supplementary_Data/MAG_phylum_marker_carriers.tsv"),show_col_types=FALSE)
  phyla <- d %>% distinct(phylum,n_MAG) %>% arrange(desc(n_MAG),phylum)
  cats <- c("N fixation","Hydroxylamine oxidation","Nitrate/nitrite redox","Denitrification",
    "Nitrite-to-ammonium reduction","Assimilatory nitrate reduction","Urea hydrolysis",
    "Organic-acid weathering","P acquisition","Catecholate siderophore synthesis","EPS and aggregation",
    "Sulfur oxidation","Sulfur redox, direction unresolved","Hydrogenase homologs","Putative CO oxidation","CBB carbon-fixation markers")
  stopifnot(nrow(d)==112L,setequal(d$category_label,cats),nrow(phyla)==7L,
    max(abs(d$carrier_pct-100*d$carrier_n/d$n_MAG))<1e-8)
  labels <- c("N fixation","Hydroxylamine oxidation","Nitrate/nitrite redox","Denitrification",
    "Nitrite-to-ammonium reduction","Assimilatory nitrate reduction","Urea hydrolysis",
    "Organic-acid weathering","P acquisition","Catecholate siderophore synthesis","EPS and aggregation",
    "Sulfur oxidation","Sulfur redox (direction unresolved)","Hydrogenase homologues","Putative CO oxidation","CBB-associated carbon fixation")
  ys <- c(seq(17,44.6,length.out=7),seq(57.2,71,length.out=4),seq(83.6,102,length.out=5))
  xs <- seq(66,174,length.out=7)
  d$x_mm <- xs[match(d$phylum,phyla$phylum)]
  d$top_mm <- ys[match(d$category_label,cats)]
  fillmax <- max(d$mean_marker_hits_per_1000_ORFs)
  colourfun <- scales::gradient_n_pal(pal_seq)
  d$fill <- colourfun(d$mean_marker_hits_per_1000_ORFs/fillmax)
  d$diameter_mm <- 4.6*sqrt(d$carrier_pct/100)
  w <- 180; h <- 150
  grDevices::cairo_pdf(file.path(NATIVE,"Figure_5_B.pdf"),width=w/25.4,height=h/25.4,family=FONT)
  grid.newpage(); grid.rect(gp=gpar(fill="white",col=NA))
  tx <- function(x,y,label,size=9,bold=FALSE,just="left",rot=0,col="#222222") {
    grid.text(label,x=unit(x,"mm"),y=unit(h-y,"mm"),just=just,rot=rot,
      gp=gpar(fontfamily=FONT,fontsize=size,fontface=if(bold)"bold" else "plain",col=col,lineheight=1.0))
  }
  tx(1,4,"B",15.5,TRUE)
  for (m in seq_len(3)) {
    top <- c(11,51.2,77.6)[m]
    tx(2,top,c("Nitrogen transformations","Mineral-resource acquisition","Alternative-energy chemolithotrophy")[m],9.6,TRUE)
    grid.lines(x=unit(c(2,178),"mm"),y=unit(rep(h-top-2.6,2),"mm"),gp=gpar(col="#C4CBD0",lwd=.45))
    idx <- list(1:7,8:11,12:16)[[m]]
    for (x in xs) grid.lines(x=unit(c(x,x),"mm"),y=unit(h-c(min(ys[idx])-2.3,max(ys[idx])+2.3),"mm"),gp=gpar(col="#E3E7E9",lwd=.35))
  }
  for (i in seq_along(cats)) {
    tx(59.5,ys[i],labels[i],8.6,just="right")
    grid.lines(x=unit(c(63,177),"mm"),y=unit(rep(h-ys[i],2),"mm"),gp=gpar(col="#EBEEF0",lwd=.35))
  }
  for (i in which(d$carrier_n>0)) grid.circle(x=unit(d$x_mm[i],"mm"),y=unit(h-d$top_mm[i],"mm"),
    r=unit(d$diameter_mm[i]/2,"mm"),gp=gpar(fill=d$fill[i],col="#48545D",lwd=.55))
  for(i in seq_len(nrow(phyla))) {
    tx(xs[i]+.6,108,paste0(phyla$phylum[i],"\n(n = ",phyla$n_MAG[i],")"),8.5,just="right",rot=55)
  }
  tx(2,135,"Carrier prevalence within phylum (%)",9,TRUE)
  for(i in seq_along(c(25,50,75,100))) {
    val <- c(25,50,75,100)[i]; x <- c(8,29,50,71)[i]
    grid.circle(x=unit(x,"mm"),y=unit(h-143,"mm"),r=unit(2.3*sqrt(val/100),"mm"),gp=gpar(fill="white",col="#48545D",lwd=.55))
    tx(x+4,143,as.character(val),8.5)
  }
  tx(101,133.8,"Mean accepted marker hits",8.8,TRUE)
  tx(101,138,"per 1,000 ORFs",8.8,TRUE)
  for(i in 1:160) grid.rect(x=unit(101+(i-.5)*.3,"mm"),y=unit(h-142,"mm"),width=unit(.31,"mm"),height=unit(3,"mm"),gp=gpar(fill=colourfun((i-.5)/160),col=NA))
  for(v in c(0,.5,1,1.5)) if(v<=fillmax+.01) tx(101+48*v/fillmax,147,sprintf("%.1f",v),8.2,just="centre")
  dev.off()
  write_tsv(d,file.path(AUDIT,"Figure_5_B_plotted_values.tsv"))
  write_json(list(native_mm=c(w,h),source_rows=nrow(d),displayed_positive=sum(d$carrier_n>0),
    excluded_rows=0,area_scale_max_mm=4.6,fill_limits=c(0,fillmax),colours=pal_seq,
    source="MAG_phylum_marker_carriers.tsv",new_analysis=FALSE),file.path(AUDIT,"Figure_5_B_QA.json"),pretty=TRUE,auto_unbox=TRUE)
}

figureS1 <- function() {
  e <- new.env(parent=globalenv())
  Sys.setenv(CHILE_READABILITY_FIGURE2_ROOT=file.path(AUDIT,"S1"))
  expr <- parse(file.path(PROJECT,"script/446_render_v26_figure2_s1.R"))
  for(x in expr) {
    if(is.call(x) && identical(x[[1]],as.name("<-")) && identical(x[[2]],as.name("surface_pc_vectors"))) break
    eval(x,e)
  }
  p <- e$make_pcoa_panel(e$depth_pcoa_data,NULL,.389,"0.361","C",.20,FALSE)
  labels <- e$depth_pcoa_data %>% mutate(
    display_pit=sub("^NB","NA",site_pit),
    label_x=c(-.120,-.069,.020,.098,-.030,.119,.043,.035,-.080),
    label_y=c(.037,.026,.058,.008,.033,-.047,.015,-.037,-.140))
  dx <- labels$label_x-labels$axis1;dy <- labels$label_y-labels$axis2
  distance <- sqrt(dx^2+dy^2)
  labels$start_x <- labels$axis1+.008*dx/distance
  labels$start_y <- labels$axis2+.008*dy/distance
  labels$end_x <- labels$label_x-.011*dx/distance
  labels$end_y <- labels$label_y-.011*dy/distance
  p <- p + geom_segment(data=labels,aes(x=start_x,y=start_y,xend=end_x,yend=end_y),
      inherit.aes=FALSE,colour="#6D7880",linewidth=.24) +
    geom_text(data=labels,aes(x=label_x,y=label_y,label=display_pit),inherit.aes=FALSE,
      family=FONT,fontface="bold",size=3.05,colour="#262B30") +
    theme(plot.margin=margin(6,15,6,7))
  ggsave(file.path(NATIVE,"Figure_S1_C.pdf"),p,width=9.27,height=8.765,units="cm",device=cairo_pdf,bg="white")
  write_tsv(e$depth_pcoa_data,file.path(AUDIT,"Figure_S1_C_coordinates.tsv"))
  write_tsv(labels,file.path(AUDIT,"Figure_S1_C_label_positions.tsv"))
}

figure6 <- function() {
  src <- paste(readLines(file.path(PROJECT,"script/695_v35_repair_figure6.R")),collapse="\n")
  src <- replace_exact(src,"OUT <- file.path(PROJECT,'analysis/EM_V35_coauthor_candidate_2026-09-30')",paste0("OUT <- '",AUDIT,"'"))
  src <- replace_exact(src,"Q <- file.path(OUT,'INTERNAL_QA/figure6_build')","Q <- file.path(OUT,'figure6')")
  style <- paste0(
    'once(\'tree_rooted, layout = "rectangular", size = 0.52, colour = "grey32"\',\'tree_rooted, layout = "rectangular", size = 0.55, colour = "grey22"\')\n',
    'once(\'node = archaeal_mrca, fill = "#E5F3EE", alpha = 0.68,\',\'node = archaeal_mrca, fill = "#E5F3EE", alpha = 0.25,\')\n',
    'code <- gsub("color:#6B7280", "color:#59646E;font-size:8.5pt", code, fixed=TRUE)\n')
  src <- replace_exact(src,"# No legacy append-to-supplement or process quit; only the established standalone export.",paste0(style,"\n# Only standalone exports."))
  eval(parse(text=src),new.env(parent=globalenv()))
  file.copy(file.path(AUDIT,"figure6/figures/Figure_6.pdf"),DELIVERY,overwrite=TRUE)
}

figureS2 <- function() {
  suppressPackageStartupMessages({library(ape);library(circlize);library(ComplexHeatmap);library(cowplot);library(ggnewscale);library(ggplotify);library(ggtext);library(ggtree);library(tidyr);library(stringr);library(scales)})
  plotfile <- file.path(PROJECT,"script/471_render_geoderma_v03_figure_s2.R")
  datafile <- file.path(PROJECT,"script/212_build_v18_s3_clean.R")
  e <- new.env(parent=globalenv())
  te <- new.env(parent=e)
  te$MAG_DIR <- file.path(PROJECT,"data/processed/metagenome_output/MAGs")
  te$META_FILE <- file.path(PROJECT,"data/metadata/env_2024_meta.tsv")
  eval_region(datafile,"# Shared within-figure evidence palettes.","# ---- 4.",te)
  eval_region(datafile,"# ---- 5.","p_tree_base <-",te)
  eval_region(datafile,"# ---- 8.","tree_theme <-",te)
  e$tree_env <- te
  eval_region(plotfile,"focal_mags <-","# ---- Dedicated legends",e)
  original_a_data <- e$p_a$data
  original_b_data <- e$p_b$data
  stopifnot(sum(original_a_data$isTip)==221L,sum(original_b_data$isTip)==7L)
  for(i in seq_along(e$p_a$layers)) {
    if(any(grepl("Segment",class(e$p_a$layers[[i]]$geom)))) {
      e$p_a$layers[[i]]$aes_params$colour <- "grey38"
      e$p_a$layers[[i]]$aes_params$linewidth <- .24
    }
  }
  # Thin borders belong to each native data track, never to a guessed overlay.
  ids <- which(vapply(e$p_b$layers,function(l) is.data.frame(l$data)&&all(c("lab","value","x","y","width")%in%names(l$data)),logical(1)))
  stopifnot(length(ids)==4)
  cells <- list(); bounds <- list()
  for(i in seq_along(ids)) {
    l <- e$p_b$layers[[ids[i]]]; c <- l$data
    cells[[i]] <- c %>% mutate(value=as.character(value),track=c("Class","Completeness","Abundance","K10944")[i])
    bounds[[i]] <- data.frame(xmin=min(c$x-c$width/2),xmax=max(c$x+c$width/2),ymin=min(c$y)-.5,ymax=max(c$y)+.5)
    l$aes_params$colour <- "#CDD3D7"; l$aes_params$linewidth <- .4/ggplot2::.pt
    e$p_b$layers[[ids[i]]] <- l
  }
  e$p_b <- e$p_b + ylim(.25,8.65) + geom_rect(data=bind_rows(bounds),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),inherit.aes=FALSE,fill=NA,colour="#959FA6",linewidth=.65/ggplot2::.pt)
  for(l in e$p_b$layers) if(is.data.frame(l$data)&&"from"%in%names(l$data)&&any(grepl("K10944",as.character(l$data$from)))) l$data$y <- l$data$y+.18
  ggsave(file.path(NATIVE,"Figure_S2_A.pdf"),e$p_a,width=18.3,height=14.2,units="cm",device=cairo_pdf,bg="white")
  ggsave(file.path(NATIVE,"Figure_S2_B.pdf"),e$p_b,width=19,height=6.8,units="cm",device=cairo_pdf,bg="white")
  # Compact explicit keys retain every category and all continuous scale limits.
  legend_pdf <- function(path,w,h,body) {
    cairo_pdf(path,width=w/25.4,height=h/25.4,family=FONT);grid.newpage();grid.rect(gp=gpar(fill="white",col=NA));body(w,h);dev.off()
  }
  textmm <- function(x,y,s,size=9.4,bold=FALSE) grid.text(s,x=unit(x,"mm"),y=unit(y,"mm"),just="left",gp=gpar(fontfamily=FONT,fontsize=size,fontface=if(bold)"bold" else "plain"))
  colourbar <- function(x,y,w,colours,limits,ticks,title,values=NULL) {
    textmm(x,y+9,title,9.5,TRUE)
    fn <- scales::gradient_n_pal(colours,values=values)
    for(k in 1:140) grid.rect(x=unit(x+(k-.5)*w/140,"mm"),y=unit(y+4,"mm"),width=unit(w/140+.01,"mm"),height=unit(2.8,"mm"),gp=gpar(fill=fn((k-.5)/140),col=NA))
    for(v in ticks[ticks>=limits[1]&ticks<=limits[2]]) grid.text(format(v,trim=TRUE),x=unit(x+w*(v-limits[1])/diff(limits),"mm"),y=unit(y,"mm"),gp=gpar(fontfamily=FONT,fontsize=8.5))
  }
  legend_pdf(file.path(NATIVE,"Figure_S2_A_legend.pdf"),183,39,function(w,h) {
    colourbar(4,24,45,c("#F7FBFF","#9ECAE1","#4A6D8C"),te$abundance_limits,pretty(te$abundance_limits,n=3),"Mean relative abundance (log10 per-mille)")
    textmm(103,33,"Selected local bacterial MAGs",9.5,TRUE)
    grid.points(x=unit(105,"mm"),y=unit(27,"mm"),pch=23,size=unit(3,"mm"),gp=gpar(fill="#F2C14E",col="black"))
    textmm(109,27,"MAG122, MAG138, MAG156",9.2)
    textmm(4,18,"Phylum",9.5,TRUE)
    for(i in seq_along(e$phylum_breaks)) {
      x<-4+((i-1)%%3)*59;y<-12-floor((i-1)/3)*4.5
      grid.points(x=unit(x+1,"mm"),y=unit(y,"mm"),pch=16,size=unit(2.6,"mm"),gp=gpar(col=phylum_palette[e$phylum_breaks[i]]))
      textmm(x+4,y,unname(e$phylum_labels[e$phylum_breaks[i]]),8.8)
    }
  })
  legend_pdf(file.path(NATIVE,"Figure_S2_B_legend.pdf"),190,43,function(w,h) {
    x0<-c(4,68,132);textmm(x0[1],39,"Phylum",9.5,TRUE);textmm(x0[2],39,"Class",9.5,TRUE);textmm(x0[3],39,"K10944 (AmoA/PmoA)",9.5,TRUE)
    for(i in seq_along(te$archaeal_phylum_levels)) {
      y<-34-(i-1)*4.7;grid.points(x=unit(5.4,"mm"),y=unit(y,"mm"),pch=16,size=unit(2.6,"mm"),gp=gpar(col=phylum_palette[te$archaeal_phylum_levels[i]]));textmm(9,y,te$archaeal_phylum_levels[i],9)
    }
    for(i in seq_along(te$archaeal_class_levels)) {
      y<-34-(i-1)*4.7;grid.rect(x=unit(69.4,"mm"),y=unit(y,"mm"),width=unit(2.6,"mm"),height=unit(2.6,"mm"),gp=gpar(fill=te$archaeal_class_palette[te$archaeal_class_levels[i]],col=NA));textmm(73,y,te$archaeal_class_levels[i],9)
    }
    for(i in 1:2) {
      y<-34-(i-1)*4.7;grid.rect(x=unit(133.4,"mm"),y=unit(y,"mm"),width=unit(2.6,"mm"),height=unit(2.6,"mm"),gp=gpar(fill=te$pmo_palette[i],col="#B9C1C7",lwd=.4));textmm(137,y,c("not detected","present")[i],9)
    }
    colourbar(4,4,52,pal_seq,c(-3,1.5),c(-3,-1.5,0,1.5),"Relative abundance (log10 per-mille)")
    colourbar(105,4,52,te$completeness_colours,c(50,100),c(50,75,100),"Completeness (%)",scales::rescale(c(50,65,80,90,100),from=c(50,100)))
  })
  write_tsv(bind_rows(cells),file.path(AUDIT,"Figure_S2_B_plotted_cells.tsv"),quote="needed")
  stopifnot(identical(original_a_data,e$p_a$data),identical(original_b_data,e$p_b$data))
  write_tsv(te$mag_metrics %>% filter(Domain=="Bacteria") %>% transmute(
    MAG,phylum,focal=MAG %in% e$focal_mags,mean_relative_abundance=mean_rel_abundance,
    mean_per_mille,log10_mean_per_mille),file.path(AUDIT,"Figure_S2A_bacterial_tree_abundance.tsv"))
  write_tsv(te$ar_meta %>% select(MAG,phylum,class,Completeness,pmo_present),file.path(AUDIT,"Figure_S2B_archaeal_metadata.tsv"))
  write_tsv(as.data.frame(te$ar_abundance_matrix) %>% tibble::rownames_to_column("MAG"),file.path(AUDIT,"Figure_S2B_abundance_matrix.tsv"))
  write_json(list(bacteria=221,archaea=7,abundance_ring=te$abundance_matrix,archaeal_abundance=te$ar_abundance_matrix,
    branch_lengths_bacteria=te$bac_tree$edge.length,branch_lengths_archaea=te$ar_tree$edge.length,
    native_tree_coordinates_unchanged=TRUE,new_analysis=FALSE,ring_legend_matches_native_colour_scale=TRUE),file.path(AUDIT,"Figure_S2_QA.json"),pretty=TRUE,auto_unbox=TRUE,digits=16)
}

mode <- commandArgs(trailingOnly=TRUE)
if(!length(mode)) stop("Use script/figures.py for the approved selection")
stopifnot(all(mode %in% c("heatmaps", "figureS1", "figure6")))
for(m in mode) {message("Rendering ",m);get(m)()}
