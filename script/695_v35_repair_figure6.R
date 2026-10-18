# Repair support-to-split mapping only; preserve V19 figure geometry and all source data.
suppressPackageStartupMessages({library(ape);library(jsonlite)})
PROJECT <- Sys.getenv("CHILE_PROJECT_ROOT", unset=getwd())
OUT <- file.path(PROJECT,'analysis/EM_V35_coauthor_candidate_2026-09-30')
Q <- file.path(OUT,'INTERNAL_QA/figure6_build')
dir.create(Q,recursive=TRUE,showWarnings=FALSE)
code <- paste(readLines(file.path(PROJECT,'script/470_render_geoderma_v03_figure6.R')),collapse='\n')
once <- function(old,new) {
 stopifnot(lengths(regmatches(code,gregexpr(old,code,fixed=TRUE)))==1L)
 code <<- sub(old,new,code,fixed=TRUE)
}
src <- file.path(PROJECT,'Environmental_Microbiome/V35/Supplementary_Data/AmoA_sequences')
once('PHYLO_DIR <- file.path(INPUT_ROOT, "phylogeny")',paste0('PHYLO_DIR <- "',src,'"'))
once('SOURCE_INPUT_DIR <- file.path(INPUT_ROOT, "source_data")',paste0('SOURCE_INPUT_DIR <- "',src,'"'))
once('  file.exists(CURRENT_S3), file.exists(TREE_FILE), file.exists(IQTREE_FILE),','  file.exists(TREE_FILE), file.exists(IQTREE_FILE),')
once('ape::root(outgroup = bacterial_tips, resolve.root = TRUE)','ape::root(outgroup = bacterial_tips, resolve.root = TRUE, edgelabel = TRUE)')
once('      label = round(ufboot)', '      label = round(ufboot),\n      x = x - ifelse(node %in% c(ape::getMRCA(tree_rooted, c("QUERY_MAG005", "QUERY_MAG070", "REF_Archaeal_AmoA_A0A060HNG6")), ape::getMRCA(tree_rooted, c("QUERY_MAG005", "REF_Archaeal_AmoA_A0A060HNG6"))), 0.055 * tree_xmax, 0),\n      y = y - ifelse(node == ape::getMRCA(tree_rooted, c("QUERY_MAG005", "QUERY_MAG070", "REF_Archaeal_AmoA_A0A060HNG6")), 0.35, 0)')
binding <- '
focal_node <- p_tree$data[p_tree$data$node == archaeal_mrca, ]
parent_node <- p_tree$data[p_tree$data$node == focal_node$parent, ]
stopifnot(nrow(focal_node)==1, nrow(parent_node)==1, archaeal_support=="100/100")
stem_x <- (focal_node$x + parent_node$x)/2
p_tree <- p_tree + annotate("segment", x=stem_x, xend=stem_x,
 y=focal_node$y, yend=focal_node$y+0.79, colour="grey25", linewidth=0.35)
'
once('# ---- MAG005/MAG070 context caveat',paste0(binding,'\n# ---- MAG005/MAG070 context caveat'))
# No legacy append-to-supplement or process quit; only the established standalone export.
code <- strsplit(code,'quit(save = "no", status = 0)',fixed=TRUE)[[1]][1]
Sys.setenv(CHILE_READABILITY_FIGURE5_ROOT=Q,CHILE_NATIVE_PANEL_DIR=file.path(Q,'native_panels'))
eval(parse(text=code))
desc <- function(t,n) {
 if(n<=Ntip(t)) return(t$tip.label[n])
 sort(unlist(lapply(t$edge[t$edge[,1]==n,2],function(k)desc(t,k))))
}
old <- ladderize(root(tree_unrooted,outgroup=bacterial_tips,resolve.root=TRUE,edgelabel=FALSE),right=FALSE)
stopifnot(identical(old$edge,tree_rooted$edge),identical(old$edge.length,tree_rooted$edge.length),identical(old$tip.label,tree_rooted$tip.label))
ledger <- lapply(seq_len(nrow(tree_rooted$edge)),function(i){
 child <- tree_rooted$edge[i,2]; side <- desc(tree_rooted,child)
 support <- if(child>Ntip(tree_rooted))tree_rooted$node.label[child-Ntip(tree_rooted)] else ''
 boot <- suppressWarnings(as.numeric(sub('^.*/','',support)))
 displayed <- if(child==archaeal_mrca) support_label else if(!is.na(boot)&&boot>=75)as.character(round(boot)) else ''
 row <- p_tree$data[p_tree$data$node==child,]
 list(parent=tree_rooted$edge[i,1],child=child,length=tree_rooted$edge.length[i],side_A=sort(side),side_B=sort(setdiff(tree_rooted$tip.label,side)),support=support,displayed=displayed,x=row$x,y=row$y)
})
write_json(list(ordinary_threshold=75,focal_node=archaeal_mrca,all_edges=ledger,
 unchanged_topology_and_lengths=TRUE,tree_refitted=FALSE,source_tree=TREE_FILE,
 label_only_offset='Corrected lower 86 and restored 78 moved left by 0.055*tree_xmax; lower86 also shifted down0.35 tree-y units to clear horizontal stem; all branch/tip coordinates unchanged.',
 source_renderer='script/470_render_geoderma_v03_figure6.R',v19_leader_preserved=TRUE),
 file.path(Q,'corrected_split_display_ledger.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat('Corrected Figure6 with identical rooted edge geometry; wrote split/display ledger.\n')
