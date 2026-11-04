#!/usr/bin/env Rscript
# Reuse only the original quantification expressions, not its superseded tests.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
ROOT<-normalizePath('.');OUT<-file.path(ROOT,'analysis/C03_C18_correction_2026-10-28/exclude_both')
Q<-file.path(OUT,'qa/curated_recalculation');dir.create(Q,recursive=TRUE,showWarnings=FALSE)
ids<-c('C01',sprintf('C%02d',4:17))
script<-file.path(ROOT,'release_2026-11-05/source_code/38_curated_module_validation.R')
lines<-readLines(script,warn=FALSE);end<-which(startsWith(lines,'# ---- Uniform multivariate tests'))
expr<-parse(script,keep.source=TRUE);refs<-attr(expr,'srcref')
e<-new.env(parent=globalenv())
e$source<-function(file,...)base::source(file,local=e)
for(i in seq_along(expr)){
 if(refs[[i]][1]>=end)break
 x<-expr[[i]]
 if(is.call(x)&&identical(x[[1]],as.name('<-'))){
   lhs<-as.character(x[[2]])
   if(length(lhs)==1L && lhs=='OUT'){e$OUT<-Q;next}
   if(length(lhs)==1L && lhs=='samples'){e$samples<-ids;next}
 }
 if(refs[[i]][1]<75 && is.call(x)&&identical(x[[1]],as.name('stopifnot')))next
 eval(x,e)
}
stopifnot(identical(e$meta$sample_id,ids))
fresh<-as.data.table(e$feature_tpm)
frozen<-fread(file.path(OUT,'inputs/curated_feature_sample_TPM.tsv'))[sample_id%in%ids]
key<-c('module','tier','category','marker','marker_label','sample_id','analysis_variant')
a<-merge(fresh[,c(key,'TPM'),with=FALSE],frozen[,c(key,'TPM'),with=FALSE],by=key,suffixes=c('_new','_old'))
stopifnot(nrow(a)==nrow(fresh),nrow(a)==nrow(frozen))
delta<-max(abs(a$TPM_new-a$TPM_old));stopifnot(delta<1e-7)
fwrite(a[,.(rows=.N,max_abs_TPM_difference=max(abs(TPM_new-TPM_old))),by=.(analysis_variant,module)],file.path(Q,'validation.tsv'),sep='\t')
groups<-c('module','tier','sample_id','site','site_display','depth_label','depth_display')
mt<-fresh[analysis_variant=='primary_marker_fractional',.(module_total_TPM=sum(TPM)),by=groups]
ct<-fresh[analysis_variant=='sensitivity_category_fractional',.(category_TPM=sum(TPM)),by=c(groups,'category')]
for(pair in list(list(mt,'curated_module_total_sample_TPM.tsv','module_total_TPM'),list(ct,'curated_category_sample_TPM.tsv','category_TPM'))){
 old<-fread(file.path(OUT,'inputs',pair[[2]]),na.strings=NULL)[sample_id%in%ids]
 cols<-setdiff(names(pair[[1]]),pair[[3]])
 z<-merge(pair[[1]],old,by=cols,suffixes=c('_new','_old'))
 stopifnot(nrow(z)==nrow(old),max(abs(z[[paste0(pair[[3]],'_new')]]-z[[paste0(pair[[3]],'_old')]]))<1e-7)
 fwrite(pair[[1]],file.path(Q,pair[[2]]),sep='\t')
}
write_json(list(status='complete',libraries=ids,max_abs_TPM_difference=delta,
 note='Curated markers use annotation-name filters and within-module fractional allocation; whole-KO counts use all valid KOs. They are distinct estimands.'),
 file.path(Q,'SUCCESS.json'),pretty=TRUE,auto_unbox=TRUE)
message('Original curated quantification reproduced for all retained samples; maximum difference ',delta)
